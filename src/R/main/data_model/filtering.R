#'-----------------------------------------------------------------------------
#' ML Model Training Script
#'
#' This module contains code for training a variety of regression models on 
#' prepared genomics datasets (post feature-selection).
#' 
#' @author eczech
#'-----------------------------------------------------------------------------
options( java.parameters = "-Xmx4g" )
source('utils.R')
source('data_model/training_lib.R')
source('data_model/training_models.R')
source('data_model/training_filter.R')
source('data_model/training_viz.R')
source('data_model/filtering_lib.R')
source('~/repos/portfolio/functional/ml/R/trainer.R')
source('~/repos/portfolio/functional/ml/R/results.R')
source('~/repos/portfolio/functional/ml/R/caret_decorators.R')
source('data_model/filtering_models.R')
source('data_model/filtering_novelties.R')
source('data_model/filtering_settings.R')
source('data_pres/publication/visualizations.R')
lib('MASS')
lib('caret')
lib('doMC')
lib('iterators')
lib('caretEnsemble') # library(devtools); unload(inst("caretEnsemble")); install_local('/home/eczech/repos/misc/caretEnsemble');
lib('ROCR')
lib('plotly')
SEED <- 1024

## Choose dataset to use for modeling (must pick one of the following)
EnableCosmic()
#EnableCtd()

# TRAINER_DATA_DIRNAME <- 'filtering_data'
TRAINER_DATA_DIRNAME <- 'filtering_data.ga'
# TRAINER_DATA_DIRNAME <- 'filtering_data.ga.ds'

# RES_CACHE_DIRNAME <- 'filter_result_data'
RES_CACHE_DIRNAME <- 'filter_result_data.ga'
# RES_CACHE_DIRNAME <- 'filter_result_data.ga.ds'

PREPROC <- c('zv', 'center', 'scale')

RESULT_CACHE <- Cache(dir=file.path(CACHE_DIR, RES_CACHE_DIRNAME), project=RESPONSE_TYPE)
select <- dplyr::select

##### Data Loading #####

# Load training data for this response type
d.prep <- GetTrainingData(TRAIN_CACHE, RESPONSE_TYPE, RESPONSE_SELECTOR, min.mutations=3) %>% select(-tissue)

# Frequency of origin
# table(TransformOriginMostFrequent$convert(d.prep[,'origin'])) %>% prop.table

# Frequency of response
# table(DichotomizeOutcome(d.prep[,'response'], RESPONSE_THRESH))

# Load a dataset equivalent to the above but with missing response labels (ultimately to be predicted)
d.predict <- GetPredictionData(TRAIN_CACHE, RESPONSE_TYPE, PREDICTION_SELECTOR, names(d.prep)) %>% select(-tissue)

if (!all(names(d.prep) == names(d.predict)))
  stop('Column names in training data do not equal those in dataset to be predicted')

##### Data Partitioning #####

# # Downsampling (This did not go well -- all models predict at 50%)
# set.seed(SEED)
# ds <- downSample(d.prep %>% select(-response), DichotomizeOutcome(d.prep$response, threshold = RESPONSE_THRESH), list=T)
# mask <- d.prep[,'tumor_id'] %in% ds$x[,'tumor_id']
# d.tr <- list(X=(ds$x %>% select(-tumor_id)), y=d.prep[,'response'][mask], y.bin=ds$y)

set.seed(SEED)
idx.tr <- createDataPartition(d.prep[,'response'], p=.8)[[1]]

split.data <- function(data, idx, type, N){
  d <- data[idx,]; X <- d %>% select(-response, -tumor_id)
  y <- d[,'response']; y.bin <- DichotomizeOutcome(y, threshold = RESPONSE_THRESH)
  n <- length(y); n.pos <- sum(y.bin == 'pos')
  summary <- data.frame(n=n, pct.of.total=n/N, pos=n.pos, pos.pct=n.pos/n, type=type)
  list(X=X, y=y, y.bin=y.bin, summary=summary)
}
d.tr <- split.data(d.prep, idx.tr, 'training', nrow(d.prep))
d.ho <- split.data(d.prep, -idx.tr, 'holdout', nrow(d.prep))

##### Model Trainer #####

trainer <- Trainer(cache.dir=file.path(CACHE_DIR, TRAINER_DATA_DIRNAME), 
                   cache.project=RESPONSE_TYPE, seed=SEED)
trainer$generateFoldIndex(d.tr$y, CreateFoldIndex)
fold.data.gen <- GetFeatScoringFoldGen(PREPROC, RESPONSE_THRESH, feat.limit=5000, n.core=6)
trainer$generateFoldData(d.tr$X, d.tr$y, fold.data.gen, FilteringDataSummarizer())

##### Training Constants ######

# origin.trans <- TransformOriginSolidLiquid
origin.trans <- TransformOriginMostFrequent
n.feats <- c(c(0,1,2,3,4,5,6,7,8,9), seq(10, 50, by=10), 100, 200, 300)
#n.feats <- c(40)
#n.feats <- c(c(2))
origin.name <- 'wmostfreqorigin'

get.model.definition <- function(...){
  lapply(n.feats, function(x){ 
    list(model=GetModelForTransform(x, ...), max.feats=x)
  })
}

models.def <- list()
models.def$svm <- get.model.definition(
  'svm', n.core=3, origin.transform=origin.trans, origin.name=origin.name,
  method=GetFilteringSVMModel(), tuneLength=20, preProcess='zv', tol=1
)
models.def$enet <- get.model.definition(
  'enet', n.core=3, origin.transform=origin.trans, origin.name=origin.name,
  method='glmnet', tuneLength=20, preProcess='zv'
)
models.def$enetwt <- get.model.definition(
  'enetwt', n.core=3, origin.transform=origin.trans, origin.name=origin.name,
  weight.fun=function(y) ifelse(y == 'pos', 3, 1),
  method='glmnet', tuneLength=10, preProcess='zv'
)
models.def$xgb <- get.model.definition(
  'xgb', n.core=3, origin.transform=origin.trans, origin.name=origin.name,
  method='xgbTree', tuneLength=8, preProcess='zv'
)
models.def$gbm <- get.model.definition(
  'gbm', n.core=3, origin.transform=origin.trans, origin.name=origin.name,
  method='gbm', tuneLength=8, preProcess='zv', verbose=F
)
models.def$rf <- get.model.definition(
  'rf', n.core=3, origin.transform=origin.trans, origin.name=origin.name,
  method='rf', tuneLength=4, preProcess='zv'
)
# models.def$etree <- get.model.definition(
#   'etree', n.core=1, origin.transform=origin.trans, origin.name=origin.name,
#   method='extraTrees', tuneLength=4, preProcess='zv'
# )
# models.def$nnet <- get.model.definition(
#   'nnet', n.core=3, origin.transform=origin.trans, origin.name=origin.name,
#   method='nnet', tuneLength=6, preProcess='zv', trace=F
# )
models.def$knn <- get.model.definition(
  'knn', n.core=3, origin.transform=origin.trans, origin.name=origin.name,
  method='knn', tuneLength=5, preProcess='zv'
)
models.def$sda <- get.model.definition(
  'sda', n.core=3, origin.transform=origin.trans, origin.name=origin.name,
  method='sda', tuneGrid=expand.grid(diagonal = c(F, T), lambda = seq(.001, .999, length = 8)), 
  preProcess=c('zv', 'center', 'scale'), verbose=F
)
models.def$pls <- get.model.definition(
  'pls', n.core=3, origin.transform=origin.trans, origin.name=origin.name,
  method='pls', tuneLength=5,
  preProcess=c('zv', 'center', 'scale')
)
models.def$pam <- get.model.definition(
  'pam', n.core=3, origin.transform=origin.trans, origin.name=origin.name,
  method='pam', tuneLength=5,
  preProcess=c('zv', 'center', 'scale')
)
models.def$scrda <- get.model.definition(
  'scrda', n.core=3, origin.transform=origin.trans, origin.name=origin.name,
  method=GetSCRDAModel(5), 
  tuneLength=5,
  preProcess=c('zv', 'center', 'scale')
)
models.def$mars <- get.model.definition(
  'mars', n.core=8, origin.transform=origin.trans, origin.name=origin.name,
  method='earth', tuneLength=6, preProcess=c('zv', 'center', 'scale')
)

# models.def$hdrda <- get.model.definition(
#   'hdrda', n.core=3, origin.transform=origin.trans, origin.name=origin.name,
#   method=GetHDRDAModel(), tuneLength=5,
#   preProcess=c('zv', 'center', 'scale')
# )
# Add knn, lda, mda, pam, sda



# env <- new.env()
# load('/home/eczech/genomics_data_cache/filtering_data.ga/cosmic/model_mars_2_wmostfreqorigin.Rdata', envir=env)

# Takes too long
# models.def$ensavg <- get.model.definition(
#   'ensavg', n.core=6, origin.transform=origin.trans, origin.name=origin.name, allow.parallel=F,
#   method=m.ens.avg, tuneLength=1, preProcess='zv'
# )

ec <- T
get.trained.models <- function(models.def){
  lapply(names(models.def), function(m){
    model <- models.def[[m]]
    lapply(model, function(m.part){ 
      if (m == 'pls' && m.part$max.feats == 1){
        loginfo('Skipping PLS with max feats 1 since it seems to get stuck indefinitely')
        return(NULL)
      }
        
      #trainer$train(m.part$model, enable.cache=ec)
      tryCatch({
        trainer$train(m.part$model, enable.cache=ec)
      }, error=function(e){
        loginfo('Ignoring model %s for feature count %s due to error', m, m.part$max.feats)
        return(NULL)
      })
    }) %>% setNames(sapply(model, function(m) m$model$name)) %>% .[!sapply(., is.null)]
  }) %>% setNames(names(models.def))
}
models <- get.trained.models(models.def)

# Validate predictions from models in CV and make sure no NA's were allowed 
# in prediction matrixes attached to caret training results
ValidatePredictions(models)

# Create performance summary over all models (will be used for visualization was well as ensemble selection)
registerDoMC(3)
cv.res.all <- GetAggregateFilterCVRes(models)

# GenerateErrorProfileVis(cv.res.all)

##### Ensemble Selection #####

model.perf <- cv.res.all %>% select(model, fold, kappa, n.feats)
model.cor <- model.perf %>% filter(n.feats == 20) %>% group_by(n.feats) %>% do({
  d <- .
  cor.mat <- dcast(d , fold ~ model, value.var='kappa') %>% select(-fold)
  cor.mat <- suppressWarnings(cor(cor.mat) %>% melt() %>% setNames(c('m1', 'm2', 'cor')))
  cor.mat <- cor.mat %>% filter(m1 != m2) %>% mutate(cor=ifelse(is.na(cor), 1, cor))
  cor.mat <- cor.mat %>% mutate(m1=as.character(m1), m2=as.character(m2))
  cor.mat
}) %>% ungroup

model.cor %>% group_by(m1, m2) %>% 
  summarise(cor=mean(cor)) %>% ungroup %>%
  mutate(cor=cut(cor, breaks=3)) %>%
  ggplot(aes(x=m1, y=m2, fill=cor)) + geom_tile() + 
  scale_fill_manual(values = c("green", "yellow", "red"), guide=guide_legend(title='Error Correlation')) + 
  theme_bw() +
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  xlab('Model') + ylab('Model')

# cv.res.all %>% filter(n.feats==20) %>% group_by(model) %>% 
#   summarise(value=mean(kappa)) %>%
#   arrange(desc(value))

##### Ensemble Training #####


sub.models <- models[GetEnsembleModelNames()]
if (any(is.na(names(sub.models))))
  stop('Failed to extract models for ensembling')
ens.models.def <- list()
ens.models.def$ensglm <- GetEnsModelDefinition(
  'ensglm', sub.models, method='glm', origin.transform=origin.trans, 
  origin.name=origin.name, n.core=3
)
ens.models.def$ensavg <- GetEnsModelDefinition(
  'ensavg', sub.models, method=GetEnsembleAveragingModel(), origin.transform=origin.trans, 
  origin.name=origin.name, n.core=3
)
ens.models <- get.trained.models(ens.models.def)
for (m in names(ens.models)) models[[m]] <- ens.models[[m]]

# Again, validate that predictions from all models (including ensembles now) are valid
ValidatePredictions(models)

# Again, create performance summary over all models (including ensembles)
registerDoMC(3)
cv.res.all <- GetAggregateFilterCVRes(models)

##### Performance #####

# Create visualizations for genomics conference
perf.plots <- cv.res.all %>% 
  filter((model != 'pam' | n.feats > 5) & (model != 'mars')) %>% 
  GeneratePerfProfileVis(RESPONSE_TYPE)

cv.res.model <- 'ensglm'
get.cv.res <- function(cv.res.model){
  SummarizeTrainingResults(
    models[[cv.res.model]], T, fold.summary=ResSummaryFun('roc'), model.summary=ResSummaryFun('roc')
  ) 
}
cv.res <- get.cv.res(cv.res.model)

# cv.perf.key <- paste('cv_model_perf', origin.name, sep='_')
# RESULT_CACHE$store(cv.perf.key, cv.res)
# cv.res <- RESULT_CACHE$load('cv_model_perf_norigin')

PlotFoldConfusion(cv.res)
PlotFoldMetric(cv.res, 'acc')
PlotFoldMetric(cv.res, 'cacc', order.by.score=F) 
PlotFoldMetric(cv.res, 'kappa', order.by.score=F)
PlotFoldMetric(cv.res, 'nir')
PlotFoldMetric(cv.res, 'auc')
PlotFoldMetric(cv.res, 'mcp')
PlotFoldMetric(cv.res, 'ptp')
PlotFoldMetric(cv.res, 'cacc') +
  ggtitle(sprintf('%s %s Accuracy Over Baseline', toupper(RESPONSE_TYPE), cv.res.model)) + 
  ggsave(sprintf('~/repos/musc_genomics/src/R/main/data_pres/images/filtering/%s/cacc_%s.png', RESPONSE_TYPE, cv.res.model))

# Look at selected feature frequencies
cv.top.feats <- ExtractTopFeatures(models)
# cv.top.feats %>% filter(fold == 2) %>% group_by(n.feats) %>% tally
RESULT_CACHE$store('cv_top_feats', cv.top.feats)

# cv.top.n.feat <- cv.top.feats %>% filter(n.feats == 20) %>% group_by(feature) %>% tally %>%
#   setNames(c('feature', 'frequency')) %>% arrange(desc(frequency)) %>% data.frame
# print(cv.top.n.feat, row.names=F)
# GenerateTopFeatVis(cv.top.n.feat) + ylim(-8, 5) + coord_flip()
  

##### Top Feature Subset Accuracies #####

top.feat.ct <- 20
top.model.names <- GetModelNamesByFeatCount(models, top.feat.ct, origin.name)
top.models <- GetModelsByFeatCount(models, top.feat.ct, origin.name)

top.model.res <- SummarizeTrainingResults(
  top.models, T, fold.summary=ResSummaryFun('roc'), model.summary=ResSummaryFun('roc')
)
PlotFoldMetric(top.model.res, 'cacc', extract.model.names=T) +
  ylab('Accuracy Over Baseline') +
  ggtitle(sprintf('%s Accuracy Over Baseline w/ %s Features', toupper(RESPONSE_TYPE), top.feat.ct)) + 
  ggsave(sprintf('~/repos/musc_genomics/src/R/main/data_pres/images/filtering/%s/cacc_%s_feats.png', RESPONSE_TYPE, top.feat.ct))

PlotFoldMetric(top.model.res, 'kappa', extract.model.names=T) + 
  ylab('Kappa') +
  ggtitle(sprintf('%s Kappa w/ %s Features', toupper(RESPONSE_TYPE), top.feat.ct)) + 
  ggsave(sprintf('~/repos/musc_genomics/src/R/main/data_pres/images/filtering/%s/kappa_%s_feats.png', RESPONSE_TYPE, top.feat.ct))

PlotFoldMetric(top.model.res, 'acc', extract.model.names=T) +
  ylab('Raw Accuracy') +
  ggtitle(sprintf('%s Raw Accuracy w/ %s Features', toupper(RESPONSE_TYPE), top.feat.ct)) + 
  ggsave(sprintf('~/repos/musc_genomics/src/R/main/data_pres/images/filtering/%s/acc_%s_feats.png', RESPONSE_TYPE, top.feat.ct))

##### Holdout Fit #####

ho.fold.data.gen <- GetFeatScoringFoldGen(PREPROC, RESPONSE_THRESH, feat.limit=5000, n.core=1)
ho.fit.key <- paste('holdout_fit', origin.name, sep='_')
ho.dat.key <- paste('holdout_data', origin.name, sep='_')
# trainer$getCache()$invalidate(ho.fit.key)
# trainer$getCache()$invalidate(ho.dat.key)

# Remove models that cause fitting issues in holdout
ho.models.def <- models.def[!names(models.def) %in% c('pls')]

# Fit models on holdout data
ho.fit <- RunHoldoutFit(trainer, ho.models.def, d.tr, d.ho, ho.fit.key, ho.dat.key, ho.fold.data.gen, n.feat=GetHoldoutFeatCt())
ho.models <- RestackHoldoutFit(ho.fit)

### Ensemble Holdout ###

ho.ens.models.def <- GetEnsembleModelDefFromSubModelFit(
  ho.models, GetEnsembleModelNames(), origin.trans, origin.name, n.feat=GetHoldoutFeatCt()
)

ho.ens.fit.key <- paste('holdout_fit_ens', origin.name, sep='_')
# trainer$getCache()$invalidate(ho.ens.fit.key)
ho.ens.fit <- RunHoldoutFit(trainer, ho.ens.models.def, d.tr, d.ho, ho.ens.fit.key, ho.dat.key, ho.fold.data.gen)
ho.ens.models <- RestackHoldoutFit(ho.ens.fit)
for (m in names(ho.ens.models)) ho.models[[m]] <- ho.ens.models[[m]]

ValidatePredictions(ho.models)

##### Holdout Results #####

ho.top.feats <- ExtractTopFeatures(ho.models)
# ho.top.feats %>% group_by(n.feats) %>% tally
RESULT_CACHE$store('holdout_top_feats', ho.top.feats)

ho.fit.top <- GetModelsByFeatCount(ho.models, GetOptimalFeatCt(), origin.name)
ho.res <- SummarizeTrainingResults(list(ho.fit.top), T, model.summary=ResSummaryFun('lift'))
# ho.perf.key <- paste('ho_model_perf', origin.name, sep='_')
# RESULT_CACHE$store(ho.perf.key, ho.res)
# ho.res <- RESULT_CACHE$load('ho_model_perf_wmostfreqorigin')

PlotHoldOutConfusion(ho.res)
PlotHoldOutLift(ho.res)
PlotHoldOutMetric(ho.res, 'auc') 
PlotHoldOutMetric(ho.res, 'acc') 
PlotHoldOutMetric(ho.res, 'cacc') 
PlotHoldOutMetric(ho.res, 'kappa') 
PlotHoldOutMetric(ho.res, 'mcp') 


#### Holdout Gain/Sensitivity Analysis

ho.lift.model <- GetOptimalModel()
fit <- ho.models[[ho.lift.model]][[sprintf('%s.%s.wmostfreqorigin', ho.lift.model, GetOptimalFeatCt())]]

# Create gain/sensitivity data for holdout set
n.pos <- sum(fit$y.test == 'pos')
n.pred <- length(fit$y.pred$prob)
rocr.pred <- prediction(fit$y.pred$prob, fit$y.test)
rocr.perf <- performance(rocr.pred, 'tpr', 'rpp') # gain chart
d.ho.lift <- data.frame(x=rocr.perf@x.values[[1]], y=rocr.perf@y.values[[1]]) %>%
  mutate(x=x * n.pred, y=y*n.pos) %>%
  mutate(sensitivity=y/x) %>% 
  filter(x > 0)

# Save results above in result cache
ho.lift.key <- sprintf('ho_lift_data_%s', ho.lift.model)
RESULT_CACHE$store(ho.lift.key, d.ho.lift)
# d.ho.lift <- RESULT_CACHE$load(ho.lift.key)

PlotHoldoutGainFinal(d.ho.lift) +
  ggsave(sprintf('~/repos/musc_genomics/src/R/main/data_pres/images/filtering/%s/ho_lift_%s.png', RESPONSE_TYPE, ho.lift.model))

PlotHoldoutSensitivityFinal(d.ho.lift) +
  ggsave(sprintf('~/repos/musc_genomics/src/R/main/data_pres/images/filtering/%s/ho_sens_%s.png', RESPONSE_TYPE, ho.lift.model))


##### Calibration #####

cal.data <- get.cv.res('rf')$predictions
#cal.data <- ho.res$predictions
cal.data %>% group_by(model) %>% do({
  d <- .
  p <- seq(0, 1, by=.1)
  d %>% mutate(bin=cut(y.pred.prob, breaks=p, include.lowest=T)) %>%
    group_by(bin) %>% summarise(pct.pos=sum(y.test == 'pos')/n(), pct.neg=sum(y.test == 'neg')/n(), n=n())
}) %>% ggplot(aes(x=bin, y=pct.pos)) + geom_bar(stat='identity') + facet_wrap(~model)

##### True Out-of-sample Predictions #####

uk.fold.data.gen <- GetFeatScoringFoldGen(PREPROC, RESPONSE_THRESH, feat.limit=5000, n.core=1)
uk.fit.key <- paste('unknown_fit', origin.name, sep='_')
uk.dat.key <- paste('unknown_data', origin.name, sep='_')
# trainer$getCache()$invalidate(uk.fit.key)
# trainer$getCache()$invalidate(uk.dat.key)

# Remove models that cause fitting issues in holdout
uk.models.def <- models.def[!names(models.def) %in% c('pls')]

# Create filter to remove models with irrelevant feature subset sizes
uk.n.feat <- GetUnlabeledFeatCt()

# Append cross-validation data to holdout data as the training data set for unknown samples
d.tr.uk <- list(X=rbind(d.tr$X, d.ho$X), y=c(d.tr$y, d.ho$y))
# Create the true "holdout" dataset on unlabeled data (with arbitrary response values since 
# performance measures taken against those response will be ignored)
d.ho.uk <- list(X=d.predict[,names(d.tr.uk$X)], y=rep(0, nrow(d.predict)))

# Refit all models on all labeled samples, and use those fit models to make predictions
# for samples with no known label
uk.fit <- RunHoldoutFit(trainer, uk.models.def, d.tr.uk, d.ho.uk, uk.fit.key, uk.dat.key, uk.fold.data.gen, n.feat=uk.n.feat)
uk.models <- RestackHoldoutFit(uk.fit)

## Unknown Label Prediction Ensembles ##

# Create ensemble model definitions using the base models already fit above
uk.ens.models.def <- GetEnsembleModelDefFromSubModelFit(uk.models, GetEnsembleModelNames(), origin.trans, origin.name, n.feat=uk.n.feat)

uk.ens.fit.key <- paste('unknown_fit_ens', origin.name, sep='_')
# trainer$getCache()$invalidate(uk.ens.fit.key)

# Train ensemble models for unlabeled samples
uk.ens.fit <- RunHoldoutFit(trainer, uk.ens.models.def, d.tr.uk, d.ho.uk, uk.ens.fit.key, uk.dat.key, uk.fold.data.gen)
uk.ens.models <- RestackHoldoutFit(uk.ens.fit)
for (m in names(uk.ens.models)) uk.models[[m]] <- uk.ens.models[[m]]

## Save top features ##

# Extract the features selected when training over the entire known, labeled dataset for this source
uk.top.feats <- ExtractTopFeatures(uk.models)
# uk.top.feats %>% group_by(n.feats) %>% tally
RESULT_CACHE$store('unknown_top_feats', uk.top.feats) # Store the top feature data frame for later use

##### Identify Out-of-Sample Novelties #####

# Generate a vector containing the names of the top features chosen over all labeled samples
uk.subset.feats <- uk.top.feats %>% filter(n.feats == GetOptimalFeatCt()) %>% .$feature

# Create novelty/anomaly training dataset using top features 
X.tr.nov <- d.tr.uk$X[,uk.subset.feats] %>% select(-origin) 

# Create dataset for which anomaly predictions will be made (this includes all samples
# with no known label)
X.ho.nov <- d.ho.uk$X[,uk.subset.feats] %>% select(-origin)

# Collect novelty info based on tumor origin
uk.ho.origin.ind <- NoveltyAttributesByOrigin(d.predict, d.tr.uk, d.ho.uk, origin.trans)

# Collect novelty info using two-class SVM with anomaly simulation
nov.two.class <- NoveltyAttributesTwoClassSVM(X.tr.nov, X.ho.nov, SEED, N.sim.factor=3, n.core=3)
# summary(nov.two.class$novelty.model$resample$Accuracy)
# PlotUnlabeledHoldout(X.ho.nov, nov.two.class$novelty.predictions) +
#   ggsave(sprintf('~/repos/musc_genomics/src/R/main/data_pres/images/filtering/%s/uknown_sample_parallel_coordinates.png', RESPONSE_TYPE))

# Collect novelty info using one-class SVM
nov.one.class <- NoveltyAttributesOneClassSVM(X.tr.nov, X.ho.nov, SEED)

# Collect novelty info using euclidean distance
nov.dist <- NoveltyAttributesByDist(X.tr.nov, X.ho.nov)

# Aggregate all anamoly/novelty information into a single data frame
uk.meta <- uk.ho.origin.ind %>%
  mutate(is.outlier.one.class=!nov.one.class$novelty.predictions) %>%
  mutate(is.outlier.two.class=percent_rank(1-nov.two.class$novelty.predictions)) %>%
  mutate(mean.distance.rank=percent_rank(nov.dist))
# uk.meta %>% ggplot(aes(y=is.outlier.two.class, x=mean.distance.rank)) + geom_point() + geom_smooth()

##### Prediction Export #####

# Aggregate top features across all modeling rounds
all.top.feats <- GetTopFeatures(RESULT_CACHE, GetOptimalFeatCt())

# Extract the models fit across all labeled samples that also pertain only
# to the feature subset size of interest
pred.models <- GetModelsByFeatCount(uk.models, GetOptimalFeatCt(), origin.name)

# For each model selected above, extract the respective predictions on unlabeled
# samples and aggregate all of those predictions (from each model into a single data frame)
pred.exp <- lapply(names(pred.models), function(m){
  model <- pred.models[[m]]
  if (!length(model$fit$top.feats) == GetOptimalFeatCt()) stop(sprintf(
      'Final prediction aggregation for model "%s" was expecting %s top features but found %s', 
      m, final.feat.ct(), length(model$fit$top.feats)
  ))
  data.frame(as.character(model$y.pred$class), model$y.pred$prob) %>%
    setNames(c(sprintf('%s.class', m), sprintf('%s.prob', m)))
})
pred.exp <- do.call(cbind, pred.exp)
if (any(is.na(pred.exp))) stop('Found NA values in final predictions')

# Calculate the most frequently predicted class for each unlabeled sample
avg.pred.class <- pred.exp %>% select(ends_with('class')) %>% 
  apply(1, function(x) names(sort(table(x),decreasing=TRUE))[1])

# Calculate the average predicted probability of sensitivity for each unlabeled sample
avg.pred.prob <- pred.exp %>% select(ends_with('prob')) %>% apply(1, mean)

# Attach the most frequently predicted classes and average probabilities to the prediction
# data frame assembled so far 
pred.exp <- pred.exp %>% mutate(overall.class=avg.pred.class, overall.prob=avg.pred.prob)

# Plot correlation in predicted probabilities across all models
# pred.exp %>% select(ends_with('prob')) %>% cor %>% melt %>%
#   mutate(value=cut(value, breaks=c(-1, .7, .8, .9, 1), labels=c('low', '>.7', '>.8', '>.9'))) %>%
#   ggplot(aes(x=Var1, y=Var2, fill=value)) + geom_tile() +
#   theme(axis.text.x = element_text(angle=25, hjust=1))

# Attaches the metadata associated with each unlabeled sample (including the likelihood of those
# samples being significantly different from known training data) to the prediction information from
# each model and saves the resulting information in a CSV file
save.predictions <- function(metadata, predictions, type, sort.model='ensglm'){
  exp.file <- sprintf(
    '~/repos/musc_genomics/data/predictions/%s/unknown_sample_%s_predictions.csv', 
    RESPONSE_TYPE, type
  )
  exp.data <- cbind(metadata, predictions %>% select(ends_with(type))) %>% 
    arrange_(sprintf('desc(%s)', paste(sort.model, type, sep='.'))) %>% 
    select(tumor_id, origin, everything())
  write.csv(exp.data, file=exp.file, row.names=F)
  exp.data
}

# Create and save separate data frames for both the predicted classes and probabilities of each unlabeled sample
pred.exp.prob <- save.predictions(uk.meta, pred.exp, 'prob')
pred.exp.class <- save.predictions(uk.meta, pred.exp, 'class')

# pred.exp.final[,c(names(pred.exp.final)[1:7], rev(names(pred.exp.final))[1:3])] %>% head(10)
# pred.exp.final %>% ggplot(aes(x=origin, y=ensglm.prob)) + geom_boxplot() + coord_flip()

##### Rectified Predictions #####

get.rectified.accuracy <- function(
  y.pred, y.test, 
  p.lo.seq=seq(.025, .5, by=.025), 
  p.hi.seq=seq(.5, .975, by=.025)){
  
  foreach(p.lo=p.lo.seq, .combine=rbind) %:%
    foreach(p.hi=p.hi.seq, .combine=rbind) %do%{
      p.rec <- factor(sapply(y.pred, function(p){
        if (p <= p.lo) 'neg'
        else if (p >= p.hi) 'pos'
        else 'na'
      }))
      n.na <- sum(p.rec == 'na')
      n.pos <- sum(p.rec == 'pos')
      n.neg <- sum(p.rec == 'neg')
      n.pos.correct <- sum(p.rec == 'pos' & y.test == 'pos')
      n.neg.correct <- sum(p.rec == 'neg' & y.test == 'neg')
      n.miss <- length(p.rec) - n.pos.correct - n.neg.correct - n.na
      acc.pos <- n.pos.correct/n.pos
      acc.neg <- n.neg.correct/n.neg
      acc.all <- sum(p.rec == 'pos' & y.test == 'pos') + sum(p.rec == 'neg' & y.test == 'neg')
      acc.all <- acc.all / (n.pos + n.neg)
      data.frame(
        n=length(p.rec), n.na=n.na, n.pos=n.pos, n.neg=n.neg, 
        n.pos.correct=n.pos.correct, n.neg.correct=n.neg.correct, n.miss=n.miss,
        acc.pos=acc.pos, acc.neg=acc.neg, acc.all=acc.all,
        p.lo=p.lo, p.hi=p.hi
      )
    }
}

data.frame(pred=fit$y.pred$prob, true=fit$y.test) %>%
  mutate(bin=cut(pred, breaks=seq(0, 1, by=.25))) %>%
  group_by(bin) %>% summarise(acc=sum(true == 'pos')/n(), n=n())

ho.acc.rec <- get.rectified.accuracy(
  fit$y.pred$prob, fit$y.test,
  p.lo.seq=seq(.08, .4, by=.01), p.hi.seq=seq(.55, .92, by=.01)
)

ho.acc.rec %>% filter(abs(p.lo + p.hi - 1) < .001) %>%
  mutate(window=.5 - p.lo) %>%
  mutate(window=sprintf('.5 +/- %s', str_pad(window, 4, pad='0', side='right'))) %>%
  mutate(n.miss=n - n.pos.correct - n.neg.correct) %>%
  mutate(pct.na=n.na/n) %>%
  select(pct.na, acc.pos, acc.neg, window) %>%
  melt(id.vars=c('window')) %>%
  ggplot(aes(x=factor(window), y=value, color=variable)) + geom_point() +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))


##### Simple Filter Models #####

# source('~/repos/portfolio/functional/ml/R/parallel.R')
# registerCores(3, log.file = '/tmp/musc_genomics_training.log')
registerDoMC(3)

RunXGBFilter <- function(
  X, y, max.feats=50, tune.length=10,
  outer.method='repeatedcv', outer.number=10, outer.repeats=1,
  inner.method='cv', inner.number=8){
  seed <- 1024
  GetXGBFilter <- function(limit){
    model.args <- list(
      method='xgbTree',
      metric='Accuracy', 
      preProcess=c('zv'),
      tuneLength=tune.length,
      trControl=trainControl(
        method=inner.method, number=inner.number, classProbs=T, 
        returnData=T, savePredictions='final',
        allowParallel=F, verboseIter=T
      )
    )
    GetLimitFilter(seed, limit, model.args)
  }
  set.seed(seed)
  sbfctrl <- sbfControl(
    functions=GetXGBFilter(max.feats), 
    method=outer.method, number=outer.number, repeats=outer.repeats, 
    saveDetails=T, verbose=T, allowParallel=T)
  set.seed(seed)
  sbf(X, y, sbfControl = sbfctrl)
}


set.seed(SEED)
d.samp <- d.prep %>% sample_frac(1)
# d.samp <- d.prep
X <- d.samp %>% select(-tumor_id, -response, -origin) %>%
  mutate(origin=ifelse(toupper(origin) == 'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 1, 0))
y <- DichotomizeOutcome(d.samp[,'response'], threshold=RESPONSE_THRESH)
# table(x=X[,'origin'], y=y)
# fisher.test(X[,'origin'], y)$p.value

# Benchmarking runs
by.max.feats <- lapply(seq(5, 50, by=5), function(i){
  cat(sprintf('Running model for max feat %s\n', i))
  t.start <- proc.time()
  res.xgb <- RunXGBFilter(
    X, y, outer.method='cv', outer.number=10, outer.repeats=1,
    inner.method='cv', inner.number=5, max.feats=i, tune.length=8)
  time.diff <- proc.time() - t.start  
  list(max.feat=i, time=time.diff, res=res.xgb)
})
res.xgb.w.o <- res.xgb

data.frame(r1=res.xgb$resample$Accuracy, r2=res.xgb.w.o$resample$Accuracy) %>%
  melt(id.vars=NULL) %>%
  ggplot(aes(x=variable, y=value)) + geom_boxplot()

# 298.063 seconds for sample frac .3 (123 rows)
# 377 seconds for sample frac .3
# 311 seconds for sample frac 1, w/ domc(3)
# 629 seconds for sample frac 1, w/ domc(3)

# Actual run
# res.xgb <- RunXGBFilter(X, y)

