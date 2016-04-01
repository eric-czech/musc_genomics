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
source('data_model/filtering_models.R')
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
#EnableCosmic()
EnableCtd()

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
d.prep <- GetTrainingData(TRAIN_CACHE, RESPONSE_TYPE, RESPONSE_SELECTOR, min.mutations=3)

# Load a dataset equivalent to the above but with missing response labels (ultimately to be predicted)
d.predict <- GetPredictionData(TRAIN_CACHE, RESPONSE_TYPE, PREDICTION_SELECTOR, names(d.prep))

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

GetPrepFun <- function(origin.transform=NULL){
  function(X){
    if (!'origin' %in% names(X)) X
    else {
      if (is.null(origin.transform))
        stop('Origin found in feature set but transform given was null')
      lvls <- as.character(0:(origin.transform$max.val))
      X <- X %>% mutate(origin=factor(origin.transform$convert(origin), levels=lvls))
      as.data.frame(model.matrix(~. - 1, X))
    }
  }
}

BIAS_FEATS <- c('bias1', 'bias2')
GetRandVars <- function(n){
  set.seed(SEED)
  data.frame(bias1=rnorm(n), bias2=rnorm(n))
}

GetFilterModel <- function(name, max.feats, n.core=3, 
                           origin.transform=NULL, ...){
  prep.fun <- GetPrepFun(origin.transform)
  pred.fun <- function(fit, d, i){ 
    X.test <- d$X.test
    if (length(fit$top.feats) == length(BIAS_FEATS) && all(fit$top.feats == BIAS_FEATS)){
      X.test <- GetRandVars(nrow(X.test))
    }
    X <- X.test[, fit$top.feats, drop=F] %>% prep.fun
    list(
      prob=predict(fit, X, type='prob')[,1], 
      class=predict(fit, X, type='raw')
    )
  }
  list(
    name=name, predict=pred.fun, test=bin.test,
    train=function(d, idx, i){
      registerDoMC(n.core)
      # d$feat.scores contains scores for all 36k features
      feats <- d$feat.scores
      if (is.null(origin.transform))
        feats <- feats %>% filter(feature != 'origin')
      
      # Temporarily filter to only GE and MU features
      # feats <- feats %>% filter(str_detect(feature, '^mu\\.') | str_detect(feature, '^ge\\.'))
      
      X.train <- d$X.train
      if (max.feats == 0){
        top.feats <- BIAS_FEATS
        X.train <- GetRandVars(nrow(X.train))
      } else {
        top.feats <- feats %>% arrange(score) %>% head(max.feats) %>% .$feature 
      }
      X <- prep.fun(X.train[,top.feats, drop=F])
      y <- d$y.train.bin
      
      fit <- train(X, y, ...)
      fit$top.feats <- top.feats
      trim_model(fit) # Remove massive, embedded '.Environment' attributes
    }
  )
}

GetModelForTransform <- function(max.feats, model.name, n.core=3, k=5, allow.parallel=T, 
                                 origin.transform=NULL, origin.name=NULL, ...){
  if (is.null(origin.name))
    origin.name <- ifelse(is.null(origin.transform), 'norigin', 'worigin')
  name <- sprintf('%s.%s.%s', model.name, max.feats, origin.name)
  GetFilterModel(
    name, max.feats, n.core=n.core,
    origin.transform=origin.transform,
    metric='Accuracy', 
    ...,
    trControl=trainControl(
      method='cv', number=k, classProbs=T, 
      returnData=F, savePredictions='final',
      allowParallel=allow.parallel, verboseIter=F
    )
  )
}

GetEnsembleModelForTransform <- function(max.feats, sub.models, model.name, method='glm', 
                                         origin.transform=NULL, origin.name=NULL, ...){
  if (is.null(origin.name))
    origin.name <- ifelse(is.null(origin.transform), 'norigin', 'worigin')
  name <- sprintf('%s.%s.%s', model.name, max.feats, origin.name)
  prep.fun <- GetPrepFun(origin.transform)
  pred.fun <- function(fit, d, i){ 
    X.test <- d$X.test
    if (length(fit$top.feats) == length(BIAS_FEATS) && all(fit$top.feats == BIAS_FEATS)){
      X.test <- GetRandVars(nrow(X.test))
    }
    X <- X.test[, fit$top.feats, drop=F] %>% prep.fun
    tryCatch({
      list(
        prob=predict(fit, X, type='prob'), 
        class=predict(fit, X, type='raw')
      )
    }, error=function(e) browser())
  }
  fit.prep.fun <- function(fit, model.results){
    fit$top.feats <- model.results[[1]]$top.feats
    fit
  }
  res <- c(
    list(name=name, predict=pred.fun, test=bin.test),
    GetFitEnsembleTrain(
      sub.models, fit.prep.fun=fit.prep.fun, 
      method=method, metric='Accuracy', trControl=trainControl(
        method='none', number=1, 
        classProbs=T, savePredictions='final', returnData=T
      )
    ) 
  )
  res
}


# Static settings


# origin.trans <- TransformOriginSolidLiquid
origin.trans <- TransformOriginMostFrequent
n.feats <- c(c(0,1,2,3,4,5,6,7,8,9), seq(10, 50, by=10), 100, 200, 300)
#n.feats <- c(0, 1, 3, 10, 20, 50, 100)
#n.feats <- c(0,1,2,3)
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
  method='svmRadial', tuneLength=20, preProcess='zv'
)
models.def$enet <- get.model.definition(
  'enet', n.core=3, origin.transform=origin.trans, origin.name=origin.name,
  method='glmnet', tuneLength=20, preProcess='zv'
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
  method=GetSCRDAModel(5), tuneLength=5,
  preProcess=c('zv', 'center', 'scale')
)

# models.def$hdrda <- get.model.definition(
#   'hdrda', n.core=3, origin.transform=origin.trans, origin.name=origin.name,
#   method=GetHDRDAModel(), tuneLength=5,
#   preProcess=c('zv', 'center', 'scale')
# )
# Add knn, lda, mda, pam, sda

# models.def$mars <- get.model.definition(
#   'mars', n.core=8, origin.transform=origin.trans, origin.name=origin.name,
#   method='earth', tuneLength=6, preProcess=c('zv', 'center', 'scale')
# )

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
        logdebug('Skipping PLS with max feats 1 since it seems to get stuck indefinitely')
        return(NULL)
      }
        
      tryCatch({
        trainer$train(m.part$model, enable.cache=ec)  
      }, error=function(e){
        logerror('Ignoring model %s for feature count %s due to error', m, m.part$max.feats)
        return(NULL)
      })
    }) %>% setNames(sapply(model, function(m) m$model$name)) %>% .[!sapply(., is.null)]
  }) %>% setNames(names(models.def))
}
models <- get.trained.models(models.def)

# Create performance summary over all models (will be used for visualization was well as ensemble selection)
registerDoMC(3)
cv.res.all <- GetAggregateFilterCVRes(models)

##### Ensemble Selection #####

model.perf <- cv.res.all %>% select(model, fold, kappa, n.feats)
model.cor <- model.perf %>% group_by(n.feats) %>% do({
  d <- .
  cor.mat <- dcast(d , fold ~ model, value.var='kappa') %>% select(-fold)
  cor.mat <- suppressWarnings(cor(cor.mat) %>% melt() %>% setNames(c('m1', 'm2', 'cor')))
  cor.mat <- cor.mat %>% filter(m1 != m2) %>% mutate(cor=ifelse(is.na(cor), 1, cor))
  cor.mat <- cor.mat %>% mutate(m1=as.character(m1), m2=as.character(m2))
  cor.mat
})

model.cor %>% group_by(m1, m2) %>% 
  summarise(cor=mean(cor)) %>% 
  mutate(cor=cut(cor, breaks=2)) %>%
  ggplot(aes(x=m1, y=m2, fill=cor)) + geom_tile() + 
  scale_fill_manual(values = c("green", "red"))

cv.res.all %>% filter(n.feats==20) %>% group_by(model) %>% summarise(kappa=mean(kappa)) %>%
  arrange(desc(kappa))

ens.model.names <- c('rf', 'knn', 'pam', 'svm', 'xgb')

##### Ensemble Training #####

get.ens.model.definition <- function(name, sub.models, ...){
  lapply(n.feats, function(x){ 
    feat.ct.models <- GetModelsByFeatCount(sub.models, x, origin.name)
    list(model=GetEnsembleModelForTransform(x, feat.ct.models, name, ...), max.feats=x)
  })
}

sub.models <- models[ens.model.names]
if (any(is.na(names(sub.models))))
  stop('Failed to extract models for ensembling')
ens.models.def <- list()
ens.models.def$ensglm <- get.ens.model.definition(
  'ensglm', sub.models, method='glm', origin.transform=origin.trans, origin.name=origin.name
)
ens.models.def$ensavg <- get.ens.model.definition(
  'ensavg', sub.models, method=GetEnsembleAveragingModel(), origin.transform=origin.trans, origin.name=origin.name
)
ens.models <- get.trained.models(ens.models.def)
for (m in names(ens.models)) models[[m]] <- ens.models[[m]]
  

##### Performance #####

# Create visualizations for genomics conference
cv.res.all %>% filter(model != 'pam' | n.feats > 5) %>% GeneratePerfProfileVis(RESPONSE_TYPE)

    
cv.res.model <- 'ensglm'
cv.res <- SummarizeTrainingResults(
  models[[cv.res.model]], T, fold.summary=ResSummaryFun('roc'), model.summary=ResSummaryFun('roc')
)

cv.perf.key <- paste('cv_model_perf', origin.name, sep='_')
RESULT_CACHE$store(cv.perf.key, cv.res)
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
top.model <- 'svm.10.wmostfreqorigin'
lapply(models$svm[[top.model]], function(m) m$fit$top.feats) %>% unlist %>% table %>%
  as.data.frame %>% arrange(desc(Freq)) %>% setNames(c('feature', 'frequency'))

##### Top Feature Subset Accuracies #####

top.feat.ct <- 20
top.models <- GetModelsByFeatCount(models, top.feat.ct, origin.name)

top.model.res <- SummarizeTrainingResults(
  top.models, T, fold.summary=ResSummaryFun('roc'), model.summary=ResSummaryFun('roc')
)
PlotFoldMetric(top.model.res, 'kappa') + 
  ggtitle(sprintf('%s Accuracy Over Baseline w/ %s Features', toupper(RESPONSE_TYPE), top.feat.ct)) + 
  ggsave(sprintf('~/repos/musc_genomics/src/R/main/data_pres/images/filtering/%s/cacc_%s_feats.png', RESPONSE_TYPE, top.feat.ct))

##### Holdout Fit #####

ho.fit.key <- paste('holdout_fit', origin.name, sep='_')
ho.dat.key <- paste('holdout_data', origin.name, sep='_')
# trainer$getCache()$invalidate(ho.fit.key)
# trainer$getCache()$invalidate(ho.dat.key)
ho.fold.data.gen <- GetFeatScoringFoldGen(PREPROC, RESPONSE_THRESH, feat.limit=5000, n.core=1)
ho.fit <- trainer$getCache()$load(ho.fit.key, function(){
  ho.models <<- list()
  for (m in names(models.def)){
    for (m.fun in models.def[[m]]){
      ho.models[[m.fun$model$name]] <<- m.fun$model
    }
  }
  trainer$holdout(ho.models, d.tr$X, d.tr$y, d.ho$X, d.ho$y, ho.fold.data.gen, ho.dat.key) 
})

# Look at selected features
ho.fit.top <- ho.fit[top.model.names]
ho.fit.top[[1]]$fit$top.feats

ho.res <- SummarizeTrainingResults(list(ho.fit.top), T, model.summary=ResSummaryFun('lift'))
ho.perf.key <- paste('ho_model_perf', origin.name, sep='_')
RESULT_CACHE$store(ho.perf.key, ho.res)
# ho.res <- RESULT_CACHE$load('ho_model_perf_wmostfreqorigin')

PlotHoldOutConfusion(ho.res)
PlotHoldOutLift(ho.res)
PlotHoldOutMetric(ho.res, 'auc') 
PlotHoldOutMetric(ho.res, 'acc') 
PlotHoldOutMetric(ho.res, 'cacc') 
PlotHoldOutMetric(ho.res, 'kappa') 
PlotHoldOutMetric(ho.res, 'mcp') 

##### Calibration #####

cal.data <- cv.res$predictions
#cal.data <- ho.res$predictions
cal.data %>% group_by(model) %>% do({
  d <- .
  p <- seq(0, 1, by=.1)
  d %>% mutate(bin=cut(y.pred.prob, breaks=p, include.lowest=T)) %>%
    group_by(bin) %>% summarise(pct.pos=sum(y.test == 'pos')/n(), pct.neg=sum(y.test == 'neg')/n(), n=n())
}) %>% ggplot(aes(x=bin, y=pct.pos)) + geom_bar(stat='identity') + facet_wrap(~model)

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
      acc.pos <- sum(p.rec == 'pos' & y.test == 'pos')/n.pos
      acc.neg <- sum(p.rec == 'neg' & y.test == 'neg')/n.neg
      acc.all <- sum(p.rec == 'pos' & y.test == 'pos') + sum(p.rec == 'neg' & y.test == 'neg')
      acc.all <- acc.all / (n.pos + n.neg)
      data.frame(
        n=length(p.rec), n.na=n.na, n.pos=n.pos, n.neg=n.neg, 
        acc.pos=acc.pos, acc.neg=acc.neg, acc.all=acc.all,
        p.lo=p.lo, p.hi=p.hi
      )
    }
}

fit <- ho.fit$xgb.10.worigin
ho.acc.rec <- get.rectified.accuracy(
  fit$y.pred$prob, fit$y.test,
  p.lo.seq=.3, p.hi.seq=.7
)

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

