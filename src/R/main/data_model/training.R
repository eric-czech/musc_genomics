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
source('data_model/training_viz.R')
source('~/repos/portfolio/functional/ml/R/trainer.R')
source('~/repos/portfolio/functional/ml/R/results.R')
lib('MASS')
lib('caret')
lib('doMC')
lib('iterators')
lib('caretEnsemble') # library(devtools); unload(inst("caretEnsemble")); install_local('/home/eczech/repos/misc/caretEnsemble');
lib('ROCR')
lib('plotly')
SEED <- 1024

RESPONSE_TYPE <- 'cosmic' 
RESPONSE_SELECTOR <- function(d){d %>% filter(!is.na(ic_50)) %>% rename(response=ic_50) %>% select(-auc)}
RESPONSE_THRESH <- -1 

# RESPONSE_TYPE <- 'ctd' 
# RESPONSE_SELECTOR <- function(d){d %>% filter(!is.na(auc)) %>% rename(response=auc) %>% select(-ic_50)}
# RESPONSE_THRESH <- -1

PREPROC <- c('zv', 'center', 'scale')

RESULT_CACHE <- Cache(dir=file.path(CACHE_DIR, 'result_data'), project=RESPONSE_TYPE)
select <- dplyr::select

d.prep <- GetTrainingData(TRAIN_CACHE, RESPONSE_TYPE, RESPONSE_SELECTOR, min.mutations=3)

## Create training, holdout, and calibration datasets
set.seed(SEED)
idx.tr <- createDataPartition(d.prep[,'response'], p=(1/2))[[1]]
idx.ho <- createDataPartition(d.prep[,'response'][-idx.tr], p=(2/3))[[1]]

split.data <- function(data, idx, type, N){
  d <- data[idx,]; X <- d %>% select(-response, -tumor_id)
  y <- d[,'response']; y.bin <- DichotomizeOutcome(y, threshold = RESPONSE_THRESH)
  n <- length(y); n.pos <- sum(y.bin == 'pos')
  summary <- data.frame(n=n, pct.of.total=n/N, pos=n.pos, pos.pct=n.pos/n, type=type)
  list(X=X, y=y, y.bin=y.bin, summary=summary)
}
d.tr <- split.data(d.prep, idx.tr, 'training', nrow(d.prep))
d.ho <- split.data(d.prep[-idx.tr,], idx.ho, 'holdout', nrow(d.prep))
d.cb <- split.data(d.prep[-idx.tr,], -idx.ho, 'calibration', nrow(d.prep))

d.summary <- do.call('rbind', lapply(list(d.tr, d.cb, d.ho), function(x) x$summary))
if (sum(d.summary$n) != nrow(d.prep)) stop('Sum of all observations in split data does not equal total')
# d.summary  # Show the size and class relevance within each data partition

## Response Analysis (determining classification cutoffs)
PlotResponseDist(RESPONSE_TYPE, d.prep[,'response'], RESPONSE_THRESH)
RESULT_CACHE$store('response_data', d.prep[,'response'])

# X <- X[1:25,c(GetFeatures(X, 'cn')[1:500], GetFeatures(X, 'ge')[1:500], GetFeatures(X, 'mu')[1:300])]
# y <- y[1:25]

##### Define Model Trainer #####

trainer <- Trainer(cache.dir=file.path(CACHE_DIR, 'training_data'), 
                   cache.project=RESPONSE_TYPE, seed=SEED)
trainer$generateFoldIndex(d.tr$y, CreateFoldIndex)
fold.data.gen <- GetFoldDataGenerator(PREPROC, RESPONSE_THRESH, F, n.core=8, sml.num.p=.0001, 
                                      lrg.num.p=.01, sml.bin.p=.1, lrg.bin.p=.15, pca.thresh=.99)
trainer$generateFoldData(d.tr$X, d.tr$y, fold.data.gen, GetDataSummarizer())


##### Classifiication Models #####

bin.sml.models <- list()
bin.pca.models <- list() 
bin.lrg.models <- list() 

    
# ShowBestTune(bin.models$gbm)

# Complete (takes about 2 hr 30 minutes to run all of the following)
ec <- T
bin.sml.models$svm.radial.sml <- trainer$train(bin.model.svm.radial.sml, enable.cache=ec)
bin.sml.models$svm.linear.sml <- trainer$train(bin.model.svm.linear.sml, enable.cache=ec)
bin.sml.models$pls <- trainer$train(bin.model.pls.sml, enable.cache=F)
bin.sml.models$pam <- trainer$train(bin.model.pam.sml, enable.cache=ec)
bin.sml.models$knn <- trainer$train(bin.model.knn.sml, enable.cache=ec)
bin.sml.models$rf <- trainer$train(bin.model.rf.sml, enable.cache=ec)
bin.sml.models$lasso <- trainer$train(bin.model.lasso.sml, enable.cache=ec)
bin.sml.models$ridge <- trainer$train(bin.model.ridge.sml, enable.cache=ec)
bin.sml.models$enet <- trainer$train(bin.model.enet.sml, enable.cache=ec)
bin.sml.models$gbm <- trainer$train(bin.model.gbm.sml, enable.cache=ec)


source('data_model/training_models.R')
bin.sml.models$rda <- trainer$train(bin.model.rda.sml, enable.cache=F)

bin.sml.models$scrda <- trainer$train(bin.model.scrda.sml, enable.cache=F)
bin.lrg.models$scrda <- trainer$train(bin.model.scrda.lrg, enable.cache=F)

#bin.sml.models$mars <- trainer$train(bin.model.mars.sml, enable.cache=ec)

ec <- T
bin.pca.models$svm.radial.pca <- trainer$train(bin.model.svm.radial.pca, enable.cache=ec)
bin.pca.models$pls <- trainer$train(bin.model.pls.pca, enable.cache=F)
bin.pca.models$pam <- trainer$train(bin.model.pam.pca, enable.cache=ec)
bin.pca.models$knn <- trainer$train(bin.model.knn.pca, enable.cache=ec)
bin.pca.models$rf <- trainer$train(bin.model.rf.pca, enable.cache=ec)
bin.pca.models$lasso <- trainer$train(bin.model.lasso.pca, enable.cache=ec)
bin.pca.models$ridge <- trainer$train(bin.model.ridge.pca, enable.cache=ec)
bin.pca.models$enet <- trainer$train(bin.model.enet.pca, enable.cache=ec)
bin.pca.models$gbm <- trainer$train(bin.model.gbm.pca, enable.cache=ec)

# Under Construction

bin.models$rda <- trainer$train(bin.model.rda, enable.cache=F)
bin.sml.models$et <- trainer$train(bin.model.et.sml, enable.cache=ec)
bin.pca.models$et <- trainer$train(bin.model.et.pca, enable.cache=ec)

# bin.models$nb <- trainer$train(bin.model.nb)
# bin.models$rf <- trainer$train(bin.model.rf)
# bin.models$gbm <- trainer$train(bin.model.gbm)

scale <- function(x) (x - mean(x)) / sd(x)
X.ge <- d.tr$X[,GetFeatures(d.tr$X, 'ge')] %>% mutate_each(funs(scale)) 
X.cn <- d.tr$X[,GetFeatures(d.tr$X, 'cn')] %>% mutate_each(funs(scale)) 

set.seed(SEED)
registerDoMC(3)
folds <- createMultiFolds(d.tr$y.bin, k=10, times=3)
fit.scrda.ge <- train(
  X.cn, d.tr$y.bin, method=GetSCRDAModel(3), preProcess='zv',
  metric=bin.tgt.metric, tuneLength=8,
  trControl=trainControl(index=folds, savePredictions=T, classProbs=T, verboseIter=T)
)
var.imp <- varImp(fit.scrda.ge, alpha=fit.scrda.ge$bestTune$alpha, delta=fit.scrda.ge$bestTune$delta, scale=F)
table(var.imp$importance$Overall)

ShowBestTune(bin.sml.models$scrda)
params <- bin.sml.models$scrda[[3]]$fit$bestTune
var.imp <- varImp(bin.sml.models$scrda[[3]]$fit, alpha=params$alpha, delta=params$delta, scale=F)
table(var.imp$importance)

##### Classification Ensembles #####

ens.models <- list(
  pam=function(i) bin.sml.models$pam[[i]]$fit,
  #pls=function(i) bin.sml.models$pls[[i]]$fit,
  gbm=function(i) bin.sml.models$gbm[[i]]$fit,
  rf=function(i) bin.sml.models$rf[[i]]$fit,
  knn=function(i) bin.sml.models$knn[[i]]$fit,
  glmnet=function(i) bin.sml.models$enet[[i]]$fit,
  svmRadial=function(i) bin.sml.models$svm.radial.sml[[i]]$fit,
  svmLinear=function(i) bin.sml.models$svm.linear.sml[[i]]$fit
)
bin.model.ens1 <- GetEnsembleModel(ens.models, 'bin.ens1', bin.test, 
                                   bin.predict.ens.sml, method='glm',
                                   metric=bin.tgt.metric, trControl=trainControl(method='cv', savePredictions = 'final'))
bin.sml.models$bin.model.ens1 <- trainer$train(bin.model.ens1, enable.cache=F)

##### Classification Hold Out #####

sml.models <- list(
  bin.model.scrda.sml, bin.model.rf.sml, bin.model.svm.radial.sml, bin.model.pam.sml, bin.model.pls.sml, 
  bin.model.knn.sml, bin.model.enet.sml, bin.model.lasso.sml, bin.model.ridge.sml,
  bin.model.svm.linear.sml, bin.model.gbm.sml
)
if (any(sapply(sml.models, is.null))) stop('Some models are null')
# trainer$getCache()$invalidate('holdout_fit')
ho.sml.fit <- trainer$getCache()$load('holdout_fit', function(){
  trainer$holdout(sml.models, d.tr$X, d.tr$y, d.ho$X, d.ho$y, fold.data.gen, 'holdout_data') 
})

pca.models <- list(
  bin.model.rf.pca, bin.model.svm.radial.pca, bin.model.pam.pca, bin.model.pls.pca, 
  bin.model.knn.pca, bin.model.enet.pca, bin.model.lasso.pca, bin.model.ridge.pca
)
pca.models <- list(
  bin.model.svm.radial.pca
)
# trainer$getCache()$invalidate('holdout_pca_fit')
ho.pca.fit <- trainer$getCache()$load('holdout_pca_fit', function(){
  trainer$holdout(pca.models, d.tr$X, d.tr$y, d.ho$X, d.ho$y, fold.data.gen, 'holdout_data') 
})

# trainer$getCache()$invalidate('calibration_fit')
cb.fit <- trainer$getCache()$load('calibration_fit', function(){
  trainer$holdout(models, d.tr$X, d.tr$y, d.cb$X, d.cb$y, fold.data.gen, 'calibration_data') 
})


ens.models.ho <- lapply(ho.sml.fit, function(m) {function(i) m$fit}) %>% 
  setNames(sapply(ho.sml.fit, function(m) m$model))
bin.model.ens1.ho <- GetEnsembleModel(ens.models.ho, 'bin.ens1.ho', bin.test,  
                                   bin.predict.ens.sml, method='glm',
                                   metric=bin.tgt.metric, trControl=trainControl(method='cv', savePredictions = 'final'))
ens.ho.fit <- trainer$holdout(list(t1=bin.model.ens1.ho), d.tr$X, d.tr$y, d.ho$X, d.ho$y, fold.data.gen, 'holdout_data')



# ho.data <- bs.data.gen(X.ho, y.ho.bin, X.ho, y.ho.bin)

##### Classification Results ##### 

## CV Results

cv.res <- SummarizeTrainingResults(bin.lrg.models, T, fold.summary=ResSummaryFun('pr'), model.summary=ResSummaryFun('pr'))
RESULT_CACHE$store('cv_model_perf', cv.res)

# ROC curves per-fold
PlotPerFoldROC(cv.res)
PlotPerFoldPR(cv.res)

# ROC curves across folds
PlotAllFoldROC(cv.res) %>% ggplotly() %>% layout(showlegend = T) %>% plot.ly
PlotAllFoldPR(cv.res)

# AUC ranges by model
PlotFoldMetric(cv.res, 'auc')
PlotFoldMetric(cv.res, 'bacc')
PlotFoldMetric(cv.res, 'acc')
PlotFoldMetric(cv.res, 'kappa')
PlotFoldMetric(cv.res, 'sens')
PlotFoldMetric(cv.res, 'spec')
PlotFoldMetric(cv.res, 'mcp')
PlotFoldMetric(cv.res, 'acc_margin_0.1')
PlotFoldMetricByMargin(cv.res, 'acc')


## Holdout results

# Calibration checks
cal.data <- cv.res$predictions
#cal.data <- ho.res$predictions
cal.data %>% group_by(model) %>% do({
  d <- .
  p <- seq(0, 1, by=.1)
  d %>% mutate(bin=cut(y.pred.prob, breaks=p, include.lowest=T)) %>%
    group_by(bin) %>% summarise(pct.pos=sum(y.test == 'pos')/n(), pct.neg=sum(y.test == 'neg')/n(), n=n())
}) %>% ggplot(aes(x=bin, y=pct.pos)) + geom_bar(stat='identity') + facet_wrap(~model)


ho.res <- SummarizeTrainingResults(list(ho.sml.fit), T, fold.summary=NULL, model.summary=ResSummaryFun('pr'))
RESULT_CACHE$store('ho_model_perf', ho.res)


cb.res <- SummarizeTrainingResults(list(cb.sml.fit), T, fold.summary=NULL, model.summary=ResSummaryFun('lift'))
RESULT_CACHE$store('cb_model_perf', cb.res)

p.res <- ho.res
PlotHoldOutMetric(p.res, 'auc') 
PlotHoldOutMetric(p.res, 'acc')
PlotHoldOutMetric(p.res, 'kappa')
PlotHoldOutMetric(p.res, 'sens')
PlotHoldOutMetric(p.res, 'spec')
PlotHoldOutROC(p.res)
PlotHoldOutPR(p.res)
PlotHoldOutLift(p.res)

cb.acc <- foreach(fold=cb.fit, .combine=rbind) %do% {
  p.max <- .9; p.min <- .1; p.inc <- .1
  p.range <- range(fold$y.pred$prob)
  foreach(lo=seq(p.range[1], p.range[2], length=10), .combine=rbind)%do%{
    foreach(hi=seq(lo, p.range[2], length=10), .combine=rbind)%do%{
      if (lo >= hi) return(NULL)
      y.range <- sapply(fold$y.pred$prob, function(x) if (x <= lo) 'lo' else if (x >= hi) 'hi' else 'na')
      stats <- data.frame(range=y.range, y.pred=fold$y.pred$class, y.test=fold$y.test) %>%
        group_by(range) %>% do({
          d <- data.frame(.)
          data.frame(n=nrow(d), acc=sum(d$y.pred == d$y.test)/nrow(d))
        }) %>% ungroup %>%
        mutate(model=fold$model, lo=lo, hi=hi, range = as.character(range))
      stats
    }
  }
}
n.grps <- cb.acc %>% group_by(model, range, lo, hi) %>% summarise(ct=n()) %>% .$ct %>% unique %>% length
if (n.grps != 1) stop('Found unexpected group duplicates')

cb.score <- cb.acc %>% group_by(model, lo, hi) %>% do({
  d <- .
  
  acc <- data.frame(d)$acc %>% setNames(d$range)
  ct <- data.frame(d)$n %>% setNames(d$range)
  N <- sum(ct)
  is.valid <- ct['lo'] >= 5 && ct['hi'] >= 5
  
  #value <- ifelse(is.valid, .8 * as.numeric(acc['hi']) + .2 * as.numeric(acc['lo']), NA)
  value <- ifelse(is.valid, acc['hi'] * acc['lo'], NA)
  
  data.frame(value=value)
})
cb.score %>% na.omit %>% ggplot(aes(x=factor(lo), y=factor(hi), fill=value)) + 
  geom_tile() + facet_wrap(~model, scales='free')

##### Regression Models #####

reg.models <- list()

reg.models$svm.radial <- trainer$train(reg.model.svm.radial, enable.cache=F)

reg.models$pls <- trainer$train(reg.model.pls, enable.cache=F)

model.glmnet <- list(
  train=function(d){
    registerDoMC(5)
    glmnet.lambda <- GetGlmnetLambda(d$X.train.sml, d$y.train)
    train(
      d$X.train.sml, d$y.train, method='glmnet',
      tuneGrid = expand.grid(.alpha = seq(0, .001, length.out = 2), .lambda=glmnet.lambda),
      trControl = trctrl(verboseIter=T)
    )
  }, predict=predict.test.data
)
models$glmnet <- trainer.rs$train(model.glmnet)


foreach(fold=models$glmnet) %do% fold$fit$bestTune

dt <- foreach(fold=fit.glmnet, .combine=rbind)%do%{
  fold$fit$results
}
svm.grid.fit <- trainer$train(svm.grid)

cv.scores <- foreach(m=names(models), .combine=rbind) %do% {
  foreach(fold=models[[m]], .combine=rbind) %do% {
    #browser()
    y.pred <- fold$y.pred
    y.test <- fold$y.test
    r2   <- GetRegressorScore(y.test, y.pred, 'rsquared')
    rmse <- GetRegressorScore(y.test, y.pred, 'rmse')
    data.frame(model=m, rsquared=r2, rmse=rmse)
  }
}

cv.scores %>% ggplot(aes(x=model, y=rmse)) + geom_point() + geom_boxplot()
cv.scores %>% ggplot(aes(x=model, y=rsquared)) + geom_point() + geom_boxplot()




