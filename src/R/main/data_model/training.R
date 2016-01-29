#'-----------------------------------------------------------------------------
#' ML Model Training Script
#'
#' This module contains code for training a variety of regression models on 
#' prepared genomics datasets (post feature-selection).
#' 
#' @author eczech
#'-----------------------------------------------------------------------------
source('utils.R')
source('data_model/training_lib.R')
source('~/repos/portfolio/functional/ml/R/trainer.R')
lib('MASS')
lib('caret')
lib('doMC')
lib('iterators')

#registerDoMC(8)
SEED <- 1024

RESPONSE_TYPE <- 'cosmic' # This will include ctd2 as well at some point
RESPONSE_SELECTOR <- function(d){ 
  d %>% filter(!is.na(ic_50)) %>% rename(response=ic_50) %>% select(-auc)
}
select <- dplyr::select

d.prep <- GetTrainingData(TRAIN_CACHE, RESPONSE_TYPE, RESPONSE_SELECTOR, min.mutations=3)

set.seed(SEED)
d.prep.tr <- d.prep %>% sample_frac(.8, replace = F)
d.prep.ho <- d.prep %>% filter(!tumor_id %in% d.prep.tr$tumor_id)
#preproc <- c('zv', 'center', 'scale')
preproc <- c('center', 'scale')

set.seed(SEED)
#d.samp <- d.prep.tr %>% sample_frac(.2)
#X <- d.samp[,sample(c(c.ge, c.cn, c.mu), replace=F, size = 100)]; y <- d.samp[,'response']

X <- d.prep.tr %>% select(-response, -tumor_id); y <- d.prep.tr[,'response']

trctrl <- function(...) trainControl(method = "cv", number = 5, savePredictions='final', ...)
data.generator <- GetFoldDataGenerator(preproc)

set.seed(SEED)

trainer <- Trainer(cache.dir=file.path(CACHE_DIR, 'training_data'), cache.project=RESPONSE_TYPE, seed=SEED)
trainer$generateFolds(function() createFolds(y, k = 10, list = T, returnTrain = F))
trainer$generateFoldData(X, y, data.generator, GetDataSummarizer())
predict.test.data <- function(fit, d){ predict(fit, d$X.test.sml[,names(d$X.train.sml)]) }


models <- list()

model.svm.radial <- list(
  train=function(d){
    registerDoMC(3)
    train(
      d$X.train.sml, d$y.train, method='svmRadial', 
      tuneLength=25, search='grid',
      trControl = trctrl(verboseIter=T)
    )
  }, predict=predict.test.data
)
models$svm.radial <- trainer$train(model.svm.radial)

train(d$X.train.sml, d$y.train, method = "svmRadial", tuneLength = 25, trControl = trctrl())

model.pls <- list(
  train=function(d){
    registerDoMC(3)
    train(
      d$X.train.sml, d$y.train, method='pls', 
      tuneGrid=data.frame(ncomp=c(1,2,3)),
      trControl = trctrl(verboseIter=T)
    )
  }, predict=predict.test.data
)
models$pls <- trainer$train(model.pls)

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
models$glmnet <- trainer$train(model.glmnet)

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


