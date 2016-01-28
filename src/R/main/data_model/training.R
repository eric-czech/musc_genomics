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
lib('MASS')
lib('caret')
lib('doMC')
lib('iterators')

#registerDoMC(8)
SEED <- 1024
model.cache <- Cache(dir=CACHE_DIR, project='training_models')
data.cache <- Cache(dir=CACHE_DIR, project='training_data')

RESPONSE_TYPE <- 'cosmic' # This will include ctd2 as well at some point
RESPONSE_SELECTOR <- function(d){ 
  d %>% filter(!is.na(ic_50)) %>% rename(response=ic_50) %>% select(-auc)
}
select <- dplyr::select

d.prep <- GetTrainingData(data.cache, RESPONSE_TYPE, RESPONSE_SELECTOR)

set.seed(SEED)
d.prep.tr <- d.prep %>% sample_frac(.8, replace = F)
d.prep.ho <- d.prep %>% filter(!tumor_id %in% d.prep.tr$tumor_id)
trctrl <- function(...) trainControl(method = "cv", number = 3, savePredictions='final', ...) # 3
fsctrl <- function(...) sbfControl(method = "cv", number = 5, ...) # 5
#preproc <- c('zv', 'center', 'scale')
preproc <- c('center', 'scale')

set.seed(SEED)
#d.samp <- d.prep.tr %>% sample_frac(.2)
#X <- d.samp[,sample(c(c.ge, c.cn, c.mu), replace=F, size = 100)]; y <- d.samp[,'response']
X <- d.prep.tr[,c(c.ge, c.cn, c.mu)]; y <- d.prep.tr[,'response']

loginfo('Creating hyperparameter estimates')
X.preproc <- predict(preProcess(X, method=preproc), X)
glmnet.lambda <- GetGlmnetLambda(X.preproc)
svm.sigma <- GetSvmSigma(X.preproc)
trctrl <- function(...) trainControl(method = "cv", number = 5, savePredictions='final', ...)


set.seed(SEED)
GetFoldFile <- function(fold) sprintf('%s_%s', RESPONSE_TYPE, fold)
cv.res <- foreach(fold=createFolds(y, k = 10, list = T, returnTrain = F), i=icount()) %do% {
  loginfo('Running CV fold %s', i)
  #browser()
  
  loader <- function(){
    X.train.all <- X[-fold,]; y.train <- y[-fold]
    X.test <- X[fold,]; y.test <- y[fold]
    
    # Apply feature selector
    loginfo('Running feature selection')
    registerDoMC(8)
    X.dim <- dim(X.train.all)
    X.train.sml <- ApplyFeatureFilter(X.train.all, y.train, c.numeric, c.binary,
                                      numeric.score.p=.0001, binary.score.p=.15)
    
    # Removing zero-variance features
    # X.train.sml <- ApplyZeroVarianceFilter(X.train.sml)
    # X.train.all <- ApplyZeroVarianceFilter(X.train.all)
    
    # Apply preprocessing to feature subset
    loginfo('Running preprocessing')
    pp.all <- preProcess(X.train.all, method=preproc)
    pp.sml <- preProcess(X.train.sml, method=preproc)
    X.train.all <- predict(pp.all, X.train.all)
    X.train.sml <- predict(pp.sml, X.train.sml)
    X.test.all  <- predict(pp.all, X.test)
    X.test.sml  <- predict(pp.sml, X.test[,names(X.train.sml)])
    
    list(
      preproc=list(pp.all=pp.all, pp.sml=pp.sml),
      X.train.sml=X.train.sml, X.train.all=X.train.all,
      X.test.sml=X.test.sml, X.test.all=X.test.all,
      y.train=y.train, y.test=y.test
    )
  }
  d <- FetchFromDisk(GetFoldFile(i), loader, dir=GetCachePath('cv_fold_data'))

  
  gfl <- function(d, t) length(GetFeatures(d, t))
  logdebug('Dimension reduction summary after preprocessing:')
  logdebug('Mutation features   : %s --> %s', gfl(d$X.train.all, 'mu'), gfl(d$X.train.sml, 'mu'))
  logdebug('Expression features : %s --> %s', gfl(d$X.train.all, 'ge'), gfl(d$X.train.sml, 'ge'))
  logdebug('Copy Number features: %s --> %s', gfl(d$X.train.all, 'cn'), gfl(d$X.train.sml, 'cn'))
  
  loginfo('Running PLS models')
  set.seed(SEED); registerDoMC(5)
  m.pls.1 <- train(
    d$X.train.sml, d$y.train, method='pls', 
    tuneGrid=data.frame(ncomp=c(1,2,3)),
    trControl = trctrl()
  )
  
  loginfo('Running SVM models')
  set.seed(SEED); registerDoMC(3)
  m.svm.1 <- train(
    d$X.train.sml, d$y.train, method='svmRadial', 
    tuneLength=300, search='random',
    trControl = trctrl(verboseIter=T)
  )
  
  set.seed(SEED); registerDoMC(3)
  m.svm.2 <- train(
    d$X.train.sml, d$y.train, method='svmRadial', 
    tuneGrid = expand.grid(.sigma = as.numeric(svm.sigma[-2]), .C = 2^(1:6)),
    trControl = trctrl(verboseIter=T)
  )
  
  models <- list(
    svm.grid=m.svm.1,
    svm.est=m.svm.2,
    pls=m.pls.1
  )
  list(models=models, data=d)
}

foreach(r=cv.res, .combine=rbind)%do%{
  r$models$svm.grid$bestTune
}
# Best tuning params for svm over large tune length:
# sigma C
# 6  0.0002670770 8
# 5  0.0002895063 4
# 51 0.0003733021 4
# 52 0.0003325179 4
# 53 0.0002968672 4
# 54 0.0003208822 4
# 61 0.0003000269 8
# 55 0.0003545485 4
# 56 0.0002722678 4
# 62 0.0003385101 8
cv.scores <- foreach(r=cv.res, .combine=rbind)%do%{
  foreach(m=names(r$models), .combine=rbind)%do%{
    y.pred <- predict(r$models[[m]], r$data$X.test.sml[,names(r$data$X.train.sml)])
    y.true <- r$data$y.test
    r2   <- GetRegressorScore(y.true, y.pred, 'rsquared')
    rmse <- GetRegressorScore(y.true, y.pred, 'rmse')
    data.frame(model=m, rsquared=r2, rmse=rmse)
  }
}
cv.scores %>% ggplot(aes(x=model, y=rmse)) + geom_point() + geom_boxplot()
cv.scores %>% ggplot(aes(x=model, y=rsquared)) + geom_point() + geom_boxplot()
