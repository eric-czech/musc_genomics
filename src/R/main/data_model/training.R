#'-----------------------------------------------------------------------------
#' ML Model Training Script
#'
#' This module contains code for training a variety of regression models on 
#' prepared genomics datasets (post feature-selection).
#' 
#' @author eczech
#'-----------------------------------------------------------------------------
source('utils.R')
source('data_prep/prep.R')
source('data_model/training_lib.R')
lib('MASS')
lib('caret')
lib('doMC')
lib('iterators')


#registerDoMC(8)
SEED <- 1024
MODEL_DIR <- '~/genomics_data_cache/models'
RESPONSE_TYPE <- 'cosmic-only' # This will include ctd2 as well at some point
select <- dplyr::select



d <- GetPreparedData()

fields <- d$fields
c.ge <- fields$gene_expression
c.cn <- fields$copy_number
c.mu <- fields$mutations
d <- d$data

Prep <- function(d, c.numeric, c.binary, min.mutations=5){
  d %>% filter(!is.na(ic_50)) %>%                    # For now we're using COSMIC data only
    rename(response=ic_50) %>% select(-auc) %>%      # Select response field
    Filter(function(x)!all(is.na(x)), .) %>%         # Remove NA-only columns
    RemoveNA(row.threshold=.1) %>%                    # Remove rows with large % NA's
    RemoveRareMutations(c.binary, min.mutations) %>% # Remove cols for rare mutations
    mutate(response=ScaleVector(response))           # Scale response
  # nearZeroVar would make sense here as well but no features would be removed in past tests
}
d.prep <- Prep(d, c(c.ge, c.cn), c.mu)
c.mu <- GetFeatures(d.prep, 'mu')
c.cn <- GetFeatures(d.prep, 'cn')
c.ge <- GetFeatures(d.prep, 'ge')
c.numeric <- c(c.cn, c.ge)
c.binary <- c.mu
if (any(is.na(d.prep)))
  stop('Dataset contains unexpected NA values')



# d.prep <- data.frame(
#   tumor_id=c('t1', 't2', 't3', 't4', 't5', 't6'),
#   ic_50=c(-2.3, -1.2, -.3, .9, 1.5, 2.8),
#   auc=rnorm(6),
#   'cn.A1BG'=c(2.1, 1.3, .4, -.7, -1.3, -3),
#   'cn.A1BG.AS1'=c(-1, -2, -5, -3, -4, -3),
#   'ge.A1BG'=rnorm(6),
#   'ge.A1BG.AS1'=c(-6, -3, -1, 1, 3, 3.5),
#   'mu.XPO4_L499F'=c(0, 0, 0, 1, 1, 1),
#   'mu.XPO4_X498.SPLICE'=c(1, 1, 0, 0, 0, 1)
# )
# d.prep <- rbind(d.prep, d.prep)
# d.prep <- rbind(d.prep, d.prep)
# d.prep <- rbind(d.prep, d.prep)
# c.mu <- GetFeatures(d.prep, 'mu')
# c.cn <- GetFeatures(d.prep, 'cn')
# c.ge <- GetFeatures(d.prep, 'ge')
# d.prep <- Prep(d.prep, c(c.ge, c.cn), c.mu)

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

cat('Creating hyperparameter estimates')
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
    X.test.sml  <- predict(pp.sml, X.test)
    
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
    tuneLength=8, search='grid',
    trControl = trctrl(verboseIter=T)
  )
  
  set.seed(SEED); registerDoMC(3)
  m.svm.2 <- train(
    d$X.train.sml, d$y.train, method='svmRadial', 
    tuneGrid = data.frame(.sigma = as.numeric(svm.sigma[1]), .C = 2^(-3:3)),
    trControl = trctrl(verboseIter=T)
  )
  
  models <- list(
    svm.grid=m.svm.1,
    svm.est=m.svm.2,
    pls=m.pls.1
  )
  list(models=models, data=d)
}


