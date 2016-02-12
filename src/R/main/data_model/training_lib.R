#'-----------------------------------------------------------------------------
#' ML Model Training Library
#'
#' This module contains code for feature selection and preprocessing utilities
#' @author eczech
#'-----------------------------------------------------------------------------
source('utils.R')
source('data_prep/prep.R')
lib('caret')
lib('gam')
lib('glmnet')
lib('kernlab')



GetFeatureSelector <- function(num.threshold=.05, bin.threshold=.05){
  res <- caretSBF
  res$filter <- function(score, x, y) {
    if (length(unique(x)) == 2) score <= bin.threshold
    else score <= num.threshold
  }
  res$score <- function(x, y) {
    if (length(unique(x)) == 2) .anovaScore(y, factor(x))
    else .gamScore(x, y)
  }
  res
}

ScaleVector <- function(d){
  (d - mean(d, na.rm = T)) / (sd(d, na.rm=T))
}

GetGlmnetLambda <- function(X, y, alpha=.01, family='gaussian'){
  #' Calculates estimated, optimal values for lambda glmnet parameter
  init <- glmnet(as.matrix(X), y, nlambda = 100, alpha = alpha, family=family)
  lambda <- unique(init$lambda)
  lambda <- lambda[-c(1, length(lambda))]
  lambda <- lambda[1:length(lambda)]
  lambda
}

GetSvmSigma <- function(X, frac=1){
  sigest(as.matrix(X), frac = frac)
}


GetTrainingData <- function(cache, response.type, response.selector, min.mutations){
  loader <- function(){
    d <- GetPreparedData(cache)
    
    fields <- d$fields
    c.ge <- fields$gene_expression
    c.cn <- fields$copy_number
    c.mu <- fields$mutations
    d <- d$data
    
    d.prep <- response.selector(d) %>%
      Filter(function(x)!all(is.na(x)), .) %>%         # Remove NA-only columns
      RemoveNA(row.threshold=.1) %>%                   # Remove rows with large % NA's
      RemoveRareMutations(c.mu, min.mutations) %>% # Remove cols for rare mutations
      mutate(response=ScaleVector(response))           # Scale response
  
    if (any(is.na(d.prep)))
      stop('Dataset contains unexpected NA values')
    d.prep
  }
  cache$load(sprintf('data_prep_02_%s', response.type), loader)
}

DichotomizeOutcome <- function(y) {
  if (is.null(y)) return(NULL)
  factor((sign(y) + 1) * .5, levels=c(0, 1), labels=c('neg', 'pos'))
}

GetFoldDataGenerator <- function(preproc, linear.only, n.core=8, 
                                 sml.num.p=.0001, lrg.num.p=.15, sml.bin.p=.15, lrg.bin.p=.15){    
  function(X.train.all, y.train, X.test, y.test){
    # Apply feature selector
    loginfo('Running feature selection')
    registerDoMC(n.core)

    c.numeric <- c(GetFeatures(X.train.all, 'cn'), GetFeatures(X.train.all, 'ge'))
    c.binary <- GetFeatures(X.train.all, 'mu')
    
    # Note that feature selection is done using only numeric (not binary) responses
    X.train.sml <- ApplyFeatureFilter(X.train.all, y.train, c.numeric, c.binary, 'numeric', linear.only,
                                      numeric.score.p=sml.num.p, binary.score.p=sml.bin.p)
    X.train.lrg <- ApplyFeatureFilter(X.train.all, y.train, c.numeric, c.binary, 'numeric', linear.only,
                                      numeric.score.p=lrg.num.p, binary.score.p=lrg.bin.p)
    
    # Removing zero-variance features
    # X.train.sml <- ApplyZeroVarianceFilter(X.train.sml)
    # X.train.all <- ApplyZeroVarianceFilter(X.train.all)
    
    # Apply preprocessing to feature subset
    loginfo('Running preprocessing')
    pp.lrg <- preProcess(X.train.lrg, method=preproc)
    pp.sml <- preProcess(X.train.sml, method=preproc)
    X.train.lrg <- predict(pp.lrg, X.train.lrg)
    X.train.sml <- predict(pp.sml, X.train.sml)
    
    # X.test may be null if this data preprocessing call is not for resampling iteration
    if (is.null(X.test)){
      X.test.lrg <- NULL
      X.test.sml <- NULL
    } else {
      X.test.lrg  <- predict(pp.lrg, X.test[,names(X.train.lrg)])
      X.test.sml  <- predict(pp.sml, X.test[,names(X.train.sml)])
    }
    
    list(
      preproc=list(pp.lrg=pp.lrg, pp.sml=pp.sml), X.names=names(X.train.all),
      X.train.sml=X.train.sml, X.train.lrg=X.train.lrg,
      X.test.sml=X.test.sml, X.test.lrg=X.test.lrg,
      y.train=y.train, y.test=y.test,
      y.train.bin=DichotomizeOutcome(y.train), y.test.bin=DichotomizeOutcome(y.test)
    )
  }
}

CreateFoldIndex <- function(y, index, level){
  if (level == 1) # Outer folds should be returned for training set, not test set
    createFolds(y, k = 10, returnTrain = T)
  else if (level == 2) # Inner folds should be partioned the same as above
    createFolds(y[index], k=10, returnTrain = T)
  else
    stop(sprintf('Folds at level %s not supported', level))
}

GetDataSummarizer <- function(){
  function(d){
    gfl <- function(f, t) length(f[str_detect(f, sprintf('^%s.', t))])
    n.sml <- names(d$X.train.sml)
    n.lrg <- names(d$X.train.lrg)
    logdebug('Dimension reduction summary after preprocessing:')
    logdebug('Mutation features    : %s --> %s / %s', gfl(d$X.names, 'mu'), gfl(n.lrg, 'mu'), gfl(n.sml, 'mu'))
    logdebug('Expression features  : %s --> %s / %s', gfl(d$X.names, 'ge'), gfl(n.lrg, 'ge'), gfl(n.sml, 'ge'))
    logdebug('Copy Number features : %s --> %s / %s', gfl(d$X.names, 'cn'), gfl(n.lrg, 'cn'), gfl(n.sml, 'cn'))
  }
}

RESULT_METRICS <- c('auc', 'acc', 'tpr', 'tnr')
GetResultSummary <- function(d){
  pred <- prediction(d$y.pred.prob, d$y.test, label.ordering=c('neg', 'pos'))
  roc <- performance(pred, 'tpr', 'fpr') 
  auc <- performance(pred, 'auc')
  acc <- sum(d$y.pred.class == d$y.test) / length(d$y.pred.class)
  tpr <- sum(d$y.pred.class == d$y.test & d$y.test == 'pos') / sum(d$y.test == 'pos')
  tnr <- sum(d$y.pred.class == d$y.test & d$y.test == 'neg') / sum(d$y.test == 'neg')
  data.frame(
    x=roc@x.values[[1]], y=roc@y.values[[1]], 
    t=roc@alpha.values[[1]], auc=auc@y.values[[1]],
    acc=acc, tpr=tpr, tnr=tnr
  )
}

