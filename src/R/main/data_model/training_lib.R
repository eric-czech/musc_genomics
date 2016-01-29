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

.anovaScore <- function(x, y) {
  pv <- try(anova(lm(x ~ y), test = "F")[1, "Pr(>F)"], silent = TRUE)
  if (any(class(pv) == "try-error") || is.na(pv) || is.nan(pv)) 
    pv <- 0 # Return 0 on error so no potentially good predictors are lost
  pv
}

.gamScore <- function(x, y){
  browser()
  pv <- try(anova(gam::gam(y ~ s(x)))[2, "Pr(F)"], 
            silent = TRUE)
  if (any(class(pv) == "try-error")) 
    pv <- try(anova(lm(x ~ y), test = "F")[1, "Pr(>F)"], 
              silent = TRUE)
  if (any(class(pv) == "try-error") || is.na(pv) || is.nan(pv)) 
    pv <- 0 # Return 0 on error so no potentially good predictors are lost
  pv
}

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

GetGlmnetLambda <- function(X, y, alpha=.01){
  #' Calculates estimated, optimal values for lambda glmnet parameter
  init <- glmnet(as.matrix(X), y, family = 'gaussian', nlambda = 100, alpha = alpha)
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


GetFoldDataGenerator <- function(preproc){    
  function(X.train.all, y.train, X.test, y.test){
    # Apply feature selector
    loginfo('Running feature selection')
    registerDoMC(8)
    X.dim <- dim(X.train.all)
    c.numeric <- c(GetFeatures(X.train.all, 'cn'), GetFeatures(X.train.all, 'ge'))
    c.binary <- GetFeatures(X.train.all, 'mu')
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
}

GetDataSummarizer <- function(){
  function(d){
    gfl <- function(d, t) length(GetFeatures(d, t))
    logdebug('Dimension reduction summary after preprocessing:')
    logdebug('Mutation features   : %s --> %s', gfl(d$X.train.all, 'mu'), gfl(d$X.train.sml, 'mu'))
    logdebug('Expression features : %s --> %s', gfl(d$X.train.all, 'ge'), gfl(d$X.train.sml, 'ge'))
    logdebug('Copy Number features: %s --> %s', gfl(d$X.train.all, 'cn'), gfl(d$X.train.sml, 'cn'))
  }
}