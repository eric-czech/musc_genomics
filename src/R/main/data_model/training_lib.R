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

GetGlmnetLambda <- function(X.preproc, alpha=.01){
  #' Calculates estimated, optimal values for lambda glmnet parameter
  init <- glmnet(as.matrix(X.preproc), y, family = 'gaussian', nlambda = 100, alpha = alpha)
  lambda <- unique(init$lambda)
  lambda <- lambda[-c(1, length(lambda))]
  lambda <- lambda[1:length(lambda)]
  lambda
}

GetSvmSigma <- function(X.preproc, frac=1){
  sigest(as.matrix(X.preproc), frac = frac)
}


GetTrainingData <- function(cache, response.type, response.selector){
  loader <- function(){
    d <- GetPreparedData(cache)
    
    fields <- d$fields
    c.ge <- fields$gene_expression
    c.cn <- fields$copy_number
    c.mu <- fields$mutations
    d <- d$data
    
    d.prep <- response.selector(d) %>%
      #d %>% filter(!is.na(ic_50)) %>%                 # For now we're using COSMIC data only
      rename(response=ic_50) %>% select(-auc) %>%      # Select response field
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
