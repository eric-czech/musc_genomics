#'-----------------------------------------------------------------------------
#' ML Model Training Library
#'
#' This module contains code for feature selection and preprocessing utilities
#' @author eczech
#'-----------------------------------------------------------------------------
source('utils.R')
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