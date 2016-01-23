#'-----------------------------------------------------------------------------
#' ML Model Training Functions
#'
#' This module contains code for feature selection and preprocessing utilities
#' @author eczech
#'-----------------------------------------------------------------------------
source('utils.R')
lib('caret')
lib('gam')

.anovaScore <- function(x, y) {
  pv <- try(anova(lm(x ~ y), test = "F")[1, "Pr(>F)"], silent = TRUE)
  if (any(class(pv) == "try-error") || is.na(pv) || is.nan(pv)) 
    pv <- 0 # Return 0 on error so no potentially good predictors are lost
  pv
}

.gamScore <- function(x, y){
  pv <- try(anova(gam::gam(y ~ s(x)), test = "F")[2, "Pr(F)"], 
            silent = TRUE)
  if (any(class(pv) == "try-error")) 
    pv <- try(anova(lm(x ~ y), test = "F")[1, "Pr(>F)"], 
              silent = TRUE)
  if (any(class(pv) == "try-error") || is.na(pv) || is.nan(pv)) 
    pv <- 0 # Return 0 on error so no potentially good predictors are lost
  pv
}

GetFeatureSelector <- function(score.threshold=.05){
  res <- caretSBF
  res$filter <- function(score, x, y) {
    score <= score.threshold
  }
  res$score <- function(x, y) {
    if (length(unique(x)) == 2) .anovaScore(y, factor(x))
    else .gamScore(x, y)
  }
  res
}