#'-----------------------------------------------------------------------------
#' Filtering Model Final Settings
#'
#' This script contains settings for parameters that are meant to represent 
#' finalized versions of choices made in model training
#' 
#' @author eczech
#'-----------------------------------------------------------------------------

#' @title Final ensemble model subset
GetEnsembleModelNames <- function() c('rf', 'scrda', 'svm')

#' @title Final feature subset size
#final.feat.ct <- function() 20
GetOptimalFeatCt <- function() 20

#' @title Feature subset size for which to train holdout models over
GetHoldoutFeatCt <- function() c(10, 20, 30, 40)

#' @title Feature subset size for which to train models over entire known, labeled dataset
GetUnlabeledFeatCt <- function() c(10, 20, 30, 40)

#' @title Model name chosen to make final predictions
GetOptimalModel <- function() 'ensglm'