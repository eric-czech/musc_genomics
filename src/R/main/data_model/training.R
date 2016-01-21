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
lib('MASS')
lib('caret')

#d <- GetPreparedData()
d <- Boston


sbf(formula, data, sbfControl = sbfControl(), ...)

