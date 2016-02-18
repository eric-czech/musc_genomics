#'-----------------------------------------------------------------------------
#' Pre-Training Analysis Utilities
#'
#' This module contains code used to examine prepared sensitivity data to 
#' determine things like optimal cutoffs for training classifiers.
#' 
#' @author eczech
#'-----------------------------------------------------------------------------
source('utils.R')
source('data_model/training_lib.R')
source('data_model/training_viz.R')
lib('dplyr')
lib('ggplot2')

# Examine CTD2 and COSMIC response distributions prior to pre-processing

d.ctd2 <- GetCTD2V2Data() 
d.cosmic <- GetCOSMICData()

rbind(
  data.frame(type='COSMIC IC50', response=d.cosmic$ic_50),
  data.frame(type='CTD2 AUC', response=d.ctd2$auc)
) %>% ggplot(aes(x=response)) +
  geom_histogram(binwidth=.2, alpha=.3) + 
  geom_density(aes(y=.2 * ..count..)) + theme_bw() + 
  facet_wrap(~type, scales='free') + ggtitle('Sensitivity Distributions (on original scales)')
