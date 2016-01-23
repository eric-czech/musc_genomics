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


d <- GetPreparedData()

fields <- d$fields
c.ge <- fields$gene_expression
c.cn <- fields$copy_number
c.mu <- fields$mutations
d <- d$data

Prep <- function(d, c.numeric, c.binary, min.mutations=3){
  d %>% filter(!is.na(ic_50)) %>%                    # For now, we're using COSMIC data
    rename(response=ic_50) %>% select(-auc) %>%      # Select response field
    Filter(function(x)!all(is.na(x)), .) %>%         # Remove NA-only columns
    RemoveRareMutations(c.binary, min.mutations) %>% # Remove cols for rare mutations
    mutate_each_(funs(ScaleNumeric), c.numeric)      # Scale numeric features
}
d.prep <- Prep(d, c(c.ge, c.cn), c.mu)

# d <- Boston 

ctl <- sbfControl(functions = rfSBF)
sbf(medv ~ ., d)

f.score <- GetUnivariateScores(d, 'response', c(c.cn, c.ge), c.mu)
rm.numeric <- f.score %>% filter(type=='numeric' & score > .1) %>% .$feature
rm.binary  <- f.score %>% filter(type=='numeric' & score > .1) %>% .$feature



