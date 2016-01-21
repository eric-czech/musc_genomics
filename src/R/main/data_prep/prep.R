#'-----------------------------------------------------------------------------
#' Genomics data importation script
#'
#' This script will collect data from the cBioPortal database and merge this
#' information with sensitivity data in the COSMIC and CTD2 databases.
#' 
#' @author eczech
#'-----------------------------------------------------------------------------
source('utils.R')
source('data_import/import.R')
source('data_prep/prep_lib.R')
lib('dplyr')
lib('reshape2')

# options(enable.cache=T)
# options(enable.cache=F)

GetPreparedData <- function(raw.data=NULL, scale.numeric.fields=T, min.mutations=3){
  if (!is.null(raw.data)) d <- raw.data
  else d <- GetRawData()
  
  # Scaling function that ignores NA values
  scale.v <- function(x) (x - mean(x, na.rm=T)) / sd(x, na.rm=T)
  
  # Determine genomic data fields
  c.mu <- names(d)[str_detect(names(d), '^mu.')]
  c.ge <- names(d)[str_detect(names(d), '^ge.')]
  c.cn <- names(d)[str_detect(names(d), '^cn.')]
  
  # Determine remaining fields like tumor id, AUC, and IC 50
  c.id <- setdiff(names(d), c(c.mu, c.ge, c.cn))
  
  # Scale all numeric fields, if requested
  if (scale.numeric.fields) d <- d %>% mutate_each_(funs(scale.v), c(c.ge, c.cn))
  
  # Dummy encode mutation values
  d.mu <- foreach(x=c.mu, .combine=cbind) %do% PrepareMutation(d[,x])
  
  # Verify that feature name cleaning did not result in collisions
  if (max(d.mu) > 1) stop('Encoding of mutations resulted in feature collisions (should not be possible)')
  
  # Verify that no records were lost or gained through encoding
  if (nrow(d.mu) != nrow(d)) stop('Number of rows in encoded mutation data do not match expected value')
  
  # Restrict data to containing only mutations that occur at least min.mutations times
  frequent.mutations <- which(apply(d.mu, 2, sum) >= min.mutations)
  d.mu <- d.mu %>% select(one_of(names(d.mu)[frequent.mutations]))
  
  # Merge mutation data with all other data and reset mutation field names
  d <- d %>% select(-one_of(c.mu)) %>% cbind(d.mu %>% setNames(paste0('mu.', names(d.mu))))
  c.mu <- names(d)[str_detect(names(d), '^mu.')]
  
  # Return results
  list(data=d, fields=list(copy_number=c.cn, gene_expression=c.ge, mutations=c.mu))
}

d <- GetPreparedData()

# mutate(ic_50_s=scale.v(ic_50), auc_s=scale.v(auc))

# d[,c('tumor_id', 'ic_50', 'auc')] %>% head


