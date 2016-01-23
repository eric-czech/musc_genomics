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
lib('doMC')
lib('gam')
registerDoMC(12)

# options(enable.cache=T)
# options(enable.cache=F)

# Scaling function that ignores NA values
ScaleNumeric <- function(x) (x - mean(x, na.rm=T)) / sd(x, na.rm=T)

GetPreparedData <- function(raw.data=NULL, min.mutations=3){
  loader <- function(){
    if (!is.null(raw.data)) d <- raw.data
    else d <- GetRawData()
    
    # Determine genomic data fields
    c.mu <- names(d)[str_detect(names(d), '^mu.')]
    c.ge <- names(d)[str_detect(names(d), '^ge.')]
    c.cn <- names(d)[str_detect(names(d), '^cn.')]
    
    # Determine remaining fields like tumor id, AUC, and IC 50
    c.id <- setdiff(names(d), c(c.mu, c.ge, c.cn))
    
    # Dummy encode mutation values
    d.mu <- foreach(x=c.mu, .combine=cbind) %dopar% PrepareMutation(d[,x], x)
    
    # Verify that feature name cleaning did not result in collisions
    # *Note that it was originally discovered that gene 'CDCP2' had a mutation value
    # equal to 'M409Hfs*45,M409Hfs*45' which is clearly a reporting error and
    # one reason why this check is necessary
    if (max(d.mu) > 1) stop('Encoding of mutations resulted in feature collisions (should not be possible)')
    
    # Verify that no records were lost or gained through encoding
    if (nrow(d.mu) != nrow(d)) stop('Number of rows in encoded mutation data do not match expected value')
    
    # Restrict data to containing only mutations that occur at least min.mutations times
    d.mu <- RemoveRareMutations(d.mu, names(d.mu), min.mutations)
    
    # Merge mutation data with all other data and reset mutation field names
    d <- d %>% select(-one_of(c.mu)) %>% cbind(d.mu)
    c.mu <- names(d)[str_detect(names(d), '^mu.')]
    
    # Return results
    list(data=d, fields=list(copy_number=c.cn, gene_expression=c.ge, mutations=c.mu))
  }
  FetchFromDisk('data_prep_01', loader) 
}


GetUnivariateScores <- function(d, response, numeric.features, binary.features){
  y <- d[,response]
  nf <- foreach(feat=numeric.features, .combine=c) %do% gamScores(d[,feat], y)
  bf <- foreach(feat=binary.features, .combine=c) %do% anovaScores(y, factor(d[,feat]))
  rbind(
    data.frame(type='numeric', score=nf, feature=numeric.features),
    data.frame(type='binary', score=bf, feature=binary.features)
  )
}

# d <- GetPreparedData()

# mutate(ic_50_s=scale.v(ic_50), auc_s=scale.v(auc))

# d[,c('tumor_id', 'ic_50', 'auc')] %>% head


