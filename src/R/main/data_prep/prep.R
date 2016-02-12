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

ApplyScaling<- function(d){
  foreach(f=names(d), .combine=cbind) %dopar% ScaleNumeric(d[,f])
}

RemoveNA <- function(d, row.threshold=.1){
  if (max(table(d$tumor_id)) != 1)
    stop('Found non-unique tumor id(s)')
  cts <- d %>% group_by(tumor_id) %>% do({
    r <- .
    data.frame(tumor_id=r$tumor_id[1], pct_na=sum(is.na(r))/ncol(r), stringsAsFactors=F)
  })
  rm.tumor <- cts %>% filter(pct_na >= row.threshold) %>% .$tumor_id
  rm.warn <- cts %>% filter(pct_na > 0 & pct_na < row.threshold) %>% .$tumor_id
  
  if (length(rm.warn) > 0)
    warning(sprintf(
      'The following tumors have NA predictor values but not enough to meet given 
      fractional (by-row) threshold of "%s": %s', row.threshold, paste(rm.warn, collapse=', ')))
  if (length(rm.tumor) > 0)
    loginfo('\nRemoving %s tumor records with high number of NA values', length(rm.tumor))
  d %>% filter(!tumor_id %in% rm.tumor)
}

GetPreparedData <- function(cache, raw.data=NULL, min.mutations=3){
  loader <- function(){
    if (!is.null(raw.data)) d <- raw.data
    else d <- GetRawData()
    
    # Determine genomic data fields
    c.mu <- GetFeatures(d, 'mu')
    c.ge <- GetFeatures(d, 'ge')
    c.cn <- GetFeatures(d, 'cn')
    
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
    c.mu <- GetFeatures(d, 'mu')
    
    # Return results
    list(data=d, fields=list(copy_number=c.cn, gene_expression=c.ge, mutations=c.mu))
  }
  cache$load('data_prep_01', loader)
  #FetchFromDisk('data_prep_01', loader) 
}

GetFeatures <- function(d, type){
  names(d)[str_detect(names(d), sprintf('^%s.', type))]
}

ApplyFeatureFilter <- function(X, y, numeric.features, binary.features, response.type, linear.only,
                               numeric.score.p=.05, binary.score.p=.05){
  f.scores <- GetUnivariateScores(X, y, numeric.features, binary.features, response.type, linear.only)
  rm.numeric <- f.scores %>% filter(type=='numeric' & score > numeric.score.p) %>% .$feature
  rm.binary  <- f.scores %>% filter(type=='binary' & score > binary.score.p) %>% .$feature  
  X %>% select(-one_of(c(rm.binary, rm.numeric)))
}


ApplyZeroVarianceFilter <- function(d){
  zero.var <- foreach(f=names(d), .combine=cbind) %dopar%{
    length(unique(d[,f])) == 1
  }
  d[,names(d)[!zero.var]]
}




