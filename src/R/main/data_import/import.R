#'-----------------------------------------------------------------------------
#' Genomics data importation script
#'
#' This script will collect data from the cBioPortal database and merge this
#' information with sensitivity data in the COSMIC and CTD2 databases.
#' 
#' @author eczech
#'-----------------------------------------------------------------------------
library(devtools)
source_url('https://cdn.rawgit.com/eric-czech/portfolio/master/demonstrative/R/common/utils.R')
source('data_import/import_lib.R')
lib('cgdsr')
lib('foreach')
lib('dplyr')
lib('iterators')
lib('reshape2')

# options(enable.cache=T)
# options(enable.cache=F)

#---------------------------#
# Load cBioPortal CCLE Data #
#---------------------------#

d.biop <- GetBioPortalData()
d.ctd2 <- GetCTD2V2Data()

d.cosmic <- GetCOSMICData()

scale.v <- function(x) (x - mean(x, na.rm=T)) / sd(x, na.rm=T)

d <- d.biop %>% 
  left_join(d.ctd2, by='tumor_id') %>% 
  left_join(d.cosmic, by='tumor_id') %>%
  select(tumor_id, ic_50, auc, everything())

# Determine genomic data fields
c.mu <- names(d)[str_detect(names(d), '^mu.')]
c.ge <- names(d)[str_detect(names(d), '^ge.')]
c.cn <- names(d)[str_detect(names(d), '^cn.')]

# Determine remaining fields like tumor id, AUC, and IC 50
c.id <- setdiff(names(d), c(c.mu, c.ge, c.cn))

# Scale all numeric fields
# sapply(d[,c.cn], class) %>% table
d.trans <- d %>% mutate_each_(funs(scale.v), c(c.ge, c.cn))


PrepareMutation <- function(d){
  if (length(d) == 0)
    return(NULL)
  
  # Fetch a single feature value to be used for all-NA rows
  proto <- str_split(d, ',') %>% unlist %>% na.omit %>% head(1)
  
  # TODO: Add regex filtering for chars in mutation names
  
  foreach(v=str_split(d, ','), i=icount(), .combine=rbind) %do% {
    if (length(v) == 1 && is.na(v))
      return(data.frame(feature=proto, value=1, i=i))
    if (length(v) > 1)
      v <- c(v, paste(v, collapse=':'))
    data.frame(feature=v, value=1, i=i)
  } %>% 
    dcast(i ~ feature, value.var='value', fun.aggregate=length) %>% 
    select(-i)
}

# c('SDF,ADF', 'SDF,ADF,YYY', 'XXX', NA) %>% setNames(c('a', 'b', 'c')) %>% PrepareMutation()

d.mu <- foreach(x=c.mu, .combine=cbind) %do% {
  PrepareMutation(d.trans[,x])
}


# mutate(ic_50_s=scale.v(ic_50), auc_s=scale.v(auc))

# d[,c('tumor_id', 'ic_50', 'auc')] %>% head

