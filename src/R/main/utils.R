library(devtools)

##### Common Libraries #####

source('~/repos/portfolio/functional/common/R/utils.R')
#source_url('http://cdn.rawgit.com/eric-czech/portfolio/master/functional/common/R/utils.R')

# Load these packages in this order to avoid namespace conflicts and warning messages from plyr package
lib('MASS')
lib('plyr')
lib('dplyr')

source('~/repos/portfolio/functional/common/R/cache.R')
#source_url('http://cdn.rawgit.com/eric-czech/portfolio/master/functional/common/R/cache.R')

source('~/repos/portfolio/functional/ml/R/metrics.R')
#source_url('http://cdn.rawgit.com/eric-czech/portfolio/master/functional/ml/R/metrics.R')

##### Environment Settings #####
lib('logging')
basicConfig(loglevels['DEBUG'])
CACHE_DIR <- '~/genomics_data_cache'

RAW_CACHE <- Cache(dir=CACHE_DIR, project='raw_data')
TRAIN_CACHE <- Cache(dir=CACHE_DIR, project='training_data')
# RAW_CACHE$disable()

##### Data Settings #####

RESPONSE_TYPE <- NULL
RESPONSE_SELECTOR <- NULL
RESPONSE_THRESH <- NULL
SELECTION_THRESH <- NULL

EnableCosmic <- function(){
  RESPONSE_TYPE <<- 'cosmic' 
  RESPONSE_SELECTOR <<- function(d){d %>% filter(!is.na(ic_50)) %>% rename(response=ic_50) %>% select(-auc)}
  RESPONSE_THRESH <<- -1 
  SELECTION_THRESH <<- .001
}

EnableCtd <- function(){
  RESPONSE_TYPE <<- 'ctd'
  RESPONSE_SELECTOR <<- function(d){d %>% filter(!is.na(auc)) %>% rename(response=auc) %>% select(-ic_50)}
  RESPONSE_THRESH <<- -1
  SELECTION_THRESH <<- .0001
}



###### Visualization Utilities #####

plot.ly <- function(p) { p %>% plotly::config(showLink = F, displayModeBar = F) }

