library(devtools)
source_url('http://cdn.rawgit.com/eric-czech/portfolio/master/functional/common/R/utils.R')
source_url('http://cdn.rawgit.com/eric-czech/portfolio/master/functional/common/R/cache.R')
source_url('http://cdn.rawgit.com/eric-czech/portfolio/master/functional/ml/R/metrics.R')
#source('~/repos/portfolio/functional/common/R/utils.R')
lib('logging')
basicConfig(loglevels['DEBUG'])
CACHE_DIR <- '~/genomics_data_cache'