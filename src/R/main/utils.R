library(devtools)
#source_url('http://cdn.rawgit.com/eric-czech/portfolio/master/functional/common/R/utils.R')
source_url('http://cdn.rawgit.com/eric-czech/portfolio/master/functional/ml/R/metrics.R')
#source_url('http://cdn.rawgit.com/eric-czech/portfolio/master/functional/common/R/cache.R')
source('~/repos/portfolio/functional/common/R/utils.R')
source('~/repos/portfolio/functional/common/R/cache.R')
#source('~/repos/portfolio/functional/ml/R/metrics.R')

lib('logging')
basicConfig(loglevels['DEBUG'])
CACHE_DIR <- '~/genomics_data_cache'

RAW_CACHE <- Cache(dir=CACHE_DIR, project='raw_data')
TRAIN_CACHE <- Cache(dir=CACHE_DIR, project='training_data')
MODEL_CACHE <- Cache(dir=CACHE_DIR, project='training_models')
# RAW_CACHE$disable()