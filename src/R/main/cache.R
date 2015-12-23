#'-----------------------------------------------------------------------------
#' Shared Data Cache
#'
#' This module hosts a single cache object containing lazyily-loaded objects
#' 
#' @author eczech
#'-----------------------------------------------------------------------------
source('utils.R')

DATA_CACHE <- list()

CacheEnabled <- function() coalesce(options()$enable.cache, T)

Fetch <- function(key, loader){
  #' Fetches an object from cache, loading it if necessary
  #' 
  #' Note that caching can be disabled using global option 'enable.cache'
  #' (e.g. options(enable.cache=F) -- Disables all caching)
  #' 
  #' Args: 
  #'    key: Name of object
  #'    loader: Function used to load and store object, if not already cached
  #' Returns:
  #'    The desired object
  if (!key %in% names(DATA_CACHE) || !CacheEnabled()){
    DATA_CACHE[[key]] <<- loader()
  }
  DATA_CACHE[[key]]
}

DEFAULT_CACHE_PATH <- file.path('~', 'genomics_data_cache')

FetchFromDisk <- function(filename, loader, dir=DEFAULT_CACHE_PATH){
  fpath <- file.path(dir, paste0(filename, '.Rdata'))
  if (!file.exists(dir))
    dir.create(dir)
  if (file.exists(fpath) && CacheEnabled()){
    env <- new.env()
    load(fpath, envir=env)
    env$res
  } else {
    res <- loader()
    save(res, file=fpath)
    res
  }
}
