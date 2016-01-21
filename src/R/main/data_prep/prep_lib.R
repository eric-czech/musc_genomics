#'-----------------------------------------------------------------------------
#' Function library for data preparation actions
#'
#' These functions will help augment different data encoding schemes (mostly
#' centered around how mutations are handled) as well as expose helper functions
#' for scaling and any other actions taken prior to modeling.
#' 
#' @author eczech
#'-----------------------------------------------------------------------------
source('utils.R')
lib('dplyr')
lib('foreach')
lib('iterators')
lib('stringr')

CleanMutation <- function(m){
  paste(str_extract_all(m, '\\w')[[1]], collapse='') %>% str_replace('_', '') %>% toupper
}

PrepareMutation <- function(d){
  if (length(d) == 0)
    return(NULL)
  
  # Fetch a single feature value to be used for all-NA rows
  proto <- str_split(d, ',') %>% unlist %>% na.omit %>% head(1) %>% CleanMutation
  
  foreach(v=str_split(d, ','), i=icount(), .combine=rbind) %do% {
    if (length(v) == 1 && is.na(v))
      return(data.frame(feature=proto, value=0, i=i))
    v <- sapply(v, CleanMutation)
    if (length(v) > 1)
      v <- c(v, paste(v, collapse=':'))
    data.frame(feature=v, value=1, i=i)
  } %>% 
    dcast(i ~ feature, value.var='value', fun.aggregate=sum) %>% 
    select(-i)
}

# c('SDF,ADF', 'SDF,ADF,YYY', 'XXX', NA) %>% setNames(c('a', 'b', 'c')) %>% PrepareMutation() 
