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

CleanMutations <- function(m){
  str_replace_all(m, '[^\\w,]', '\\.') %>% 
    str_replace_all('_', '\\.') %>% str_trim %>% toupper
}

PrepareMutation <- function(d, gene){
  if (length(d) == 0)
    return(NULL)
  
  d <- CleanMutations(d)
  
  # Fetch a single feature value to be used for all-NA rows
  proto <- str_split(d, ',') %>% unlist %>% na.omit %>% head(1) 
  
  foreach(v=str_split(d, ','), i=icount(), .combine=rbind) %do% {
    if (length(v) == 1 && is.na(v))
      return(data.frame(feature=proto, value=0, i=i))
    v <- unique(v)
    if (length(v) > 1)
      v <- c(v, paste(v, collapse=':'))
    data.frame(feature=v, value=1, i=i)
  } %>% 
    mutate(feature=paste0(gene, '_', feature)) %>%
    dcast(i ~ feature, value.var='value', fun.aggregate=sum) %>% 
    select(-i)
}

RemoveRareMutations <- function(d, c.mu, min.mutations=3){
  rare.mutations <- which(apply(d[,c.mu], 2, sum) < min.mutations)
  if (length(rare.mutations) == 0) d
  else d %>% select(-one_of(c.mu[rare.mutations]))
}

GetUnivariateScores <- function(X, y, numeric.features, binary.features){
  nf <- foreach(feat=numeric.features, .combine=c) %dopar% gamScores(X[,feat], y)
  bf <- foreach(feat=binary.features, .combine=c) %dopar% anovaScores(y, factor(X[,feat]))
  rbind(
    data.frame(type='numeric', score=nf, feature=numeric.features),
    data.frame(type='binary', score=bf, feature=binary.features)
  ) %>% mutate(feature=as.character(feature))
}
#f.scores <- GetUnivariateScores(d.prep, 'response', c(c.cn, c.ge), c.mu)


# c('SDF,ADF', 'SDF,ADF,YYY', 'XXX', NA) %>% setNames(c('a', 'b', 'c')) %>% PrepareMutation() 
