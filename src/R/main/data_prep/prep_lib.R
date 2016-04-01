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
  rare.mutations <- which(apply(d[, c.mu], 2, sum) < min.mutations)
  if (length(rare.mutations) == 0) d
  else d %>% select(-one_of(c.mu[rare.mutations]))
}


.NotValid <- function(pv) any(class(pv) == "try-error") || is.na(pv) || is.nan(pv)

.GetFeatureScore <- function(y, x, response.type, feature.type, 
                             linear.only=T, feature.name='Unspecified'){
  # Modes = c('binary-response', ')
  if (!response.type %in% c('numeric', 'factor', 'binary')){
    stop(sprintf('Response type "%s" is not valid.', response.type))
  }

  # Return p-value of 1 immediately for zero variance predictors
  if (length(unique(x)) == 1) return(1) 
  
  pv <- NA
  if (response.type == 'numeric'){
    if (!linear.only && feature.type == 'numeric'){
      # This is necessary to avoid conflics with older gam package;
      # see http://stackoverflow.com/questions/20693865/s-in-mgcv-doesnt-work-when-vgam-package-loaded
      # *Note the original error encountered was "$ operator is invalid for atomic vectors" with valid x and y
      s<-mgcv:::s
      pv <- try(anova(mgcv::gam(y ~ s(x)))$s.pv, silent=T)
      if (.NotValid(pv)) {
        warning('Error occurred in Gam fit')
        browser()
      } 
    }
    if (.NotValid(pv)) pv <- try(anova(lm(y ~ x), test = "F")[1, "Pr(>F)"], silent=T)
  } else {
    if (feature.type == 'numeric') pv <- try(anova(lm(x ~ y), test = "F")[1, "Pr(>F)"], silent=T)
    else pv <- try(chisq.test(x, y)$p.value, silent=T)
  }
  if (.NotValid(pv)) {
    warning(sprintf('Error occurred in determining p-value for feature "%s"', feature.name))
    browser()
    pv <- 0 # Return 0 on error so no potentially good predictors are lost
  }
  pv
}

GetUnivariateScores <- function(X, y, numeric.features, binary.features, response.type, linear.only){
  nf <- foreach(feat=numeric.features, .combine=c, .errorhandling='stop') %dopar% 
    .GetFeatureScore(y, X[,feat], response.type, 'numeric', linear.only=linear.only)
  bf <- foreach(feat=binary.features, .combine=c, .errorhandling='stop') %dopar% 
    .GetFeatureScore(y, X[,feat], response.type, 'binary', linear.only=linear.only)
  rbind(
    data.frame(type='numeric', score=nf, feature=numeric.features),
    data.frame(type='binary', score=bf, feature=binary.features)
  ) %>% mutate(feature=as.character(feature))
}

#f.scores <- GetUnivariateScores(d.prep, 'response', c(c.cn, c.ge), c.mu)


# c('SDF,ADF', 'SDF,ADF,YYY', 'XXX', NA) %>% setNames(c('a', 'b', 'c')) %>% PrepareMutation() 
