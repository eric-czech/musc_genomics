#'-----------------------------------------------------------------------------
#' Function library for filter methods
#'
#' This module contains utilities for scoring features and creating subset
#' models in resampling.
#' 
#' @author eczech
#'-----------------------------------------------------------------------------

GetFeatScoringFoldGen <- function(preproc, y.tresh, feat.limit=5000, n.core=4){    
  function(X.train.all, y.train, X.test, y.test){
    loginfo('Running feature scoring')
    registerDoMC(n.core)
    
    y.train.bin <- DichotomizeOutcome(y.train, y.tresh)
    y.test.bin  <- DichotomizeOutcome(y.test, y.tresh)
    
    loginfo('Creating feature scores')
    feat.scores <- foreach(feat=names(X.train.all), .combine=rbind)%dopar%{
      if (feat == 'origin') score <- 0
      else score <- FeatureScore(X.train.all[,feat], y.train.bin)
      data.frame(feature=feat, score=score, stringsAsFactors = F)
    }
    
    loginfo('Subsetting to top features')
    top.feats <- feat.scores %>% arrange(score) %>% 
      head(feat.limit) %>% as.data.frame %>% .$feature %>% as.character
    X.train <- X.train.all[,top.feats]
    
    loginfo('Applying feature preprocessing')
    pp <- preProcess(X.train, method=preproc)
    X.train <- predict(pp, newdata=X.train)
    X.test  <- predict(pp, X.test[,names(X.train)])
    
    list(
      feat.scores=feat.scores, feat.limit=feat.limit,
      preproc=pp, X.names=names(X.train.all), top.feats=top.feats,
      X.train=X.train, X.test=X.test,
      y.train=y.train, y.test=y.test,
      y.train.bin=y.train.bin, y.test.bin=y.test.bin
    )
  }
}

FilteringDataSummarizer <- function(){
  function(d){
    feat.types <- d$top.feats %>% sapply(function(x){
      if (x=='origin') 'origin'
      else if (str_detect(x, '^mu\\..*')) 'mutation'
      else if (str_detect(x, '^cn\\..*')) 'copy.number'
      else if (str_detect(x, '^ge\\..*')) 'gene.expression'
      else 'unknown'
    })
    feat.tbl <- table(feat.types)
    feat.cts <- feat.tbl %>% as.numeric 
    feat.cts <- paste(names(feat.tbl), feat.cts, sep='=', collapse=', ')
    logdebug('Training dimensions: %s', dim(d$X.train))
    logdebug('Validation dimensions: %s', dim(d$X.test))
    logdebug('Feature type counts: %s', feat.cts)
  }
}

TransformOriginSolidLiquid <- function(origin){
  ifelse(origin == 'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 0, 1)
}