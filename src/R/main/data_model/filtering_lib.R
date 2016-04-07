#'-----------------------------------------------------------------------------
#' Function library for filter methods
#'
#' This module contains utilities for scoring features and creating subset
#' models in resampling.
#' 
#' @author eczech
#'-----------------------------------------------------------------------------


SplitFeatValues <- function(x, y){
  lvl <- levels(y)
  x1 <- x[y == lvl[1]]
  x2 <- x[y == lvl[2]]
  list(x1, x2)
}

GeneInAgreement <- function(x1, x2, y){
  lvl <- levels(y)
  s1 <- SplitFeatValues(x1, y)
  d1 <- mean(s1[[1]]) - mean(s1[[2]])
  s2 <- SplitFeatValues(x2, y)
  d2 <- mean(s2[[1]]) - mean(s2[[2]])
  sign(d1) == sign(d2)
}

GeneScore <- function(feat, X, y, t.test.alt='two.sided') {
  if (nlevels(y) != 2)
    stop('Response factor (y) must be a 2-level factor')
  
  x <- X[,feat]
  parts <- str_split(feat, '\\.')[[1]]
  gene <- paste(parts[2:length(parts)], collapse='.')
  type <- parts[1]
  
  # Use fisher test when x is a factor (or has low cardinality)
  if (is.factor(x) || length(unique(x)) <= 2){
    pv <- try(fisher.test(factor(x), y)$p.value, silent = TRUE)
  
  # Use t-test otherwise
  } else { 
    
    # Find the complementary feature for this gene
    if (!type %in% c('cn', 'ge'))
      stop(sprintf('Found unexpected numeric feature "%s"', feat))
    comp.feat <- ifelse(type == 'cn', 'ge', 'cn')
    comp.feat <- paste(comp.feat, gene, sep='.')
    
    agree <- F
    if (comp.feat %in% names(X)){
      comp.x <- X[,comp.feat]
      agree <- GeneInAgreement(x, comp.x, y)
    }
    
    # If gene features are in agreement with respect to direction, run
    # a t-test on the original feature
    if (agree){
      # Split x values into two groups based on response
      s <- SplitFeatValues(x, y)
      pv <- try(t.test(s[[1]], s[[2]], alternative=t.test.alt)$p.value, silent = TRUE) 
    # Otherwise, return the worst possible score
    } else {
      pv <- 1
    }
  }
  
  # Revert to worst possible score in the event of an error
  if (any(class(pv) == "try-error") || is.na(pv) || is.nan(pv)) pv <- 1
  
  pv
}

GetFeatScoringFoldGen <- function(preproc, y.tresh, feat.limit=5000, n.core=4,
                                  t.test.alt='greater'){    
  function(X.train.all, y.train, X.test, y.test){
    loginfo('Running feature scoring')
    registerDoMC(n.core)
    
    y.train.bin <- DichotomizeOutcome(y.train, y.tresh)
    y.test.bin  <- DichotomizeOutcome(y.test, y.tresh)
    
    loginfo('Creating feature scores')
    feat.scores <- foreach(feat=names(X.train.all), .combine=rbind)%dopar%{
      if (feat == 'origin') {
        score <- 0
      } else {
        #score <- FeatureScore(X.train.all[,feat], y.train.bin, t.test.alt=t.test.alt)
        score <- GeneScore(feat, X.train.all, y.train.bin, t.test.alt=t.test.alt)
      }
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

TransformOriginSolidLiquid <- list(
  convert = function(origin){
    ifelse(origin == 'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 0, 1)
  }, 
  max.val = 1
)

TransformOriginMostFrequent <- list(
  convert = function(origin){
    sapply(origin, function(x){
      if (x == 'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE') 1
      else if (x == 'LUNG') 2
      else if (x == 'BREAST') 3
      else if (x == 'CENTRAL_NERVOUS_SYSTEM') 4
      else if (x == 'SKIN') 5
      else 0
    })
  }, max.val = 5
)

GetModelNamesByFeatCount <- function(models, feat.ct, origin.name){
  m.subset <- sprintf('\\.%s\\.%s', feat.ct, origin.name)
  all.model.names <- unlist(lapply(names(models), function(m) names(models[[m]])))
  sub.model.names <- all.model.names[all.model.names %>% str_detect(m.subset)]
  if (length(sub.model.names) == 0)
    stop(sprintf('Failed to find models matching top feat count %s', feat.ct))
  top.models <- foreach(m=names(models), combine=c) %do% {
    m.subset <- names(models[[m]])[names(models[[m]]) %in% sub.model.names]
    if (length(m.subset) != 1 || sum(is.na(m.subset)) > 0)
      stop('Failed to find best subset model')
    m.subset
  }
  unlist(top.models)
}


GetModelsByFeatCount <- function(models, feat.ct, origin.name){
  top.model.names <- GetModelNamesByFeatCount(models, feat.ct, origin.name)
  top.models <- foreach(m=names(models)) %do% {
    m.subset <- names(models[[m]])[names(models[[m]]) %in% top.model.names]
    if (length(m.subset) != 1 || sum(is.na(m.subset)) > 0)
      stop('Failed to find best subset model')
    models[[m]][[m.subset]]
  } %>% setNames(names(models))
}

GetEnsModelDefinition <- function(name, sub.models, n.feat=NULL, ...){
  n.feat <- if (is.null(n.feat)) n.feats else n.feat
  lapply(n.feat, function(x){ 
    feat.ct.models <- GetModelsByFeatCount(sub.models, x, origin.name)
    list(model=GetEnsembleModelForTransform(x, feat.ct.models, name, ...), max.feats=x)
  })
}

GetEnsembleModelDefFromSubModelFit <- function(models, ens.model.names, origin.trans, origin.name, n.feat=NULL){
  sub.models <- models[ens.model.names]
  if (any(is.na(names(sub.models))))
    stop(sprintf('Failed to some or all of the following models for ensembling: %s', paste(ens.model.names, collapse=', ')))
  ens.models.def <- list()
  ens.models.def$ensglm <- GetEnsModelDefinition(
    'ensglm', sub.models, method='glm', 
    origin.transform=origin.trans, origin.name=origin.name, n.feat=n.feat
  )
  ens.models.def$ensavg <- GetEnsModelDefinition(
    'ensavg', sub.models, method=GetEnsembleAveragingModel(), 
    origin.transform=origin.trans, origin.name=origin.name, n.feat=n.feat
  )
  ens.models.def
}

ExtractTopFeatures <- function(models){
  # Extract the first model in the results to use as the representative from
  # which the top feature list will be drawn (at the top level, each list item
  # will correspond to per-feature results for one ML algo, e.g. svm).  Note
  # this this feature list, given any number of those features, will be the 
  # same for all models
  rep.model <- models[[1]]
  
  # Loop through all the per-feature-count model results contained in the
  # representative model list and extract the feature names associated with
  # that feature subset size
  res <- lapply(rep.model, function(m) {
    feats <- m$fit$top.feats
    data.frame(feature=feats, n.feats=length(feats))
  }) %>% do.call(rbind, .)
  
  if (!is.data.frame(res) || nrow(res) == 0)
    stop('Extracted feature set is not a data frame or is empty')
  rownames(res) <- NULL
  
  n.invalid <- res %>% group_by(n.feats) %>% tally %>% filter(n.feats != n) %>% nrow
  if (n.invalid > 0)
    stop('At least one extracted feature subset did not have the expected number of feature names')
  
  res
}

RunHoldoutFit <- function(trainer, models.def, d.tr, d.ho, fit.key, dat.key, dat.gen, n.feat=NULL){
  trainer$getCache()$load(fit.key, function(){
    models <- list()
    for (m in names(models.def)){
      for (m.feat in models.def[[m]]){
        m.name <- m.feat$model$name
        m.n.feat <- as.numeric(str_split(m.name, '\\.')[[1]][2])
        if (is.na(m.n.feat))
          stop(sprintf('Failed to parse number of features from model name "%s"', m.name))
        if (!is.null(n.feat) && !m.n.feat %in% n.feat)
          next
        models[[m.feat$model$name]] <- m.feat$model
      }
    }
    trainer$holdout(models, d.tr$X, d.tr$y, d.ho$X, d.ho$y, dat.gen, dat.key) 
  })
}

RestackHoldoutFit <- function(models.by.feat.ct){
  res <- list()
  for (m.feat in names(models.by.feat.ct)){
    model.class <- str_split(m.feat, '\\.')[[1]][1]
    if (is.null(model.class) || length(model.class) != 1)
      stop(sprintf('Failed to parse model class from sub model name "%s"', m.feat))
    if (!model.class %in% names(res))
      res[[model.class]] <- list()
    res[[model.class]][[m.feat]] <- models.by.feat.ct[[m.feat]]
  }
  res
}

ValidatePredictions <- function(models){
  bad.models <- c()
  for (m in names(models)){
    #if (m == 'svm' || m == 'scrda') next
    for (m.feat in names(models[[m]])){
      for (m.fold in seq_along(models[[m]][[m.feat]])){
        model <- models[[m]][[m.feat]][[m.fold]]
        if ('caretStack' %in% class(model$fit))
          next # Ignore ensemble models
        pred <- model$fit$pred
        if (is.null(pred)){
          msg <- sprintf('Predictions NULL for model "%s", subset "%s", fold "%s"', m, m.feat, m.fold)
          stop(msg)
        }
        
        if (any(is.na(pred))){
          msg <- sprintf('Found NAs in prediction matrix for model "%s", subset "%s", fold "%s"', m, m.feat, m.fold)
          bad.models <- c(bad.models, msg)
        }
      }
    }
  }
  if (length(bad.models) > 0)
    stop(sprintf('The following models + folds had NA prediction values:\n\t%s', paste(bad.models, collapse='\n\t')))
}