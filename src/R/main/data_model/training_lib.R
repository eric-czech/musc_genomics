#'-----------------------------------------------------------------------------
#' ML Model Training Library
#'
#' This module contains code for feature selection and preprocessing utilities
#' @author eczech
#'-----------------------------------------------------------------------------
source('utils.R')
source('data_prep/prep.R')
lib('caret')
lib('gam')
lib('glmnet')
lib('kernlab')
lib('ROCR')



GetFeatureSelector <- function(num.threshold=.05, bin.threshold=.05){
  res <- caretSBF
  res$filter <- function(score, x, y) {
    if (length(unique(x)) == 2) score <= bin.threshold
    else score <= num.threshold
  }
  res$score <- function(x, y) {
    if (length(unique(x)) == 2) .anovaScore(y, factor(x))
    else .gamScore(x, y)
  }
  res
}

ScaleVector <- function(d){
  (d - mean(d, na.rm = T)) / (sd(d, na.rm=T))
}

GetGlmnetLambda <- function(X, y, alpha=.01, family='gaussian'){
  #' Calculates estimated, optimal values for lambda glmnet parameter
  init <- glmnet(as.matrix(X), y, nlambda = 100, alpha = alpha, family=family)
  lambda <- unique(init$lambda)
  lambda <- lambda[-c(1, length(lambda))]
  lambda <- lambda[1:length(lambda)]
  lambda
}

GetSvmSigma <- function(X, frac=1){
  sigest(as.matrix(X), frac = frac)
}


GetTrainingData <- function(cache, response.type, response.selector, min.mutations){
  loader <- function(){
    d <- GetPreparedData(cache)

    fields <- d$fields
    c.ge <- fields$gene_expression
    c.cn <- fields$copy_number
    c.mu <- fields$mutations
    d <- d$data
    
    d.prep <- response.selector(d) %>%
      Filter(function(x)!all(is.na(x)), .) %>%         # Remove NA-only columns
      RemoveNA(row.threshold=.1) %>%                   # Remove rows with large % NA's
      RemoveRareMutations(c.mu, min.mutations) %>% # Remove cols for rare mutations
      mutate(response=ScaleVector(response))           # Scale response
    
    if (any(is.na(d.prep)))
      stop('Dataset contains unexpected NA values')
    d.prep
  }
  cache$load(sprintf('data_prep_02_%s', response.type), loader)
}

DichotomizeOutcome <- function(y, threshold) {
  if (is.null(y)) return(NULL)
  if (any(is.na(y))) stop('Outcome cannot contain NA values')
  
  # If outcome value is below or equal to 'threshold', then the drug 
  # sensitivity is high and the outcome is labeled as "pos" ("neg" otherwise)
  #factor(ifelse(y <= threshold, 'pos', 'neg'), levels=c('neg', 'pos'))
  factor(ifelse(y <= threshold, 'pos', 'neg'), levels=c('pos', 'neg'))
}

GetFoldDataGenerator <- function(preproc, y.tresh, linear.only, n.core=8, 
                                 sml.num.p=.0001, lrg.num.p=.15, 
                                 sml.bin.p=.15, lrg.bin.p=.15, pls.comp=500){    
  function(X.train.all, y.train, X.test, y.test){
    # Apply feature selector
    loginfo('Running feature selection')
    registerDoMC(n.core)
    
    c.cn <- GetFeatures(X.train.all, 'cn')
    c.ge <- GetFeatures(X.train.all, 'ge')
    c.numeric <- c(c.cn, c.ge)
    c.binary <- GetFeatures(X.train.all, 'mu')
    
    # Note that feature selection is done using only numeric (not binary) responses
    X.train.sml <- ApplyFeatureFilter(X.train.all, y.train, c.numeric, c.binary, 'numeric', linear.only,
                                      numeric.score.p=sml.num.p, binary.score.p=sml.bin.p)
    X.train.lrg <- ApplyFeatureFilter(X.train.all, y.train, c.numeric, c.binary, 'numeric', linear.only,
                                      numeric.score.p=lrg.num.p, binary.score.p=lrg.bin.p)
    
    # Apply preprocessing to feature subset
    loginfo('Running preprocessing')
    
    ## Standard preprocessing
    pp.lrg <- preProcess(X.train.lrg, method=preproc)
    pp.sml <- preProcess(X.train.sml, method=preproc)
    
    X.train.lrg <- predict(pp.lrg, X.train.lrg)
    X.train.sml <- predict(pp.sml, X.train.sml)
    
    # X.test may be null if this data preprocessing call is not for resampling iteration
    if (is.null(X.test)){
      X.test.lrg <- NULL
      X.test.sml <- NULL
    } else {
      X.test.lrg  <- predict(pp.lrg, X.test[,names(X.train.lrg)])
      X.test.sml  <- predict(pp.sml, X.test[,names(X.train.sml)])
    }
    
    list(
      preproc=list(pp.lrg=pp.lrg, pp.sml=pp.sml), 
      X.names=names(X.train.all),
      X.train.sml=X.train.sml, X.train.lrg=X.train.lrg,
      X.test.sml=X.test.sml, X.test.lrg=X.test.lrg,
      y.train=y.train, y.test=y.test,
      y.train.bin=DichotomizeOutcome(y.train, y.tresh), y.test.bin=DichotomizeOutcome(y.test, y.tresh)
    )
  }
}

# Graveyard for above:

## PCA preprocessing
#     pp.pca.ge <- preProcess(X.train.all[,c.ge], method=c(preproc, 'pca'), thresh=pca.thresh)
#     pp.pca.cn <- preProcess(X.train.all[,c.cn], method=c(preproc, 'pca'), thresh=pca.thresh)
#     pred.pca <- function(d, c.sml){ cbind(
#       predict(pp.pca.ge, d[,c.ge]) %>% setNames(paste0('ge.', names(.))),
#       predict(pp.pca.cn, d[,c.cn]) %>% setNames(paste0('cn.', names(.))),
#       d[, c.sml[c.sml %in% c.binary]]
#     )}

#     ## PLS preprocessing
#     pp.std.ge <- preProcess(X.train.all[,c.ge], method=c('center', 'scale'))
#     pp.std.cn <- preProcess(X.train.all[,c.cn], method=c('center', 'scale'))
#     

#     
#     pp.pls.ge <- plsda(X.train.ge, y.train, ncomp = pls.comp)
#     pp.pls.cn <- plsda(X.train.cn, y.train, ncomp = pls.comp)
#     
#     pred.pls <- function(d, c.sml){ cbind(
#       pls::predict.mvr(pp.pls.ge, d[,c.ge], type = "scores") %>% setNames(paste0('ge.', names(.))),
#       pls::predict.mvr(pp.pls.cn, d[,c.cn], type = "scores") %>% setNames(paste0('cn.', names(.))),
#       d[, c.sml[c.sml %in% c.binary]]
#     )}
#     X.train.pls <- pred.pls(cbind(X.train.ge, X.train.cn), names(X.train.sml))

CreateFoldIndex <- function(y, index, level){
  if (level == 1) # Outer folds should be returned for training set, not test set
    createMultiFolds(y, k = 10, times = 3)
  else if (level == 2) # Inner folds should have a single partition
    createDataPartition(y[index], p=.8)  # Previously: createFolds(y[index], k=10, returnTrain = T)
  else
    stop(sprintf('Folds at level %s not supported', level))
}

GetDataSummarizer <- function(){
  function(d){
    gfl <- function(f, t) length(f[str_detect(f, sprintf('^%s.', t))])
    n.sml <- names(d$X.train.sml)
    n.lrg <- names(d$X.train.lrg)
    n.pca <- names(d$X.train.pca)
    logdebug('Dimension reduction summary after preprocessing:')
    logdebug('Mutation features    : %s --> %s / %s / %s', 
             gfl(d$X.names, 'mu'), gfl(n.lrg, 'mu'), gfl(n.sml, 'mu'), gfl(n.pca, 'mu'))
    logdebug('Expression features  : %s --> %s / %s / %s', 
             gfl(d$X.names, 'ge'), gfl(n.lrg, 'ge'), gfl(n.sml, 'ge'), gfl(n.pca, 'ge'))
    logdebug('Copy Number features : %s --> %s / %s / %s', 
             gfl(d$X.names, 'cn'), gfl(n.lrg, 'cn'), gfl(n.sml, 'cn'), gfl(n.pca, 'cn'))
  }
}

#RESULT_METRICS <- c('auc', 'acc', 'tpr', 'tnr')
ResSummaryFun <- function(curve.type='roc') function(d) GetResultSummary(d, curve.type=curve.type)
GetResultSummary <- function(d, curve.type='lift'){
  
  # Create performance measures
  pred <- prediction(d$y.pred.prob, d$y.test)
  cmat <- confusionMatrix(d$y.pred.class, d$y.test, positive='pos')
  
  cts <- apply(as.data.frame(cmat$table), 1, function(x) c(sprintf('cm.pred.%s.true.%s', x[1], x[2]), x[3])) 
  cts <- t(as.integer(cts[2,]) %>% setNames(cts[1,]))
  
  # Extract desired measures
  auc <- performance(pred, 'auc')
  kappa <- as.numeric(cmat$overall['Kappa'])
  mcp <- coalesce(as.numeric(cmat$overall['McnemarPValue']), 1)
  acc <- as.numeric(cmat$overall['Accuracy'])
  nir <- sum(d$y.test == 'neg') / length(d$y.test)
  cacc <- acc - nir
  bacc <- as.numeric(cmat$byClass['Balanced Accuracy'])
  sens <- as.numeric(cmat$byClass['Sensitivity'])
  spec <- as.numeric(cmat$byClass['Specificity'])
  
  if (tolower(curve.type) == 'roc') {curve <- performance(pred, 'tpr', 'fpr') }
  else if (tolower(curve.type) == 'pr') {curve <- performance(pred, 'prec', 'rec')}
  else if (tolower(curve.type) == 'lift') {curve <- performance(pred, 'lift', 'rpp')}
  else if (tolower(curve.type) == 'gain') {curve <- performance(pred, 'tpr', 'rpp')}
  else {stop(sprintf('Invalid curve type %s', curve.type))}
  
  #   margin.stats <- foreach(margin=seq(0, .4, by=.05), .combine=cbind) %do%{
  #     y.pred.adj <- d$y.pred.prob %>% sapply(function(p){
  #       if (p < .5 - margin) 'neg' else if (p > .5 + margin) 'pos' else NA
  #     }) %>% factor(levels=c('neg', 'pos'))
  #     n.not.na <- sum(!is.na(y.pred.adj))
  #     acc.adj <- ifelse(n.not.na > 0, sum(y.pred.adj == d$y.test & !is.na(y.pred.adj)) / n.not.na, 0)
  #     n.pct <- n.not.na / length(y.pred.adj)
  #     n.pos.adj <- sum(!is.na(y.pred.adj) & d$y.test == 'pos'); 
  #     n.neg.adj <- sum(!is.na(y.pred.adj) & d$y.test == 'neg'); 
  #     tpr.adj <- ifelse(n.pos.adj > 0, sum(y.pred.adj == d$y.test & d$y.test == 'pos', na.rm=T) / n.pos.adj, 0)
  #     tnr.adj <- ifelse(n.neg.adj > 0, sum(y.pred.adj == d$y.test & d$y.test == 'neg', na.rm=T) / n.neg.adj, 0)
  #     res <- data.frame(acc=acc.adj, tpr=tpr.adj, tnr=tnr.adj, n.pct=n.pct)
  #     res %>% setNames(paste(names(res), margin, sep='_margin_'))
  #   }
  
  res <- data.frame(
    x=curve@x.values[[1]], y=curve@y.values[[1]], 
    t=curve@alpha.values[[1]], auc=auc@y.values[[1]],
    acc=acc, bacc=bacc, spec=spec, sens=sens, kappa=kappa, mcp=mcp,
    cacc=cacc, nir=nir, len=length(d$y.test)
  ) %>% cbind(cts)
  #   res <- cbind(res, margin.stats)
  #   if (any(is.na(res[,'acc_margin_0.3'])))
  #     browser()
  res
}

