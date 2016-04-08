

##### Constants #####
BIAS_FEATS <- c('bias1', 'bias2')

##### Data Prep Functions #####
GetPrepFun <- function(origin.transform=NULL){
  function(X){
    if (!'origin' %in% names(X)) X
    else {
      if (is.null(origin.transform))
        stop('Origin found in feature set but transform given was null')
      lvls <- as.character(0:(origin.transform$max.val))
      X <- X %>% mutate(origin=factor(origin.transform$convert(origin), levels=lvls))
      as.data.frame(model.matrix(~. - 1, X))
    }
  }
}


GetRandVars <- function(n){
  set.seed(SEED)
  data.frame(rnorm(n), rnorm(n)) %>% setNames(BIAS_FEATS)
}

GetFilterModel <- function(name, max.feats, n.core=3, 
                           origin.transform=NULL, weight.fun=NULL, ...){
  prep.fun <- GetPrepFun(origin.transform)
  pred.fun <- function(fit, d, i){ 
    X.test <- d$X.test
    if (length(fit$top.feats) == length(BIAS_FEATS) && all(fit$top.feats == BIAS_FEATS)){
      X.test <- GetRandVars(nrow(X.test))
    }
    X <- X.test[, fit$top.feats, drop=F] %>% prep.fun
    list(
      prob=predict(fit, X, type='prob')[,1], 
      class=predict(fit, X, type='raw')
    )
  }
  list(
    name=name, predict=pred.fun, test=bin.test,
    train=function(d, idx, i){
      registerDoMC(n.core)
      # d$feat.scores contains scores for all 36k features
      feats <- d$feat.scores
      if (is.null(origin.transform))
        feats <- feats %>% filter(feature != 'origin')
      
      # Temporarily filter to only GE and MU features
      # feats <- feats %>% filter(str_detect(feature, '^mu\\.') | str_detect(feature, '^ge\\.'))
      
      X.train <- d$X.train
      if (max.feats == 0){
        top.feats <- BIAS_FEATS
        X.train <- GetRandVars(nrow(X.train))
      } else {
        top.feats <- feats %>% arrange(score) %>% head(max.feats) %>% .$feature 
      }
      X <- prep.fun(X.train[,top.feats, drop=F])
      y <- d$y.train.bin
      
      weights <- NULL
      if (!is.null(weight.fun))
        weights <- weight.fun(y)
      
      fit <- train(X, y, weights=weights, ...)
      fit$top.feats <- top.feats
      trim_model(fit) # Remove massive, embedded '.Environment' attributes
    }
  )
}

GetModelForTransform <- function(max.feats, model.name, n.core=3, k=5, allow.parallel=T, 
                                 origin.transform=NULL, origin.name=NULL, weight.fun=NULL, ...){
  if (is.null(origin.name))
    origin.name <- ifelse(is.null(origin.transform), 'norigin', 'worigin')
  name <- sprintf('%s.%s.%s', model.name, max.feats, origin.name)
  GetFilterModel(
    name, max.feats, n.core=n.core,
    origin.transform=origin.transform,
    weight.fun=weight.fun,
    metric='Accuracy', 
    ...,
    trControl=trainControl(
      method='cv', number=k, classProbs=T, 
      returnData=F, savePredictions='final',
      allowParallel=allow.parallel, verboseIter=F
    )
  )
}

GetEnsembleModelForTransform <- function(max.feats, sub.models, model.name, n.core=1, method='glm', 
                                         origin.transform=NULL, origin.name=NULL, ...){
  if (is.null(origin.name))
    origin.name <- ifelse(is.null(origin.transform), 'norigin', 'worigin')
  name <- sprintf('%s.%s.%s', model.name, max.feats, origin.name)
  prep.fun <- GetPrepFun(origin.transform)
  pred.fun <- function(fit, d, i){ 
    X.test <- d$X.test
    if (length(fit$top.feats) == length(BIAS_FEATS) && all(fit$top.feats == BIAS_FEATS)){
      X.test <- GetRandVars(nrow(X.test))
    }
    X <- X.test[, fit$top.feats, drop=F] %>% prep.fun
    tryCatch({
      list(
        prob=predict(fit, X, type='prob'), 
        class=predict(fit, X, type='raw')
      )
    }, error=function(e) browser())
  }
  fit.prep.fun <- function(fit, model.results){
    fit$top.feats <- model.results[[1]]$top.feats
    fit
  }
  res <- c(
    list(name=name, predict=pred.fun, test=bin.test),
    GetFitEnsembleTrain(
      sub.models, fit.prep.fun=fit.prep.fun, n.core=n.core,
      method=method, metric='Accuracy', trControl=trainControl(
        method='none', number=1, 
        classProbs=T, savePredictions='final', returnData=T
      )
    ) 
  )
  res
}



GetFilteringSVMModel <- function(){
  m <- Decorate.Classifier.NAPredictions('svmRadial', c('pos', 'neg'))
  m$method <- 'svmRadial'
  m
}

GetFitEnsembleTrain <- function(model.results, fit.prep.fun=NULL, n.core=1, ...){
  list(
    train=function(d, idx, i){
      registerDoMC(n.core)
      models <- lapply(model.results, function(m){
        if ('fit' %in% names(m)) m$fit
        else m[[i]]$fit
      }) %>% setNames(names(model.results))
      class(models) <- "caretList"
      fit <- caretStack(models, ...)
      if (!is.null(fit.prep.fun))
        fit <- fit.prep.fun(fit, models)
      fit
    }
  )
}


.GetAvgFilterEnsemble <- function(){
  caret.list.args <- list(
    trControl=trainControl(
      method='cv', number=3, classProbs=T, 
      returnData=F, savePredictions='final',
      allowParallel=T, verboseIter=F
    ),
    tuneList=list(
      rf=caretModelSpec(method='rf', tuneLength=4),
      svm=caretModelSpec(method='svmRadial', tuneLength=15),
      enet=caretModelSpec(method='glmnet', tuneLength=15)
    )
  )
  caret.stack.args <- list(
    method=GetEnsembleAveragingModel(),
    trControl=trainControl(method='none', savePredictions = 'final')
  )
  
  GetCaretEnsembleModel(caret.list.args, caret.stack.args)
}


m.ens.avg <- .GetAvgFilterEnsemble()
