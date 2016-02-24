#'-----------------------------------------------------------------------------
#' ML Model Specifications
#'
#' This module contains a variety of functions and definitions used by the
#' ML training script to fit and determine model performance.
#' 
#' @author eczech
#'-----------------------------------------------------------------------------
source('data_model/custom_models.R')

##### Constants #####

# bin.tgt.metric <- 'Accuracy' 
bin.tgt.metric <- 'Kappa' # Switching to Kappa due to imbalance
reg.tgt.metric <- 'RMSE'

##### Utility Functions #####

ShowBestTune <- function(model){ sapply(model, function(m) m$fit$bestTune )}

reg.predict.sml <- function(fit, d, i){ predict(fit, d$X.test.sml[,names(d$X.train.sml)]) }
reg.predict.lrg <- function(fit, d, i){ predict(fit, d$X.test.lrg[,names(d$X.train.lrg)]) }
reg.predict.ens.sml  <- function(fit, d, i){ predict(fit, newdata=d$X.test.sml[,names(d$X.test.sml)]) }

bin.predict.sml <- function(fit, d, i){ list(
  prob=predict(fit, d$X.test.sml[,names(d$X.train.sml)], type='prob')[,1], 
  class=predict(fit, d$X.test.sml[,names(d$X.train.sml)], type='raw')
)}
bin.predict.rfe.sml <- function(fit, d, i){ 
  preds <- predict(fit, d$X.test.sml[,names(d$X.train.sml)])
  list(prob=preds[,2], class=preds[,1])
}
bin.predict.sml.ge <- function(fit, d, i){ 
  cols <- GetFeatures(d$X.train.sml, 'ge')
  list(
    prob=predict(fit, d$X.test.sml[,cols], type='prob')[,1], 
    class=predict(fit, d$X.test.sml[,cols], type='raw')
  )
}
bin.predict.pca <- function(fit, d, i){ list(
  prob=predict(fit, d$X.test.pca[,names(d$X.train.pca)], type='prob')[,1], 
  class=predict(fit, d$X.test.pca[,names(d$X.train.pca)], type='raw')
)}
bin.predict.lrg <- function(fit, d, i){ list(
  prob=predict(fit, d$X.test.lrg[,names(d$X.train.lrg)], type='prob')[,1],
  class=predict(fit, d$X.test.lrg[,names(d$X.train.lrg)], type='raw')
)}
bin.predict.ens.sml  <- function(fit, d, i){ list(
  prob=predict(fit, newdata=d$X.test.sml[,names(d$X.test.sml)], type='prob'),
  class=predict(fit, newdata=d$X.test.sml[,names(d$X.test.sml)], type='raw')
)}

predict.browser <- function(fit, d, i){ browser() }

bin.test <- function(d) d$y.test.bin
reg.test <- function(d) d$y.test
bin.train.sml <- function(d) d$X.train.sml

feature.selector.ge <- function(d) d[, GetFeatures(d, 'ge')]
feature.selector.cn <- function(d) d[, GetFeatures(d, 'cn')]
feature.selector.mu <- function(d) d[, GetFeatures(d, 'mu')]
bin.train.sml.ge <- function(d) feature.selector.ge(d$X.train.sml)
bin.train.sml.cn <- function(d) feature.selector.cn(d$X.train.sml)
bin.train.sml.mu <- function(d) feature.selector.mu(d$X.train.sml)

bin.train.lrg <- function(d) d$X.train.lrg
bin.train.pca <- function(d) d$X.train.pca

get.weight.fun <- function(positive.weight){
  function(y) ifelse(y == 'pos', positive.weight, 1)
}
bin.weights <- get.weight.fun(5)


# Note that this is called by caret with non-predicted data using an arbitrary
# sample of the training data (10 rows) to get performance measure names 1 time per train
ClassSummary <- function(data, ...) {
  if (!is.factor(data$obs)) 
    stop("Outcome is not a factor (ClassSummary does not apply)")
  if (length(levels(data$obs)) != 2)
    stop("Outcome must have two levels (ClassSummary does not apply)")
  if (!all(levels(data$pred) == levels(data$obs))) 
    stop("levels of observed and predicted data do not match")
  # res <- twoClassSummary(data, ...)['ROC']
  cmat <- confusionMatrix(data$pred, data$obs, positive = 'pos')
  res <- c(
    coalesce(cmat$overall['Kappa'], -1),
    cmat$overall['Accuracy'],
    coalesce(cmat$byClass['Sensitivity'], 0),
    coalesce(cmat$byClass['Specificity'], 0),
    Precision=coalesce(as.numeric(cmat$byClass[c('Pos Pred Value')]), 0)
  )
  res
}

bin.trctrl <- function(index, returnData=F, ...) trainControl(
  method='cv', summaryFunction=ClassSummary,
  savePredictions='final', index=index, returnData=returnData,
  ...)

reg.trctrl <- function(index, ...) trainControl(
  summaryFunction=function(...) defaultSummary(...),
  method = "cv", number = 5, savePredictions='final', index=index, returnData=F,
  ...)


bin.model <- function(name, n.core, train.fun, pred.fun, weight.fun=NULL, ...){
  list(
    name=name, predict=pred.fun, test=bin.test,
    train=function(d, idx, i){
      registerDoMC(n.core)
      weights <- NULL
      if (!is.null(weight.fun))
        weights <- weight.fun(d$y.train.bin)
      train(
        train.fun(d), d$y.train.bin, ...,
        metric=bin.tgt.metric, weights=weights,
        trControl = bin.trctrl(idx, classProbs=T)
      )
    }
  )
}

bin.rfe.size.fun <- function(ncol){
  c(1, 5, 10, 25, 50, 100)
}
rfeFuncs <- caretFuncs
rfeFuncs$rank <- function (object, x, y) {
  vimp <- varImp(object, scale = FALSE)$importance
  if (object$modelType == "Regression") {
    vimp <- vimp[order(vimp[, 1], decreasing = TRUE), , drop = FALSE]
  }
  else {
    if (all(levels(y) %in% colnames(vimp))) {
      avImp <- apply(vimp[, levels(y), drop = TRUE], 1, 
                     mean)
      vimp$Overall <- avImp
    }
  }
  vimp$var <- rownames(vimp)
  vimp
}
bin.rfe.model <- function(name, n.core, train.fun, pred.fun, sizes.fun=bin.rfe.size.fun, weight.fun=NULL, ...){
  list(
    name=name, predict=pred.fun, test=bin.test,
    train=function(d, idx, i){
      registerDoMC(n.core)
      sizes = sizes.fun(ncol(train.fun(d)))
      weights <- NULL
      if (!is.null(weight.fun))
        weights <- weight.fun(d$y.train.bin)
      rfectrl <- rfeControl(
        functions=rfeFuncs, index=idx,
        returnResamp="final", verbose = T, allowParallel = F
      )
      rfe(
        train.fun(d), d$y.train.bin, ...,
        rfeControl=rfectrl, metric=bin.tgt.metric, 
        weights=weights, sizes=sizes, 
        trControl=bin.trctrl(NULL, classProbs=T, returnData=T)
      )
    }
  )
}

##### Model Wrapper Functions #####

GetElasticNetModel <- function(name, alpha.start, alpha.grid, n.core, train.fun, pred.fun, weight.fun=NULL){
  list(
    name=name, predict=pred.fun, test=bin.test,
    train=function(d, idx, ...){
      registerDoMC(n.core)
      
      glmnet.lambda <- GetGlmnetLambda(train.fun(d), d$y.train.bin, alpha=alpha.start, family='binomial')
      weights <- NULL
      if (!is.null(weight.fun))
        weights <- weight.fun(d$y.train.bin)
      train(
        train.fun(d), d$y.train.bin, 
        method='glmnet', preProcess='zv', family='binomial', metric=bin.tgt.metric, 
        weights=weights, tuneGrid = expand.grid(.alpha = alpha.grid, .lambda=glmnet.lambda),
        trControl = bin.trctrl(idx, classProbs=T)
      )
    }
  )
}

# GetDataSubsetModel
GetPartitionedEnsembleModel <- function(name, model.ge, model.cn, model.mu, test.selector,
                                        method='glm', tuneLength=1){
  list(
    name=name, test=test.selector,
    train=function(d, idx, i, ...){
      m <- list()
      # Num core selection should happen in individual models
      m$ge <- model.ge$train(d, idx, i, ...)
      m$cn <- model.cn$train(d, idx, i, ...)
      m$mu <- model.mu$train(d, idx, i, ...)
      trControl <- trainControl(
        method='cv', summaryFunction=ClassSummary,
        savePredictions='final'
      )
      class(m) <- "caretList"
      caretStack(m, method=method, tuneLength=tuneLength, trControl=trControl)
    }, predict=function(fit, d, i){ 
      list(
        prob=predict(fit, newdata=d$X.test.sml[,names(d$X.test.sml)], type='prob'),
        class=predict(fit, newdata=d$X.test.sml[,names(d$X.test.sml)], type='raw')
      )
    }
  )
}

GetEnsembleModel <- function(models, name, test.selector, pred.fun, method, ...){
  args <- list(...)
  list(
    name=name, test=test.selector,
    train=function(d, idx, i, ...){
      registerDoMC(1)
      m <- lapply(models, function(m) m(i))
      class(m) <- "caretList"
      trControl <- trainControl(
        method='cv', summaryFunction=ClassSummary,
        savePredictions='final'
      )
      a <- args
      a$trControl <- trControl
      a$method <- method
      a$all.models <- m
      do.call(function(...) { caretStack(...)}, rev(a))
    }, predict=function(fit, d, i){ 
      pred.fun(fit, d, i)
    }
  )
}

##### RFE Models #####



##### Classification Models #####


### ENET ###
alpha.grid <- c(.001, seq(.1, .9, length.out=13), .999)
bin.model.lasso.sml <- GetElasticNetModel('lasso.sml', 1, 1, 5, bin.train.sml, bin.predict.sml)
bin.model.lasso.wt.sml <- GetElasticNetModel('lasso.wt.sml', 1, 1, 5, bin.train.sml, bin.predict.sml, bin.weights)
bin.model.ridge.sml <- GetElasticNetModel('ridge.sml', 0, 0, 5, bin.train.sml, bin.predict.sml)
bin.model.ridge.wt.sml <- GetElasticNetModel('ridge.wt.sml', 0, 0, 5, bin.train.sml, bin.predict.sml, bin.weights)
bin.model.enet.sml <- GetElasticNetModel('enet.sml', .5, alpha.grid, 5, bin.train.sml, bin.predict.sml)
bin.model.enet.wt.sml <- GetElasticNetModel('enet.wt.sml', .5, alpha.grid, 5, bin.train.sml, bin.predict.sml, bin.weights)

bin.model.lasso.pca <- GetElasticNetModel('lasso.pca', 1, 1, 5, bin.train.pca, bin.predict.pca)
bin.model.ridge.pca <- GetElasticNetModel('ridge.pca', 0, 0, 5, bin.train.pca, bin.predict.pca)
bin.model.enet.pca <- GetElasticNetModel('enet.pca', .5, alpha.grid, 5, bin.train.pca, bin.predict.pca)

### SCRDA ###
bin.model.scrda.sml <- bin.model(
  'scrda.sml', 5, bin.train.sml, bin.predict.sml, 
  method=GetSCRDAModel(10), preProcess='zv', tuneLength=15
)

bin.model.scrda.lrg <- bin.model(
  'scrda.lrg', 3, bin.train.lrg, bin.predict.lrg, 
  method=GetSCRDAModel(3), preProcess='zv', tuneLength=10
)

test.rfe.size.fun <- function(ncol){
  c(10, 25, 50, 75, 100, 200)
}
bin.model.scrda.rfe.sml <- bin.rfe.model(
  'scrda.rfe', 10, bin.train.sml, bin.predict.rfe.sml, sizes.fun = test.rfe.size.fun,
  method=GetSCRDAModel(8, var.imp=F), tuneLength=8
)

## HDRDA ###
hdrda.grid <- rbind(
  expand.grid(lambda=seq(0, 1, len = 6), gamma=c(0, 10^seq.int(-4, 4, len = 7)), shrinkage='ridge', stringsAsFactors=F),
  expand.grid(lambda=seq(0, 1, len = 6), gamma=seq(0, 1, len = 8), shrinkage='convex', stringsAsFactors=F)
)
bin.model.hdrda.sml <- bin.model(
  'hdrda.sml', 10, bin.train.sml, bin.predict.sml, 
  method=GetHDRDAModel(), preProcess='zv', tuneGrid=hdrda.grid#tuneLength=10
)

bin.model.hdrda.lrg <- bin.model(
  'hdrda.lrg', 3, bin.train.lrg, bin.predict.lrg, 
  method=GetHDRDAModel(), preProcess='zv', tuneLength=5
)



### RDA ###
bin.model.rda.sml <- bin.model(
  'rda.sml', 1, bin.train.sml, bin.predict.sml, 
  method=GetRDAModel(), preProcess='zv', tuneLength=1
)

bin.model.slda.sml <- bin.model(
  'slda.sml', 1, bin.train.sml.ge, bin.predict.sml.ge, 
  method='PenalizedLDA', preProcess='zv', tuneLength=1
)




### SVM ###
bin.model.svm.radial.sml <- bin.model(
  'svm.radial.sml', 6, bin.train.sml, bin.predict.sml,
  method='svmRadial', preProcess='zv', tuneLength=25
)
bin.model.svm.radial.pca <- bin.model(
  'svm.radial.pca', 6, bin.train.pca, bin.predict.pca, 
  method='svmRadial', preProcess='zv', tuneLength=25
)
bin.model.svm.linear.sml <- bin.model(
  'svm.linear.sml', 6, bin.train.sml, bin.predict.sml,
  method='svmLinear', preProcess='zv', tuneLength=35
)
bin.model.svm.linear.pca <- bin.model(
  'svm.linear.pca', 6, bin.train.pca, bin.predict.pca, 
  method='svmLinear', preProcess='zv', tuneLength=35
)

bin.model.svm.rfe.wt.sml <- bin.rfe.model(
  'svm.wt.rfe', 2, bin.train.sml, bin.predict.sml, weight.fun=bin.weights,
  method='svmRadialWeights', preProcess='zv', tuneLength=15
)
  
  
### PLS ###
bin.model.pls.sml <- bin.model(
  'pls.sml', 3, bin.train.sml, bin.predict.sml, 
  method=GetPLSModel(), preProcess='zv', 
  tuneGrid=data.frame(ncomp=c(1,2,3,4,5,6,9,12,15,18,25))
)
bin.model.pls.pca <- bin.model(
  'pls.pca', 3, bin.train.pca, bin.predict.pca,
  method=GetPLSModel(), preProcess='zv', 
  tuneGrid=data.frame(ncomp=c(1,2,3,4,5,6,9,12,15,18,25))
)

### MARS ###
bin.model.mars.sml <- bin.model(
  'mars.sml', 1, bin.train.sml, bin.predict.sml, 
  method='earth', preProcess='zv', tuneLength=10
)
bin.model.mars.pca <- bin.model(
  'mars.pca', 3, bin.train.pca, bin.predict.pca,
  method='earth', preProcess='zv', tuneLength=10
)

### RDA ###
bin.model.rda.sml <- bin.model(
  'rda.sml', 3, bin.train.sml, bin.predict.sml, 
  method=GetRDAModel(), preProcess='zv', 
  tuneGrid=expand.grid(.lambda=c(.1, .5, .9), .gamma=c(.1, .5, .9))
)
bin.model.rda.pca <- bin.model(
  'rda.pca', 3, bin.train.pca, bin.predict.pca,
  method=GetRDAModel(), preProcess='zv', 
  tuneGrid=expand.grid(.lambda=c(.1, .5, .9), .gamma=c(.1, .5, .9))
)

### KNN ###
bin.model.knn.sml <- bin.model(
  'knn.sml', 5, bin.train.sml, bin.predict.sml, 
  method='knn', preProcess='zv', tuneLength=15
)
bin.model.knn.pca <- bin.model(
  'knn.pca', 5, bin.train.pca, bin.predict.pca,
  method='knn', preProcess='zv', tuneLength=15
)
test.knn.size.fun <- function(ncol){
  c(10, 50, 100)
}
bin.model.knn.rfe.sml <- bin.rfe.model(
  'knn.rfe', 5, bin.train.sml, bin.predict.rfe.sml, sizes.fun = test.knn.size.fun,
  method='knn', tuneLength=2
)


### PAM ###
bin.model.pam.sml <- bin.model(
  'pam.sml', 5, bin.train.sml, bin.predict.sml, 
  method='pam', preProcess='zv', 
  # See https://github.com/topepo/caret/blob/master/models/files/pam.R for details on how 
  # grid is chosen by running pamr.train which returns a vector of threshold values determined
  # by the package.  By default, the package returns 30 values and caret rescales those values
  # to have the same range but of length 'tuneLength'.  Using 30 for tuneLength maintains 
  # consistency with choice by package, but removes first and last values (not sure why)
  tuneLength=30
)
bin.model.pam.pca <- bin.model(
  'pam.pca', 5, bin.train.pca, bin.predict.pca,
  method='pam', preProcess='zv', tuneLength=30
)

### RF ###
bin.model.rf.sml <- bin.model(
  'rf.sml', 5, bin.train.sml, bin.predict.sml, 
  method='rf', preProcess='zv', tuneLength=4
)
bin.model.rf.pca <- bin.model(
  'rf.pca', 5, bin.train.pca, bin.predict.pca,
  method='rf', preProcess='zv', tuneLength=8
)
# bin.model.rf.pca <-   list(
#   name='rf.pca', predict=bin.train.pca, test=bin.test,
#   train=function(d, idx, i){
#     registerDoMC(n.core)
#     train(
#       bin.train.pca(d), d$y.train.bin, method='rf', preProcess='zv', tuneLength=8,
#       metric=bin.tgt.metric,
#       trControl = bin.trctrl(idx, classProbs=T)
#     )
#   }
# )

### GBM ###
bin.model.gbm.sml <- bin.model(
  'gbm.sml', 5, bin.train.sml, bin.predict.sml, 
  method='gbm', preProcess='zv', bag.fraction=.3, tuneLength=5, verbose=F
)
bin.gbm.grid <- expand.grid(
  interaction.depth = c(1,2,3),
  n.trees=c(50, 250),
  shrinkage=10^(seq(-5, -1, len=4)),
  n.minobsinnode=c(5, 12)
)
bin.model.gbm.wt.sml <- bin.model(
  'gbm.wt.sml', 5, bin.train.sml, bin.predict.sml, weight.fun=bin.weights,
  method='gbm', preProcess='zv', tuneGrid=bin.gbm.grid, verbose=F, bag.fraction=.3
)
bin.model.gbm.pca <- bin.model(
  'gbm.pca', 5, bin.train.pca, bin.predict.pca,
  method='gbm', preProcess='zv', bag.fraction=.3, tuneLength=5, verbose=F
)

### C5.0 ###

bin.model.c50.wt.sml <- bin.model(
  'c50.wt.sml', 5, bin.train.sml, bin.predict.sml, weight.fun=bin.weights,
  method='C5.0', preProcess='zv', tuneLength=2
)

### ET ###
bin.model.et.sml <- bin.model(
  'et.sml', 1, bin.train.sml, bin.predict.sml, 
  method='extraTrees', preProcess='zv', tuneLength=3
)
bin.model.et.pca <- bin.model(
  'et.pca', 1, # Do not increase or no progess ever seems to be made 
  bin.train.pca, bin.predict.pca,
  method='extraTrees', preProcess='zv', tuneLength=3
)




##### Regression Models #####

reg.model.svm.radial <- list(
  name='reg.svm.radial', predict=reg.predict.sml, test=reg.test,
  train=function(d, idx, ...){
    registerDoMC(3)
    train(
      d$X.train.sml, d$y.train, 
      method='pam', preProcess='zv', metric=reg.tgt.metric, 
      tuneLength=15, search='grid',
      trControl = reg.trctrl(idx)
    )
  }
)

reg.model.pls <- list(
  name='reg.pls', predict=reg.predict.sml, test=reg.test,
  train=function(d, idx, ...){
    registerDoMC(1)
    train(
      d$X.train.sml, d$y.train, 
      method='pls', preProcess='zv', metric=reg.tgt.metric, 
      tuneGrid=data.frame(ncomp=c(1,2,3,6,9,12,15,18)),
      trControl = reg.trctrl(idx)
    )
  }
)