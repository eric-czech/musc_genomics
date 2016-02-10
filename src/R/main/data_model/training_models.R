#'-----------------------------------------------------------------------------
#' ML Model Specifications
#'
#' This module contains a variety of functions and definitions used by the
#' ML training script to fit and determine model performance.
#' 
#' @author eczech
#'-----------------------------------------------------------------------------

##### Constants #####

bin.tgt.metric <- 'Accuracy'
reg.tgt.metric <- 'RMSE'

##### Utility Functions #####

ShowBestTune <- function(model){ sapply(model, function(m) m$fit$bestTune )}

predict.reg.data.sml <- function(fit, d, i){ predict(fit, d$X.test.sml[,names(d$X.train.sml)]) }
predict.reg.data.lrg <- function(fit, d, i){ predict(fit, d$X.test.lrg[,names(d$X.train.lrg)]) }
predict.reg.ens.sml  <- function(fit, d, i){ predict(fit, newdata=d$X.test.sml[,names(d$X.test.sml)]) }
predict.bin.data.sml <- function(fit, d, i){ predict(fit, d$X.test.sml[,names(d$X.train.sml)], type='prob')[,2] }
predict.bin.data.lrg <- function(fit, d, i){ predict(fit, d$X.test.lrg[,names(d$X.train.lrg)], type='prob')[,2] }
predict.bin.ens.sml  <- function(fit, d, i){ predict(fit, newdata=d$X.test.sml[,names(d$X.test.sml)]) }
predict.browser <- function(fit, d, i){ browser() }

test.bin <- function(d) d$y.test.bin
test.reg <- function(d) d$y.test

bin.trctrl <- function(index, ...) trainControl(
  summaryFunction=function(...) c(twoClassSummary(...), defaultSummary(...)),
  method = "cv", number = 5, savePredictions='final', index=index, returnData=F,
  ...)

reg.trctrl <- function(index, ...) trainControl(
  summaryFunction=function(...) defaultSummary(...),
  method = "cv", number = 5, savePredictions='final', index=index, returnData=F,
  ...)


GetPLSModel <- function(){
  m <- getModelInfo(model = "pls", regex = FALSE)[[1]]
  m$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {   
    out <- if(is.factor(y)){
      plsda(x, y, method = "oscorespls", ncomp = param$ncomp, model=F, ...)
    } else {
      dat <- if(is.data.frame(x)) x else as.data.frame(x)
      dat$.outcome <- y
      plsr(.outcome ~ ., data = dat, method = "oscorespls", ncomp = param$ncomp, ...)
    }
    out
  }
  m
}

GetRDAModel <- function(){
  m <- getModelInfo(model = "rda", regex = FALSE)[[1]]
  m$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {  
    klaR:::rda(x, y, gamma = param$gamma, lambda = param$lambda, fold=3, ...)
  }
  m
}

GetSCRDAModel <- function(){
  m <- list(type = "Classification", library = "rda", loop = NULL)
  m$parameters <- data.frame(parameter = c("C", "sigma"),
                             class = rep("numeric", 2),
                             label = c("Cost", "Sigma"))
  rda::rda(colon.x[, tr.index], colon.y[tr.index])
}

GetElasticNetModel <- function(name, alpha.start, alpha.grid, n.core){
  list(
    name=name, predict=predict.bin.data.sml, test=test.bin,
    train=function(d, idx, ...){
      registerDoMC(n.core)
      glmnet.lambda <- GetGlmnetLambda(d$X.train.sml, d$y.train.bin, alpha=alpha.start, family='binomial')
      train(
        d$X.train.sml, d$y.train.bin, 
        method='glmnet', preProcess='zv', family='binomial', metric=bin.tgt.metric, 
        tuneGrid = expand.grid(.alpha = alpha.grid, .lambda=glmnet.lambda),
        trControl = bin.trctrl(idx, classProbs=T)
      )
    }
  )
}

GetEnsembleModel <- function(models, name, test.selector, pred.fun){
  list(
    name=name, test=test.selector,
    train=function(d, idx, i, ...){
      m <- lapply(models, function(m) m(i))
      class(m) <- "caretList"
      caretEnsemble(m)
    }, predict=function(fit, d, i){ 
      pred.fun(fit, d, i)
    }
  )
}


##### Classification Models #####


bin.model.lasso <- GetElasticNetModel('bin.lasso', 1, 1, 1)
bin.model.ridge <- GetElasticNetModel('bin.ridge', 0, 0, 1)
bin.model.enet <- GetElasticNetModel('bin.enet', .5, c(.001, seq(.1, .9, length.out=13), .999), 1)

bin.model.svm.radial.sml <- list(
  name='bin.svm.radial.sml', predict=predict.bin.data.sml, test=test.bin,
  train=function(d, idx, ...){
    registerDoMC(6)
    train(
      d$X.train.sml, d$y.train.bin,
      method='svmRadial', preProcess='zv', metric=bin.tgt.metric,
      tuneLength=15, trControl = bin.trctrl(idx, classProbs=T)
    )
  }
)

bin.model.pls <- list(
  name='bin.pls', predict=predict.bin.data.sml, test=test.bin,
  train=function(d, idx, ...){
    registerDoMC(3)
    train(
      d$X.train.sml, d$y.train.bin, 
      method=GetPLSModel(), preProcess='zv', metric=bin.tgt.metric, 
      trControl = bin.trctrl(idx, classProbs=T),
      tuneGrid=data.frame(ncomp=c(1,2,3,6,9,12,15,18,25))
    )
  }
)

bin.model.rda <- list(
  name='bin.rda', predict=predict.bin.data.sml, test=test.bin,
  train=function(d, idx, ...){
    registerDoMC(3)
    train(
      d$X.train.sml, d$y.train.bin, 
      method=GetRDAModel(), preProcess='zv', metric=bin.tgt.metric, 
      trControl=bin.trctrl(idx, classProbs=T),
      tuneGrid=expand.grid(.lambda=c(.1, .5, .9), .gamma=c(.1, .5, .9))
    )
  }
)

bin.model.knn.pca <- list(
  name='bin.knn.pca', predict=predict.bin.data.sml, test=test.bin,
  train=function(d, idx, ...){
    registerDoMC(3)
    train(
      d$X.train.sml, d$y.train.bin, 
      method='knn', preProcess=c('zv', 'pca'), metric=bin.tgt.metric, 
      tuneLength=15, trControl = bin.trctrl(idx, classProbs=T, preProcOptions = list(thresh = .9))
    )
  }
)


bin.model.knn <- list(
  name='bin.knn', predict=predict.bin.data.sml, 
  test=test.bin,
  train=function(d, idx, ...){
    registerDoMC(1)
    train(
      d$X.train.sml, d$y.train.bin, 
      method='knn', preProcess='zv', metric=bin.tgt.metric, 
      tuneLength=15, trControl = bin.trctrl(idx, classProbs=T)
    )
  }
)

bin.model.pam <- list(
  name='bin.pam', predict=predict.bin.data.sml, test=test.bin,
  train=function(d, idx, ...){
    registerDoMC(1)
    train(
      d$X.train.sml, d$y.train.bin, 
      method='pam', preProcess='zv', metric=bin.tgt.metric, 
      tuneGrid=data.frame(.threshold=c(0, 10^seq(-8, 1, length.out=10))),
      trControl = bin.trctrl(idx, classProbs=T)
    )
  } 
)

# bin.model.nb <- list(
#   name='bin.nb', predict=predict.bin.data,
#   train=function(d, idx, ...){
#     registerDoMC(1)
#     train(
#       d$X.train.sml, setLevels(d$y.train, c('neg', 'pos')), method='nb', 
#       metric=bin.tgt.metric, tuneLength=2, trControl = bin.trctrl(idx, classProbs=T)
#     )
#   }
# )

bin.model.rf <- list(
  name='bin.rf', predict=predict.bin.data.sml, test=test.bin,
  train=function(d, idx, ...){
    registerDoMC(8)
    train(
      d$X.train.sml, d$y.train.bin, 
      method='rf', preProcess='zv', metric=bin.tgt.metric, 
      tuneLength=5,
      #tuneGrid=data.frame(.mtry=c())
      trControl = bin.trctrl(idx, classProbs=T, verboseIter=T)
    )
  }
)


bin.model.gbm <- list(
  name='bin.gbm', predict=predict.bin.data.sml, test=test.bin,
  train=function(d, idx, ...){
    registerDoMC(5)
    train(
      d$X.train.sml, d$y.train.bin, 
      method='gbm', preProcess='zv', bag.fraction=.3, 
      metric=bin.tgt.metric, tuneLength=5,
      trControl = bin.trctrl(idx, classProbs=T, verboseIter=T)
    )
  }
)

bin.model.et <- list(
  name='bin.et', predict=predict.bin.data.sml, test=test.bin,
  train=function(d, idx, ...){
    registerDoMC(4)
    train(
      d$X.train.sml, d$y.train.bin, 
      method='extraTrees', preProcess='zv', metric=bin.tgt.metric, 
      tuneLength=5, trControl = bin.trctrl(idx, classProbs=T, verboseIter=T)
    )
  }
)


##### Regression Models #####

reg.model.svm.radial <- list(
  name='reg.svm.radial', predict=predict.reg.data.sml, test=test.reg,
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
  name='reg.pls', predict=predict.reg.data.sml, test=test.reg,
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