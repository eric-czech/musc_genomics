
GetDataSubsetModel <- function(caretModel, subset.selector){
  m <- caretModel
  m$fit <- function(x, y, wts, param, lev, last, classProbs, ...){
    caretModel$fit(subset.selector(x), y, wts, param, lev, last, classProbs, ...)  
  }
  m$predict <- function(modelFit, newdata, submodels = NULL){
    caretModel$predict(modelFit, subset.selector(newdata), submodels)  
  }
  m$prob <- function(modelFit, newdata, submodels = NULL){
    caretModel$prob(modelFit, subset.selector(newdata), submodels)  
  }
  m
}

GetHDRDAModel <- function(){
  list(
    method="hdrda",
    label = "High Dimensional Regularized Discriminant Analysis",
    library = "sparsediscrim",
    loop = NULL,
    type = c('Classification'),
    parameters = data.frame(
      parameter = c('lambda', 'gamma', 'shrinkage'),
      class = c('numeric', 'numeric', 'character'),
      label = c('Lambda', 'Gamma', 'Shrinkage Type')
    ),
    grid = function(x, y, len = NULL, search = "grid") {
      # See recommended hyperparameter settings at:
      # https://github.com/ramhiser/sparsediscrim/blob/master/R/hdrda.r#L315
      if(search == "grid") {
        lambda <- seq(0, 1, len = len)
        gamma_ridge <- c(0, 10^seq.int(-2, 4, len = len-1))
        gamma_convex <- seq(0, 1, len = len)
      } else {
        lambda <- runif(len, min = 0, max = 1)
        gamma_ridge <- runif(len, min = 10^-2, max = 10^4)
        gamma_convex <- runif(len, min = 0, max = 1)
      }
      out <- rbind(
        expand.grid(lambda=lambda, gamma=gamma_ridge, shrinkage='ridge'),
        expand.grid(lambda=lambda, gamma=gamma_convex, shrinkage='convex')
      )
      out$shrinkage <- as.character(out$shrinkage)
      out
    },
    fit = function(x, y, wts, param, lev, last, classProbs, ...) {
      sparsediscrim::hdrda(x, y, lambda=param$lambda, gamma=param$gamma, shrinkage_type=param$shrinkage)
    },
    predict = function(modelFit, newdata, submodels = NULL) {
      r <- predict(modelFit, newdata)
      r$class
    },
    prob = function(modelFit, newdata, submodels = NULL) {
      r <- predict(modelFit, newdata)
      r$posterior
    },
    levels = function(x) x$obsLevels,
    tags = c("Discriminant Analysis", "Linear Classifier"),
    sort = function(x) x[order(x$lambda, x$gamma, x$shrinkage),]
  )
}


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
  m$method <- 'pls-custom'
  m
}

GetRDAModel <- function(){
  m <- getModelInfo(model = "rda", regex = FALSE)[[1]]
  m$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {  
    klaR:::rda(x, y, gamma = param$gamma, lambda = param$lambda, fold=3, ...)
  }
  m$method <- 'rda-custom'
  m
}

GetSCRDAModel <- function(max.delta=3, var.imp=T){
  m <- list(type = "Classification", library = "rda", loop = NULL)
  params <- data.frame(
    parameter=c('alpha', 'delta', 'len'),
    class=c('numeric', 'numeric', 'numeric'),
    label=c('#Alpha', '#Delta', '#Params')
  )
  get.predictions <- function(modelFit, newdata, submodels, type, transform){
    m <- modelFit
    x <- newdata
    
    if (!is.matrix(x))
      x <- as.matrix(x)
    x <- t(x)
    
    params <- rbind(modelFit$loop.params, submodels)
    alpha <- sort(unique(params$alpha))
    delta <- sort(unique(params$delta))
    pred <- rda::predict.rda(m, m$x.train, m$y.train, x, alpha=alpha, delta=delta, type=type)
    
    if (is.null(submodels))
      return(transform(pred))
    
    stopifnot(length(alpha) == dim(pred)[1])
    stopifnot(length(delta) == dim(pred)[2])
    
    foreach(i=1:length(alpha), .combine=c)%:%foreach(j=1:length(delta))%do%{
      if (type == 'posterior') transform(pred[i, j, , ]) else transform(pred[i, j, ])
    }
  }
  m <- list(
    method="scrda",
    label = "SCRDA",
    library = c("rda"),
    type = "Classification",
    parameters = params,
    grid = function(x, y, len = NULL, search = "grid") {
      grid <- expand.grid(alpha=seq(0, 0.99, len=len), delta=seq(0, max.delta, len=len))
      cbind(grid, data.frame(len=len))
    },
    sort = function(x) {
      x[order(x$alpha, x$delta),]
    },
    loop = function (grid) {
      grid <- grid[order(grid$alpha, grid$delta), , drop = F]
      loop <- grid[1, , drop = F]
      submodels <- list(grid[-1, , drop = F])
      list(loop = loop, submodels = submodels)
    },
    fit = function(x, y, wts, param, lev, last, classProbs, ...) {
      if (!is.factor(y))
        stop('Response must be a factor for SCRDA models')
      if (is.data.frame(x))
        x.names <- names(x)
      else
        x.names <- NULL
      if (!is.matrix(x))
        x <- as.matrix(x)
      x <- t(x)
      y <- as.integer(y)
      
      alpha <- seq(0, 0.99, len=param$len)
      delta <- seq(0, max.delta, len=param$len)
      m <- rda::rda(x, y, alpha=alpha, delta=delta, genelist = T)
      m$x.names <- x.names
      m$x.train <- x
      m$y.train <- y
      m$loop.params <- param
      m
    },
    predict = function(modelFit, newdata, submodels = NULL) {
      transform <- function(pred) {
        sapply(as.integer(pred), function(p) {
          if (p == 1) modelFit$obsLevel[1]
          else if (p == 2) modelFit$obsLevel[2]
          else NA
        })
      }
      get.predictions(modelFit, newdata, submodels, 'class', transform)
    },
    prob = function(modelFit, newdata, submodels = NULL) {
      transform <- function(pred) {
        dimnames(pred) <- unname(dimnames(pred))
        dimnames(pred)[[2]] <- modelFit$obsLevels
        pred[is.na(pred)] <- 0
        pred
      }
      get.predictions(modelFit, newdata, submodels, 'posterior', transform)
    },
    #varImp = ,
    predictors = function(x, ...) {
      browser()
      rownames(x$projection)
    },
    levels = function(x) {
      browser()
      x$obsLevels
    }
  )
  if (var.imp){
    m$varImp <- function(object, estimate = NULL, ...) {
      params <- list(...)
      if (is.null(params$alpha) || is.null(params$delta))
        stop('Variable importance for SCRDA is only possible when specifying alpha and delta arguments')
      
      var.imp <- predict.rda(object, x=object$x.train, y=object$y.train,
                             xnew=NULL, alpha=params$alpha, delta=params$delta, type='nonzero')
      var.imp <- data.frame(Overall=var.imp)
      if (!is.null(object$x.names))
        rownames(var.imp) <- object$x.names
      var.imp
    }
  }
  m
}


MOST_FREQUENT <- function(x, lvl) names(sort(table(x), decreasing=TRUE)[1])
CLASS_FREQUENCY <- function(x, lvl) sum(x == lvl[1])/length(x)
SPLIT_MEAN_PROB_ON_.5 <- function(x, lvl) ifelse(mean(x) < .5, lvl[2], lvl[1])
MEAN_PROB <- function(x, lvl) mean(x)

ENS_AVG_DEFAULT_CONVERTERS <- list(
  class.to.class=function(x, lvl) names(sort(table(x), decreasing=TRUE)[1]), 
  class.to.prob=function(x, lvl) sum(x == lvl[1])/length(x), 
  prob.to.class=function(x, lvl) ifelse(mean(x) < .5, lvl[2], lvl[1]),
  prob.to.prob=function(x, lvl) mean(x)
)

GetEnsembleAveragingModel <- function(
  class.to.class=MOST_FREQUENT, 
  class.to.prob=CLASS_FREQUENCY, 
  prob.to.class=SPLIT_MEAN_PROB_ON_.5,
  prob.to.prob=MEAN_PROB){
  list(
    label = "Averaging Model",
    library = NULL,
    loop = NULL,
    type = c("Classification"),
    parameters = data.frame(parameter = "parameter", class = "character", label = "parameter"),
    grid = function(x, y, len = NULL, search = "grid") data.frame(parameter="none"),
    fit = function(x, y, wts, param, lev, last, classProbs, ...) {
      list(lev=lev)
    },
    predict = function(modelFit, newdata, submodels = NULL) {
      all.factors <- all(sapply(newdata, is.factor)) || all(sapply(newdata, is.character))
      all.numeric <- all(sapply(newdata, is.numeric))
      if (!all.factors && !all.numeric) stop(paste(
        'Aggregating model only works if training data contains',
        'all factors OR all numeric values'
      ))
      
      if (all.factors) p <- apply(newdata, 1, class.to.class, modelFit$lev)
      else p <- apply(newdata, 1, prob.to.class, modelFit$lev)
      
      factor(p, levels=modelFit$lev)
    },
    prob = function(modelFit, newdata, submodels = NULL) {
      all.factors <- all(sapply(newdata, is.factor)) || all(sapply(newdata, is.character))
      all.numeric <- all(sapply(newdata, is.numeric))
      if (!all.factors && !all.numeric) stop(paste(
        'Aggregating model can only make probability predictions if',
        'new data contains all factors OR all numeric values'
      ))
      
      if (all.factors) p <- apply(newdata, 1, class.to.prob, modelFit$lev)
      else p <- apply(newdata, 1, prob.to.prob, modelFit$lev)
      
      p <- cbind(p, 1-p)
      dimnames(p)[[2]] <- modelFit$obsLevels
      p
    },
    varImp = NULL,
    predictors = function(x, ...) NULL,
    levels = function(x) if(any(names(x) == "obsLevels")) x$obsLevels else NULL,
    sort = NULL
  )
}

GetEnsembleQuantileModel <- function(){
  
  class.to.prob  <- function(x, lvl) sum(x == lvl[1])/length(x)
  class.to.class <- function(x, lvl) names(sort(table(x), decreasing=TRUE)[1])
  prob.to.class  <- function(x, lvl, p) ifelse(quantile(x, p) < .5, lvl[2], lvl[1])
  prob.to.prob   <- function(x, lvl, p) quantile(x, p)
  
  list(
    label = "Quantile Ensemble Model",
    library = NULL,
    loop = NULL,
    type = c("Classification"),
    parameters = data.frame(parameter = "quantile", class = "numeric", label = "Quantile"),
    grid = function(x, y, len = NULL, search = "grid") {
      if (search == "grid"){
        data.frame(quantile=seq(0, 1, len=len))
      } else {
        data.frame(quantile=runif(len))
      }
    },
    fit = function(x, y, wts, param, lev, last, classProbs, ...) {
      list(lev=lev, quantile=param$quantile)
    },
    predict = function(modelFit, newdata, submodels = NULL) {
      all.factors <- all(sapply(newdata, is.factor)) || all(sapply(newdata, is.character))
      all.numeric <- all(sapply(newdata, is.numeric))
      if (!all.factors && !all.numeric) stop(paste(
        'Aggregating model only works if training data contains',
        'all factors OR all numeric values'
      ))
      
      if (all.factors) p <- apply(newdata, 1, class.to.class, modelFit$lev)
      else p <- apply(newdata, 1, prob.to.class, modelFit$lev, modelFit$quantile)
      
      factor(p, levels=modelFit$lev)
    },
    prob = function(modelFit, newdata, submodels = NULL) {
      all.factors <- all(sapply(newdata, is.factor)) || all(sapply(newdata, is.character))
      all.numeric <- all(sapply(newdata, is.numeric))
      if (!all.factors && !all.numeric) stop(paste(
        'Aggregating model can only make probability predictions if',
        'new data contains all factors OR all numeric values'
      ))
      
      if (all.factors) p <- apply(newdata, 1, class.to.prob, modelFit$lev)
      else p <- apply(newdata, 1, prob.to.prob, modelFit$lev, modelFit$quantile)
      
      p <- cbind(p, 1-p)
      dimnames(p)[[2]] <- modelFit$obsLevels
      p
    },
    varImp = NULL,
    predictors = function(x, ...) NULL,
    levels = function(x) if(any(names(x) == "obsLevels")) x$obsLevels else NULL,
    sort = NULL
  )
}

GetCaretEnsembleModel <- function(caret.list.args, caret.stack.args){
  list(
    label = "Caret Ensemble Model",
    library = NULL,
    loop = NULL,
    type = c("Classification"),
    parameters = data.frame(parameter = "parameter", class = "character", label = "parameter"),
    grid = function(x, y, len = NULL, search = "grid") data.frame(parameter="none"),
    fit = function(x, y, wts, param, lev, last, classProbs, ...) {
      train.args <- c(list(x=x, y=y), caret.list.args)
      cl <- do.call('caretList', train.args)
      
      train.args <- caret.stack.args
      train.args$all.models <- cl
      cs <- do.call('caretStack', train.args)
      # attach things
      cs
    },
    predict = function(modelFit, newdata, submodels = NULL) {
      predict(modelFit, newdata=newdata, type='raw')
    },
    prob = function(modelFit, newdata, submodels = NULL) {
      p <- predict(modelFit, newdata=newdata, type='prob')
      p <- cbind(p, 1-p)
      dimnames(p)[[2]] <- modelFit$obsLevels
      p
    },
    varImp = NULL,
    predictors = function(x, ...) NULL,
    levels = function(x) if(any(names(x) == "obsLevels")) x$obsLevels else NULL,
    sort = NULL
  )
}

GetCaretEnsembleModelFunctional <- function(caret.list.args.fun, caret.stack.args.fun){
  list(
    label = "Caret Ensemble Model",
    library = NULL,
    loop = NULL,
    type = c("Classification"),
    parameters = data.frame(parameter = "parameter", class = "character", label = "parameter"),
    grid = function(x, y, len = NULL, search = "grid") data.frame(parameter="none"),
    fit = function(x, y, wts, param, lev, last, classProbs, ...) {
      caret.list.args <- caret.list.args.fun(x, y, wts, param, lev, last, classProbs, ...)
      cl <- do.call('caretList', caret.list.args)
      
      caret.stack.args <- caret.stack.args.fun(x, y, wts, param, lev, last, classProbs, ...)
      caret.stack.args$all.models <- cl
      cs <- do.call('caretStack', caret.stack.args)
      cs
    },
    predict = function(modelFit, newdata, submodels = NULL) {
      predict(modelFit, newdata=newdata, type='raw')
    },
    prob = function(modelFit, newdata, submodels = NULL) {
      p <- predict(modelFit, newdata=newdata, type='prob')
      p <- cbind(p, 1-p)
      dimnames(p)[[2]] <- modelFit$obsLevels
      p
    },
    varImp = NULL,
    predictors = function(x, ...) NULL,
    levels = function(x) if(any(names(x) == "obsLevels")) x$obsLevels else NULL,
    sort = NULL
  )
}
