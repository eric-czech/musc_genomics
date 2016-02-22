
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

GetSCRDAModel <- function(max.delta=3){
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
        pred
      }
      get.predictions(modelFit, newdata, submodels, 'posterior', transform)
    },
    varImp = function(object, estimate = NULL, ...) {
      params <- list(...)
      if (is.null(params$alpha) || is.null(params$delta))
        stop('Variable importance for SCRDA is only possible when specifying alpha and delta arguments')
      
      var.imp <- predict.rda(object, x=object$x.train, y=object$y.train,
                             xnew=NULL, alpha=params$alpha, delta=params$delta, type='nonzero')
      var.imp <- data.frame(Overall=var.imp)
      if (!is.null(object$x.names))
        rownames(var.imp) <- object$x.names
      var.imp
    },
    predictors = function(x, ...) {
      browser()
      rownames(x$projection)
    },
    levels = function(x) {
      browser()
      x$obsLevels
    }
  )
}


# getPLSWrappedModel <- function(model){
#   #' Wraps a given classifier in a PLS preprocessor
#   #' Taken from: http://stackoverflow.com/questions/21092895/how-to-custom-a-model-in-caret-to-perform-pls-classifer-two-step-classificaton
#   param.df <- cbind(model$parameters, data.frame(parameter='ncomp', class='numeric', label='#Components'))
#   m <- list(
#     label = "PLS-RF",
#     library = c("pls", "randomForest"),
#     type = "Classification",
#     parameters = data.frame(p),
#     grid = function(x, y, len = NULL) {
#       grid <- model$grid(x, y, len=len)
#       grid <- expand.grid(grid, ncomp = seq(1, min(ncol(x) - 1, len), by = 1))
#       grid
#     },
#     loop = NULL,
#     fit = function(x, y, wts, param, lev, last, classProbs, ...) { 
#       ## First fit the pls model, generate the training set scores,
#       ## then attach what is needed to the random forest object to 
#       ## be used later
#       pre <- plsda(x, y, ncomp = param$ncomp)
#       scores <- pls:::predict.mvr(pre, x, type = "scores")
#       mod <- model$fit(scores, y, wts, param, lev, last, classProbs, ...)
#       mod$projection <- pre$projection
#       mod
#     },
#     predict = function(modelFit, newdata, submodels = NULL) {       
#       scores <- as.matrix(newdata)  %*% modelFit$projection
#       model$predict(modelFit, scores, submodels=submodels, type='raw')
#     },
#     prob = function(modelFit, newdata, submodels = NULL) {
#       scores <- as.matrix(newdata)  %*% modelFit$projection
#       model$predict(modelFit, scores, submodels=submodels, type='prob')
#     },
#     varImp = function(object, estimate = NULL, ...) {
#       model$varImp(object, estimate=estimate, ...)
#     },
#     predictors = function(x, ...) rownames(x$projection),
#     levels = function(x) x$obsLevels,
#     sort = function(x) x[order(x[,1]),])
# }