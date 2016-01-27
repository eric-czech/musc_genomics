library(caret)
library(foreach)

PreProcessor <- setClass(
  "PreProcessor",
  slots = c(method = "character", pps = "list"),
  prototype=list(method = c()),
  validity=function(object){
    return(TRUE)
  }
)

setGeneric(name="Fit",def=function(object, data){standardGeneric("Fit")})
setMethod(
  f="Fit",
  signature="PreProcessor",
  definition=function(object, data){
    pps <- foreach(f=names(data)) %dopar%{
      pp <- preProcess(subset(data, select=f), method=object@method)
      list(col=f, pp=pp)
    }
    names(pps) <- sapply(pps, function(x) x$col)
    object@pps <- pps
    return(object)
  }
)

setGeneric(name="Transform",def=function(object, data){standardGeneric("Transform")})
setMethod(
  f="Transform",
  signature="PreProcessor",
  definition=function(object, data){
    browser()
    trans.data <- foreach(f=names(data), .combine=cbind) %dopar%{
      keys <- names(object@pps)
      if (!f %in% keys) NULL
      else predict(object@pps[[f]]$pp, data[,f])
    }
    return(trans.data)
  }
)

scaler <- PreProcessor(method=c('center', 'scale'))
tmp <- Fit(scaler, d)
Transform(tmp, d)

d <- data.frame(x=1:10, y=2:11)
