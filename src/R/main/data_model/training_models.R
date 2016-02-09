GetPLSModel <- function(){
  plsModel <- getModelInfo(model = "pls", regex = FALSE)[[1]]
  plsModel$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {   
    out <- if(is.factor(y)){
      plsda(x, y, method = "oscorespls", ncomp = param$ncomp, model=F, ...)
    } else {
      dat <- if(is.data.frame(x)) x else as.data.frame(x)
      dat$.outcome <- y
      plsr(.outcome ~ ., data = dat, method = "oscorespls", ncomp = param$ncomp, ...)
    }
    out
  }
  plsModel
}

