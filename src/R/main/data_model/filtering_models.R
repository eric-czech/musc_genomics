

GetFitEnsembleTrain <- function(model.results, fit.prep.fun=NULL, ...){
  list(
    train=function(d, idx, i){
      models <- lapply(model.results, function(m) m[[i]]$fit) %>% setNames(names(model.results))
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
