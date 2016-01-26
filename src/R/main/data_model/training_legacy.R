#'-----------------------------------------------------------------------------
#' Defunct Training Routines and Code Samples
#'
#' This module contains examples of model training with caret that involve 
#' feature selection as well as hyperparameter tuning.  Unfortunately, the
#' caret interface was not flexible to support these things in conjunction with
#' ensembles so an explicit resampling loop was written instead that then calls
#' simpler caret methods like 'train', rather than 'sbf' which attempts to do
#' more.
#' 
#' @author eczech
#'-----------------------------------------------------------------------------

TrainModels <- function(d, response){
  # PLS 
  set.seed(SEED); registerDoMC(8)
  m.pls.1 <- sbf(
    X, y, method='pls', preProc = preproc, 
    tuneGrid=data.frame(ncomp=c(1,2,3)),
    trControl = trctrl(), sbfControl = fsctrl(functions = fs.l)
  )
  
  
  # GLMNET
  set.seed(SEED); registerDoMC(5)
  lambda <- GetGlmnetLambda(X.preproc)
  m.glmnet.1 <- sbf(
    X, y, method='glmnet', preProc = c(preproc), 
    tuneGrid=expand.grid(.alpha = seq(0, .12, by=.01), .lambda = lambda),
    #tuneGrid=expand.grid(.alpha = 10^seq(-1,-5,length=5), .lambda=10^seq(1,-5,length=100)),
    #tuneLength=25,
    trControl = trctrl(), 
    sbfControl = fsctrl(functions = fs.s),
    standardize=F
  )
  save(m.glment.1, file=file.path(MODEL_DIR, 'm.glmnet.1.Rdata'))
  
  
  # GLMNET w/ PCA
  set.seed(SEED)
  m.glmnet.2 <- train(
    X, y, method='glmnet', preProc = c(preproc, 'pca'),
    tuneGrid=expand.grid(.alpha = seq(0, .12, by=.01), .lambda = lambda),
    trControl = trctrl(preProcOptions = list(thresh = 0.3)), 
    type.gaussian="naive", standardize=F
  )
  save(m.glment.2, file=file.path(MODEL_DIR, 'm.glmnet.1.Rdata'))
  
  # GBM
  set.seed(SEED); registerDoMC(1)
  gbmGrid <-  expand.grid(interaction.depth = c(1, 2, 3),
                          n.trees = c(10, 50, 100),
                          shrinkage = 0.1,
                          n.minobsinnode = c(5, 10))
  m.gbm.1 <- sbf(
    X, y, method='gbm', 
    tuneGrid=gbmGrid,
    trControl = trctrl(), 
    sbfControl = fsctrl(functions = fs.s)
  )
  
  # SVM
  library(kernlab)
  set.seed(SEED); registerDoMC(1)
  #X <- X[1:10,]; y <- y[1:10]
  sigDist <- sigest(as.matrix(X.preproc), frac = 1)
  svmTuneGrid <- data.frame(.sigma = sigDist[1], .C = 2^(-2:6))
  m.svm.1 <- sbf(
    X, y, method='svmRadial', preProc = c("center", "scale"), 
    tuneGrid = svmTuneGrid,
    trControl = trctrl(verboseIter=T), 
    sbfControl = fsctrl(functions = fs.s, verbose=T)
  )
  
  set.seed(SEED); registerDoMC(1)
  m.svm.2 <- train(
    X, y, method='svmRadial', preProc = c("center", "scale"), 
    tuneGrid = data.frame(.sigma = as.numeric(sigDist[1]), .C = 2^(-2:2)),
    trControl = trctrl(verboseIter=T),
    maxiter=10
  )
  
  list(pls1=m.pls.1, pls1=m.pls.2, glmnet1=m.glmnet.1, glmnet2=m.glmnet.2)
}