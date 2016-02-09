#'-----------------------------------------------------------------------------
#' ML Model Training Script
#'
#' This module contains code for training a variety of regression models on 
#' prepared genomics datasets (post feature-selection).
#' 
#' @author eczech
#'-----------------------------------------------------------------------------
source('utils.R')
source('data_model/training_lib.R')
source('~/repos/portfolio/functional/ml/R/trainer.R')
lib('MASS')
lib('caret')
lib('doMC')
lib('iterators')
lib('caretEnsemble')
lib('ROCR')


#registerDoMC(8)
SEED <- 1024

RESPONSE_TYPE <- 'cosmic' # This will include ctd2 as well at some point
RESPONSE_SELECTOR <- function(d){ 
  d %>% filter(!is.na(ic_50)) %>% rename(response=ic_50) %>% select(-auc)
}
select <- dplyr::select

d.prep <- GetTrainingData(TRAIN_CACHE, RESPONSE_TYPE, RESPONSE_SELECTOR, min.mutations=3)

set.seed(SEED)
d.prep.tr <- d.prep %>% sample_frac(.8, replace = F)
d.prep.ho <- d.prep %>% filter(!tumor_id %in% d.prep.tr$tumor_id)
#preproc <- c('zv', 'center', 'scale')
preproc <- c('zv', 'center', 'scale')

# set.seed(SEED)
#d.samp <- d.prep.tr %>% sample_frac(.2)
#X <- d.samp[,sample(c(c.ge, c.cn, c.mu), replace=F, size = 100)]; y <- d.samp[,'response']

X <- d.prep.tr %>% select(-response, -tumor_id); 
y <- d.prep.tr[,'response']; y.bin <- factor((sign(y) + 1) * .5)

trctrl <- function(index, ...) trainControl(
  summaryFunction=function(...) c(twoClassSummary(...), defaultSummary(...)),
  method = "cv", number = 5, savePredictions='none', index=index, returnData=F,
  ...)


# Changes for next run of preprocessing
# 1. Compute both y and y.bin
# 2. Compute both large and small feature sets
# 3. Should be two trainers resulting, one for regression and one for classification
# 4. Do this to y.bin --> setLevels(d$y.train, c('neg', 'pos'))

X <- X[1:25,1:1000]
y <- y[1:25]
trainer.v1 <- Trainer(cache.dir=file.path(CACHE_DIR, 'training_data'), 
                   cache.project=paste0(RESPONSE_TYPE, '.all'), seed=SEED)
trainer.v1$generateFoldIndex(y, CreateFoldIndex)
fold.data.gen <- GetFoldDataGenerator(preproc, F, n.core=8, 
                  sml.num.p=.01, lrg.num.p=.05, sml.bin.p=.15, lrg.bin.p=.15)
trainer.v1$generateFoldData(X, y, fold.data.gen, GetDataSummarizer())


trainer.rs <- Trainer(cache.dir=file.path(CACHE_DIR, 'training_data'), 
                   cache.project=paste0(RESPONSE_TYPE, '.reg.sm'), seed=SEED)
trainer.rs$generateFoldIndex(y, CreateFoldIndex)
rs.data.gen <- GetFoldDataGenerator(preproc, 'numeric', T, n.core=8, numeric.score.p=.0001, binary.score.p=.15)
trainer.rs$generateFoldData(X, y, rs.data.gen, GetDataSummarizer())

trainer.bs <- Trainer(cache.dir=file.path(CACHE_DIR, 'training_data'), 
                      cache.project=paste0(RESPONSE_TYPE, '.bin.sm'), seed=SEED)
trainer.bs$generateFoldIndex(y.bin, CreateFoldIndex)
bs.data.gen <- GetFoldDataGenerator(preproc, 'binary', T, n.core=1, numeric.score.p=.01, binary.score.p=.15)
trainer.bs$generateFoldData(X, y.bin, bs.data.gen, GetDataSummarizer())


predict.reg.data <- function(fit, d, i){ predict(fit, d$X.test.sml[,names(d$X.train.sml)]) }
predict.bin.data <- function(fit, d, i){ predict(fit, d$X.test.sml[,names(d$X.train.sml)], type='prob')[,2] }
predict.browser <- function(fit, d, i){ browser() }
setLevels <- function(x, lvls){ levels(x) <- lvls; x }

##### Classifiication Models #####

bin.models <- list()
tgt.metric <- 'ROC'

bin.model.svm.radial.sml <- list(
  name='bin.svm.radial.sml', predict=predict.bin.data,
  train=function(d, idx, ...){
    registerDoMC(1)
    train(
      d$X.train.sml, setLevels(d$y.train, c('neg', 'pos')), method='svmRadial', 
      metric=tgt.metric, tuneLength=15, trControl = trctrl(idx, classProbs=T)
    )
  }
)
bin.models$svm.radial.sml <- trainer.bs$train(bin.model.svm.radial.sml)

# Run SVM on full set
# Run poly svm?

bin.model.pam <- list(
  name='bin.pam', predict=predict.bin.data,
  train=function(d, idx, ...){
    registerDoMC(3)
    train(
      d$X.train.sml, setLevels(d$y.train, c('neg', 'pos')), method='pam', metric=tgt.metric, 
      tuneGrid=data.frame(.threshold=c(0, 10^seq(-8, -1, length.out=10))),
      trControl = trctrl(idx, classProbs=T)
    )
  } 
)
bin.models$pam <- trainer.bs$train(bin.model.pam, enable.cache=T)


## Model components
bin.model.pls <- list(
  name='bin.pls', predict=predict.bin.data,
  train=function(d, idx, ...){
    registerDoMC(1)
    train(
      d$X.train.sml, setLevels(d$y.train, c('neg', 'pos')), 
      method=GetPLSModel(), metric=tgt.metric, 
      trControl = trctrl(idx, classProbs=T),
      tuneGrid=data.frame(ncomp=c(1,2,3,4))
    )
  }
)
bin.models$pls <- trainer.bs$train(bin.model.pls, enable.cache=F)

bin.model.knn <- list(
  name='bin.knn', predict=predict.bin.data,
  train=function(d, idx, ...){
    registerDoMC(3)
    train(
      d$X.train.sml, setLevels(d$y.train, c('neg', 'pos')), method='knn', 
      metric=tgt.metric, tuneLength=20, trControl = trctrl(idx, classProbs=T)
    )
  }
)
bin.models$knn <- trainer.bs$train(bin.model.knn)

bin.model.glmnet <- list(
  name='bin.glmnet', predict=predict.bin.data,
  train=function(d, idx, ...){
    registerDoMC(1)
    y.train <- setLevels(d$y.train, c('neg', 'pos'))
    glmnet.lambda <- GetGlmnetLambda(d$X.train.sml, y.train, family='binomial')
    train(
      d$X.train.sml, y.train, method='glmnet', family='binomial', metric=tgt.metric, 
      tuneGrid = expand.grid(.alpha = seq(0, .001, length.out = 10), .lambda=glmnet.lambda),
      trControl = trctrl(idx, classProbs=T)
    )
  }
)
bin.models$glmnet <- trainer.bs$train(bin.model.glmnet)

# bin.model.nb <- list(
#   name='bin.nb', predict=predict.bin.data,
#   train=function(d, idx, ...){
#     registerDoMC(1)
#     train(
#       d$X.train.sml, setLevels(d$y.train, c('neg', 'pos')), method='nb', 
#       metric=tgt.metric, tuneLength=2, trControl = trctrl(idx, classProbs=T)
#     )
#   }
# )
# bin.models$nb <- trainer.bs$train(bin.model.nb)

# bin.model.rf <- list(
#   name='bin.rf', predict=predict.bin.data,
#   train=function(d, idx, ...){
#     registerDoMC(1)
#     train(
#       d$X.train.sml, setLevels(d$y.train, c('neg', 'pos')), method='rf', metric=tgt.metric,
#       #tuneGrid=data.frame(ntree=c(100, 250), mtry=c(45, 45)),
#       #tuneLength=8, 
#       trControl = trctrl(idx, classProbs=T),
#       ntree=100, mtry=45
#     )
#   }
# )
# bin.models$rf <- trainer.bs$train(bin.model.rf)

# bin.model.gbm <- list(
#   name='bin.gbm', predict=predict.bin.data,
#   train=function(d, idx, ...){
#     registerDoMC(3)
#     train(
#       d$X.train.sml, setLevels(d$y.train, c('neg', 'pos')), method='gbm', metric=tgt.metric,
#       tuneLength=5, trControl = trctrl(idx, classProbs=T), verbose=F
#     )
#   }
# )
# bin.models$gbm <- trainer.bs$train(bin.model.gbm)


##### Classification Ensembles #####

bin.model.ens1 <- list(
  train=function(d, idx, ...){
    m <- list(
      bin.models$glmnet[[i]]$fit,
      bin.models$pls[[i]]$fit,
      bin.models$pam[[i]]$fit,
      bin.models$knn[[i]]$fit
    )
    class(m) <- "caretList"
    caretEnsemble(m)
  }, predict=function(fit, d, i){ 
    predict(fit, newdata=d$X.test.sml[,names(d$X.train.sml)]) 
  }
)
bin.models$bin.model.ens1 <- trainer.bs$train(bin.model.ens1)

##### Classification Hold Out #####

registerDoMC(8)
#dt <- trainer.bs$fit(bin.model.glmnet, X, y.bin, bs.data.gen)
dt <- trainer.bs$fit(bin.model.svm.radial.sml, X, y.bin, bs.data.gen)

X.ho <- d.prep.ho %>% select(-response, -tumor_id); 
y.ho <- d.prep.ho[,'response']; y.ho.bin <- factor((sign(y.ho) + 1) * .5)

X.ho.prep <- predict(dt$data$preproc$pp.sml, X.ho[,names(dt$data$X.train.sml)])
y.ho.pred <- predict(dt$fit, X.ho.prep, type='prob')[,2]

roc <- performance(prediction(y.ho.pred, as.integer(y.ho.bin)-1), 'tpr', 'fpr')
auc <- performance(prediction(y.ho.pred, as.integer(y.ho.bin)-1), 'auc')
plot(roc)
abline(a = 0, b=1)

# ho.data <- bs.data.gen(X.ho, y.ho.bin, X.ho, y.ho.bin)

##### Classification Results ##### 

cv.preds <- foreach(m=names(bin.models), .combine=rbind) %do% {
  foreach(fold=bin.models[[m]], .combine=rbind) %do% {
    y.test <- as.integer(fold$y.test) - 1
    data.frame(fold=fold$fold, y.pred=fold$y.pred, y.test=y.test, model=m)
  }
}

roc <- cv.preds %>% group_by(model, fold) %>% do({
  pred <- prediction(.$y.pred, .$y.test)
  roc <- performance(pred, 'tpr', 'fpr') 
  auc <- performance(pred, 'auc')
  data.frame(
    x=roc@x.values[[1]], y=roc@y.values[[1]], 
    t=roc@alpha.values[[1]], auc=auc@y.values[[1]]
  )
}) %>% ungroup
roc.auc <- roc %>% group_by(model, fold) %>% summarise(auc=auc[1]) %>% ungroup %>%
  group_by(model) %>% summarise(auc.mean=mean(auc), auc.se=sd(auc))
roc %>% inner_join(roc.auc, by='model') %>% 
  filter(str_detect(model, '.')) %>% 
  mutate(model.label=paste0(model, ' (', round(auc.mean, 3), ')')) %>% 
  ggplot(aes(x=x, y=y, color=factor(fold))) + 
  geom_line() + geom_abline(alpha=.5) + theme_bw() + facet_wrap(~model.label)



##### Regression Models #####

models <- list()

model.gbm <- list(
  train=function(d){
    registerDoMC(3)
    train(
      d$X.train.sml, d$y.train, method='gbm',  
      trControl = trctrl(verboseIter=T), tuneLength=10
      
    )
  }, predict=predict.test.data
)
models$gbm <- trainer.rs$train(model.gbm)

model.svm.radial <- list(
  train=function(d){
    registerDoMC(3)
    train(
      d$X.train.sml, d$y.train, method='svmRadial', 
      tuneLength=25, search='grid',
      trControl = trctrl(verboseIter=T)
    )
  }, predict=predict.test.data
)
models$svm.radial <- trainer.rs$train(model.svm.radial)

model.pls <- list(
  train=function(d){
    registerDoMC(3)
    train(
      d$X.train.sml, d$y.train, method='pls', 
      tuneGrid=data.frame(ncomp=c(1,2,3)),
      trControl = trctrl(verboseIter=T)
    )
  }, predict=predict.test.data
)
models$pls <- trainer.rs$train(model.pls)

model.glmnet <- list(
  train=function(d){
    registerDoMC(5)
    glmnet.lambda <- GetGlmnetLambda(d$X.train.sml, d$y.train)
    train(
      d$X.train.sml, d$y.train, method='glmnet',
      tuneGrid = expand.grid(.alpha = seq(0, .001, length.out = 2), .lambda=glmnet.lambda),
      trControl = trctrl(verboseIter=T)
    )
  }, predict=predict.test.data
)
models$glmnet <- trainer.rs$train(model.glmnet)


model.ens1 <- list(
  train=function(d){
    registerDoMC(3)
    glmnet.lambda <- GetGlmnetLambda(d$X.train.sml, d$y.train)
    caretList(
      d$X.train.sml, d$y.train, metric='RMSE',
      trControl=trainControl(method = "cv", number = 5, savePredictions=T),
      tuneList=list(
        glmnet=caretModelSpec(method='glmnet', tuneGrid=expand.grid(.alpha = seq(0, .001, length.out = 25), .lambda=glmnet.lambda)),
        pls=caretModelSpec(method='pls', tuneGrid=data.frame(ncomp=c(1,2,3))),
        svm.radial=caretModelSpec(method='svmRadial', tuneLength=25)
      )
    )
  }, predict=function(fit, d){ 
    predict(caretEnsemble(fit), newdata=d$X.test.sml[,names(d$X.train.sml)]) 
  }
)
models$ens.1 <- trainer.rs$train(model.ens1)

# registerDoMc(3)
# m1 <- caretList(
#   X[1:50,1:50], y[1:50], 
#   trControl=trainControl(method = "cv", number = 5, savePredictions=T),
#   methodList=c('pls', 'glmnet'), 
#   tuneList=list(
#     pls=caretModelSpec(method='pls', tuneGrid=data.frame(ncomp=c(1,2,3))),
#     glmnet=caretModelSpec(method='glmnet', tuneLength=3)
#   )
# )


foreach(fold=models$glmnet) %do% fold$fit$bestTune
dt <- foreach(fold=fit.glmnet, .combine=rbind)%do%{
  fold$fit$results
}
svm.grid.fit <- trainer$train(svm.grid)

cv.scores <- foreach(m=names(models), .combine=rbind) %do% {
  foreach(fold=models[[m]], .combine=rbind) %do% {
    #browser()
    y.pred <- fold$y.pred
    y.test <- fold$y.test
    r2   <- GetRegressorScore(y.test, y.pred, 'rsquared')
    rmse <- GetRegressorScore(y.test, y.pred, 'rmse')
    data.frame(model=m, rsquared=r2, rmse=rmse)
  }
}

cv.scores %>% ggplot(aes(x=model, y=rmse)) + geom_point() + geom_boxplot()
cv.scores %>% ggplot(aes(x=model, y=rsquared)) + geom_point() + geom_boxplot()




