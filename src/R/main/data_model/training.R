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
source('~/repos/portfolio/functional/ml/R/results.R')
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
preproc <- c('zv', 'center', 'scale')

# set.seed(SEED)
#d.samp <- d.prep.tr %>% sample_frac(.2)
#X <- d.samp[,sample(c(c.ge, c.cn, c.mu), replace=F, size = 100)]; y <- d.samp[,'response']

X <- d.prep.tr %>% select(-response, -tumor_id); 
y <- d.prep.tr[,'response']; y.bin <- DichotomizeOutcome(y)

trctrl <- function(index, ...) trainControl(
  summaryFunction=function(...) c(twoClassSummary(...), defaultSummary(...)),
  method = "cv", number = 5, savePredictions='final', index=index, returnData=F,
  ...)


# X <- X[1:25,c(GetFeatures(X, 'cn')[1:500], GetFeatures(X, 'ge')[1:500], GetFeatures(X, 'mu')[1:300])]
# y <- y[1:25]

trainer.i1 <- Trainer(cache.dir=file.path(CACHE_DIR, 'training_data'), 
                   cache.project=paste0(RESPONSE_TYPE, '.all'), seed=SEED)
trainer.i1$generateFoldIndex(y, CreateFoldIndex)
fold.data.gen <- GetFoldDataGenerator(preproc, F, n.core=8, 
                  sml.num.p=.0001, lrg.num.p=.01, sml.bin.p=.1, lrg.bin.p=.15)
trainer.i1$generateFoldData(X, y, fold.data.gen, GetDataSummarizer())

predict.reg.data.sml <- function(fit, d, i){ predict(fit, d$X.test.sml[,names(d$X.train.sml)]) }
predict.reg.data.lrg <- function(fit, d, i){ predict(fit, d$X.test.lrg[,names(d$X.train.lrg)]) }
predict.bin.data.sml <- function(fit, d, i){ predict(fit, d$X.test.sml[,names(d$X.train.sml)], type='prob')[,2] }
predict.bin.data.lrg <- function(fit, d, i){ predict(fit, d$X.test.lrg[,names(d$X.train.lrg)], type='prob')[,2] }
test.bin <- function(d) d$y.test.bin
test.reg <- function(d) d$y.test
predict.browser <- function(fit, d, i){ browser() }
setLevels <- function(x, lvls){ levels(x) <- lvls; x }

##### Classifiication Models #####

bin.models <- list()
tgt.metric <- 'Accuracy'

bin.model.svm.radial.sml <- list(
  name='bin.svm.radial.sml', predict=predict.bin.data.sml,
  train=function(d, idx, ...){
    registerDoMC(1)
    train(
      d$X.train.sml, d$y.train.bin,
      method='svmRadial', preProcess='zv', metric=tgt.metric,
      tuneLength=15, trControl = trctrl(idx, classProbs=T)
    )
  }
)
bin.models$svm.radial.sml <- trainer.i1$train(bin.model.svm.radial.sml, enable.cache=T)

# Run SVM on full set
# Run poly svm?


bin.model.pam <- list(
  name='bin.pam', predict=predict.bin.data.sml, test=test.bin,
  train=function(d, idx, ...){
    registerDoMC(1)
    train(
      d$X.train.sml, d$y.train.bin, 
      method='pam', preProcess='zv', metric=tgt.metric, 
      tuneGrid=data.frame(.threshold=c(0, 10^seq(-8, -1, length.out=10))),
      trControl = trctrl(idx, classProbs=T)
    )
  } 
)
bin.models$pam <- trainer.i1$train(bin.model.pam, enable.cache=F)


## Model components
bin.model.pls <- list(
  name='bin.pls', predict=predict.bin.data.sml, test=test.bin,
  train=function(d, idx, ...){
    registerDoMC(1)
    train(
      d$X.train.sml, d$y.train.bin, 
      method=GetPLSModel(), preProcess='zv', metric=tgt.metric, 
      trControl = trctrl(idx, classProbs=T),
      tuneGrid=data.frame(ncomp=c(1,2,3,6,9,12,15,18))
    )
  }
)
bin.models$pls <- trainer.i1$train(bin.model.pls, enable.cache=F)

bin.model.knn <- list(
  name='bin.knn', predict=predict.bin.data.sml, test=test.bin,
  train=function(d, idx, ...){
    registerDoMC(1)
    train(
      d$X.train.sml, d$y.train.bin, 
      method='knn', preProcess=c('zv', 'pca'), metric=tgt.metric, 
      tuneLength=15, trControl = trctrl(idx, classProbs=T, preProcOptions = list(pcaComp = 25))
    )
  }
)
bin.models$knn <- trainer.i1$train(bin.model.knn, enable.cache=F)


GetElasticNetModel <- function(name, alpha.start, alpha.grid, n.core){
  list(
    name=name, predict=predict.bin.data.sml, test=test.bin,
    train=function(d, idx, ...){
      registerDoMC(n.core)
      glmnet.lambda <- GetGlmnetLambda(d$X.train.sml, d$y.train.bin, alpha=alpha.start, family='binomial')
      train(
        d$X.train.sml, d$y.train.bin, 
        method='glmnet', preProcess='zv', family='binomial', metric=tgt.metric, 
        tuneGrid = expand.grid(.alpha = alpha.grid, .lambda=glmnet.lambda),
        trControl = trctrl(idx, classProbs=T)
      )
    }
  )
}

bin.model.lasso <- GetElasticNetModel('bin.lasso', 1, 1, 1)
bin.models$lasso <- trainer.i1$train(bin.model.lasso, enable.cache=F)

bin.model.ridge <- GetElasticNetModel('bin.ridge', 0, 0, 1)
bin.models$ridge <- trainer.i1$train(bin.model.ridge, enable.cache=F)

bin.model.enet <- GetElasticNetModel('bin.enet', .5, c(.001, seq(.1, .9, length.out=13), .999), 1)
bin.models$enet <- trainer.i1$train(bin.model.enet, enable.cache=F)

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
# bin.models$nb <- trainer.i1$train(bin.model.nb)

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
# bin.models$rf <- trainer.i1$train(bin.model.rf)

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
# bin.models$gbm <- trainer.i1$train(bin.model.gbm)


##### Classification Ensembles #####

GetBinEnsemble <- function(models, name, use.index){
  list(
    name=name, test=test.bin,
    train=function(d, idx, i, ...){
      m <- lapply(models, function(m) m(i))
      class(m) <- "caretList"
      caretEnsemble(m)
    }, predict=function(fit, d, i){ 
      predict(fit, newdata=d$X.test.sml[,names(d$X.test.sml)]) 
      #predict(fit, newdata=d$X.test.sml[,names(d$X.train.sml)]) 
    }
  )
}

ens.models <- list(
  pam=function(i) bin.models$pam[[i]]$fit,
  knn=function(i) bin.models$knn[[i]]$fit,
  enet=function(i) bin.models$enet[[i]]$fit
  #pls=function(i) bin.models$pls[[i]]$fit # This fails because the method type is 'custom'
)
bin.model.ens1 <- GetBinEnsemble(ens.models, 'bin.ens1')
bin.models$bin.model.ens1 <- trainer.i1$train(bin.model.ens1, enable.cache=F)

##### Classification Hold Out #####

X.ho <- d.prep.ho %>% select(-response, -tumor_id); 
y.ho <- d.prep.ho[,'response']; y.ho.bin <- DichotomizeOutcome(y.ho)

GetHoldOutPrediction <- function(model, X, y, X.ho, y.ho){
  loginfo('Creating hold out predictions for model "%s"', model$name)
  res <- trainer.i1$holdout(model, X, y, X.ho, y.ho, fold.data.gen)
  data.frame(model=model$name, y.pred=res$y.pred, y.test=res$y.test)
}


models <- list(bin.model.pam, bin.model.knn, bin.model.enet)
ho.fit <- trainer.i1$holdout(models, X, y, X.ho, y.ho, fold.data.gen)

ens.models.ho <- lapply(ho.fit, function(m) {function(i) m$fit}) %>% setNames(sapply(ho.fit, function(m) m$model))
bin.model.ens1.ho <- GetBinEnsemble(ens.models.ho, 'bin.ens1.ho')
ens.ho.fit <- trainer.i1$holdout(list(bin.model.ens1.ho), X, y, X.ho, y.ho, fold.data.gen)

ho.preds <- foreach(p=c(ho.fit, ens.ho.fit)) %do% {
  data.frame(model=p$model, y.pred=p$y.pred, y.test=p$y.test)
}

# ho.data <- bs.data.gen(X.ho, y.ho.bin, X.ho, y.ho.bin)

##### Classification Results ##### 

## CV Results

res <- SummarizeResults(bin.models, fold.summary=GetResultSummary, model.summary=GetResultSummary)

model.perf <- res$fold.summary %>% group_by(model, fold) %>%
  summarise_each(funs(head(., 1)), one_of(c('auc', 'acc.max', 'acc.cut'))) %>%
  ungroup %>% group_by(model) %>% 
  summarise_each(funs(mean, sd), one_of(c('auc', 'acc.max', 'acc.cut')))

# ROC curves per-fold
res$fold.summary %>% 
  select(model, fold, x, y, t) %>%
  inner_join(model.perf, by='model') %>%
  filter(str_detect(model, '.')) %>% 
  mutate(model.label=paste0(model, ' (', round(auc_mean, 3), ')')) %>% 
  ggplot(aes(x=x, y=y, color=factor(fold))) + 
  geom_line() + geom_abline(alpha=.5) + theme_bw() + facet_wrap(~model.label)

# ROC curves across folds
res$model.summary %>%
  select(model, x, y, t) %>%
  inner_join(model.perf, by='model') %>%
  filter(str_detect(model, '.')) %>% 
  mutate(model.label=paste0(model, ' (', round(auc_mean, 3), ')')) %>% 
  ggplot(aes(x=x, y=y, color=factor(model.label))) + 
  geom_line() + geom_abline(alpha=.5) + theme_bw() 

# AUC ranges by model
model.perf %>% arrange(auc_mean) %>%
  mutate(model=factor(model, levels=model)) %>%
  ggplot(aes(x=model, y=auc_mean, ymin=auc_mean - auc_sd, ymax=auc_mean + auc_sd, color=model)) +
  geom_pointrange() + coord_flip() + theme_bw() 


## hold-out results

ho.model.perf <- foreach(preds=ho.preds, .combine=rbind) %do% {
  GetResultSummary(preds) %>% mutate(model=preds$model[1])
}

ho.model.perf %>% group_by(model) %>% do({head(., 1)}) %>%
  select(auc, acc.max, model) %>% melt(id.vars = 'model') %>%
  ggplot(aes(x=model, y=value, fill=model)) + geom_bar(position='dodge', stat='identity') +
  facet_wrap(~variable)

ho.model.perf %>%
  select(model, x, y, t) %>%
  inner_join(ho.model.perf %>% group_by(model) %>% summarise(auc=auc[1]), by='model') %>% 
  mutate(model.label=paste0(model, ' (', round(auc, 3), ')')) %>% 
  ggplot(aes(x=x, y=y, color=factor(model.label))) + 
  geom_line() + geom_abline(alpha=.5) + theme_bw() 


##### Regression Models #####

reg.trctrl <- function(index, ...) trainControl(
  summaryFunction=function(...) defaultSummary(...),
  method = "cv", number = 5, savePredictions='final', index=index, returnData=F,
  ...)

reg.models <- list()

reg.model.svm.radial <- list(
  name='reg.svm.radial', predict=predict.reg.data.sml, test=test.reg,
  train=function(d, idx, ...){
    registerDoMC(3)
    train(
      d$X.train.sml, d$y.train, 
      method='pam', preProcess='zv', metric=tgt.metric, 
      tuneLength=15, search='grid',
      trControl = reg.trctrl(idx)
    )
  }
)
reg.models$svm.radial <- trainer.i1$train(reg.model.svm.radial, enable.cache=F)

reg.model.pls <- list(
  name='reg.pls', predict=predict.reg.data.sml, test=test.reg,
  train=function(d, idx, ...){
    registerDoMC(1)
    train(
      d$X.train.sml, d$y.train, method='pls', 
      tuneGrid=data.frame(ncomp=c(1,2,3,6,9,12,15,18)),
      trControl = reg.trctrl(idx)
    )
  }
)
reg.models$pls <- trainer.i1$train(reg.model.pls, enable.cache=F)

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

# model.gbm <- list(
#   train=function(d){
#     registerDoMC(3)
#     train(
#       d$X.train.sml, d$y.train, method='gbm',  
#       trControl = trctrl(verboseIter=T), tuneLength=10
#       
#     )
#   }, predict=predict.test.data
# )
# models$gbm <- trainer.rs$train(model.gbm)

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




