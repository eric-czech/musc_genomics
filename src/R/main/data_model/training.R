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
source('data_model/training_models.R')
source('~/repos/portfolio/functional/ml/R/trainer.R')
source('~/repos/portfolio/functional/ml/R/results.R')
lib('MASS')
lib('caret')
lib('doMC')
lib('iterators')
lib('caretEnsemble')
lib('ROCR')
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


# X <- X[1:25,c(GetFeatures(X, 'cn')[1:500], GetFeatures(X, 'ge')[1:500], GetFeatures(X, 'mu')[1:300])]
# y <- y[1:25]

trainer.i1 <- Trainer(cache.dir=file.path(CACHE_DIR, 'training_data'), 
                   cache.project=paste0(RESPONSE_TYPE, '.all'), seed=SEED)
trainer.i1$generateFoldIndex(y, CreateFoldIndex)
fold.data.gen <- GetFoldDataGenerator(preproc, F, n.core=8, 
                    sml.num.p=.0001, lrg.num.p=.01, sml.bin.p=.1, lrg.bin.p=.15)
trainer.i1$generateFoldData(X, y, fold.data.gen, GetDataSummarizer())


##### Classifiication Models #####

bin.models <- list()

# ShowBestTune(bin.models$gbm)

bin.models$svm.radial.sml <- trainer.i1$train(bin.model.svm.radial.sml, enable.cache=T)
bin.models$knn <- trainer.i1$train(bin.model.knn, enable.cache=T)
bin.models$knn.pca <- trainer.i1$train(bin.model.knn.pca, enable.cache=T)
bin.models$pam <- trainer.i1$train(bin.model.pam, enable.cache=T)
bin.models$pls <- trainer.i1$train(bin.model.pls, enable.cache=T)
bin.models$rf <- trainer.i1$train(bin.model.rf, enable.cache=T)
bin.models$lasso <- trainer.i1$train(bin.model.lasso, enable.cache=T)
bin.models$ridge <- trainer.i1$train(bin.model.ridge, enable.cache=T)
bin.models$enet <- trainer.i1$train(bin.model.enet, enable.cache=T)


bin.models$gbm <- trainer.i1$train(bin.model.gbm, enable.cache=F)
bin.models$rda <- trainer.i1$train(bin.model.rda, enable.cache=F)
bin.models$et <- trainer.i1$train(bin.model.et, enable.cache=F)



# bin.models$nb <- trainer.i1$train(bin.model.nb)
# bin.models$rf <- trainer.i1$train(bin.model.rf)
# bin.models$gbm <- trainer.i1$train(bin.model.gbm)


##### Classification Ensembles #####

ens.models <- list(
  pam=function(i) bin.models$pam[[i]]$fit,
  knn=function(i) bin.models$knn[[i]]$fit,
  enet=function(i) bin.models$enet[[i]]$fit
  #pls=function(i) bin.models$pls[[i]]$fit # This fails because the method type is 'custom'
)
bin.model.ens1 <- GetEnsembleModel(ens.models, 'bin.ens1', test.bin, predict.bin.ens.sml)
bin.models$bin.model.ens1 <- trainer.i1$train(bin.model.ens1, enable.cache=F)

##### Classification Hold Out #####

X.ho <- d.prep.ho %>% select(-response, -tumor_id); 
y.ho <- d.prep.ho[,'response']; y.ho.bin <- DichotomizeOutcome(y.ho)

GetHoldOutPrediction <- function(model, X, y, X.ho, y.ho){
  loginfo('Creating hold out predictions for model "%s"', model$name)
  res <- trainer.i1$holdout(model, X, y, X.ho, y.ho, fold.data.gen)
  data.frame(model=model$name, y.pred=res$y.pred, y.test=res$y.test)
}

models <- list(bin.model.rf, bin.model.svm.radial.sml, bin.model.pam, 
               bin.model.knn, bin.model.enet, bin.model.lasso, bin.model.ridge)
ho.fit <- trainer.i1$holdout(models, X, y, X.ho, y.ho, fold.data.gen)

ens.models.ho <- lapply(ho.fit, function(m) {function(i) m$fit}) %>% setNames(sapply(ho.fit, function(m) m$model))
bin.model.ens1.ho <- GetBinEnsemble(ens.models.ho, 'bin.ens1.ho')
ens.ho.fit <- trainer.i1$holdout(list(bin.model.ens1.ho), X, y, X.ho, y.ho, fold.data.gen)

ho.preds <- foreach(p=c(ho.fit, ens.ho.fit)) %do% {
  data.frame(model=p$model, y.pred=p$y.pred, y.test=p$y.test)
}
ho.preds <- foreach(p=ho.fit) %do% {
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
library(plotly)
p <- res$model.summary %>%
  select(model, x, y, t) %>%
  inner_join(model.perf, by='model') %>%
  filter(str_detect(model, '.')) %>% 
  mutate(model.label=paste0(model, ' (', round(auc_mean, 3), ')')) %>% 
  ggplot(aes(x=x, y=y, color=factor(model.label))) + 
  geom_line() + geom_abline(alpha=.5) + theme_bw() 
p
#ggplotly(p)

# AUC ranges by model
model.perf %>% arrange(auc_mean) %>%
  mutate(model=factor(model, levels=model)) %>%
  ggplot(aes(x=model, y=auc_mean, ymin=auc_mean - auc_sd, ymax=auc_mean + auc_sd, color=model)) +
  geom_pointrange() + coord_flip() + theme_bw() 


## Holdout results

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

reg.models <- list()

reg.models$svm.radial <- trainer.i1$train(reg.model.svm.radial, enable.cache=F)

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




