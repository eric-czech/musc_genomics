#'-----------------------------------------------------------------------------
#' ML Model Training Script
#'
#' This module contains code for training a variety of regression models on 
#' prepared genomics datasets (post feature-selection).
#' 
#' @author eczech
#'-----------------------------------------------------------------------------
options( java.parameters = "-Xmx4g" )
source('utils.R')
source('data_model/training_lib.R')
source('data_model/training_models.R')
source('data_model/training_filter.R')
source('data_model/training_viz.R')
source('data_model/filtering_lib.R')
source('~/repos/portfolio/functional/ml/R/trainer.R')
source('~/repos/portfolio/functional/ml/R/results.R')
lib('MASS')
lib('caret')
lib('doMC')
lib('iterators')
lib('caretEnsemble') # library(devtools); unload(inst("caretEnsemble")); install_local('/home/eczech/repos/misc/caretEnsemble');
lib('ROCR')
lib('plotly')
SEED <- 1024

## Choose dataset to use for modeling (must pick one of the following)
EnableCosmic()
#EnableCtd()

PREPROC <- c('zv', 'center', 'scale')

RESULT_CACHE <- Cache(dir=file.path(CACHE_DIR, 'result_data'), project=RESPONSE_TYPE)
select <- dplyr::select

d.prep <- GetTrainingData(TRAIN_CACHE, RESPONSE_TYPE, RESPONSE_SELECTOR, min.mutations=3)

##### Data Partitioning #####

set.seed(SEED)
idx.tr <- createDataPartition(d.prep[,'response'], p=.8)[[1]]

split.data <- function(data, idx, type, N){
  d <- data[idx,]; X <- d %>% select(-response, -tumor_id)
  y <- d[,'response']; y.bin <- DichotomizeOutcome(y, threshold = RESPONSE_THRESH)
  n <- length(y); n.pos <- sum(y.bin == 'pos')
  summary <- data.frame(n=n, pct.of.total=n/N, pos=n.pos, pos.pct=n.pos/n, type=type)
  list(X=X, y=y, y.bin=y.bin, summary=summary)
}
d.tr <- split.data(d.prep, idx.tr, 'training', nrow(d.prep))
d.ho <- split.data(d.prep, -idx.tr, 'holdout', nrow(d.prep))

##### Model Trainer #####

trainer <- Trainer(cache.dir=file.path(CACHE_DIR, 'filtering_data'), 
                   cache.project=RESPONSE_TYPE, seed=SEED)
trainer$generateFoldIndex(d.tr$y, CreateFoldIndex)
fold.data.gen <- GetFeatScoringFoldGen(PREPROC, RESPONSE_THRESH, feat.limit=5000, n.core=2)
trainer$generateFoldData(d.tr$X, d.tr$y, fold.data.gen, FilteringDataSummarizer())

GetPrepFun <- function(origin.transform=TransformOriginSolidLiquid){
  function(X){
    X %>% mutate(origin=origin.transform(origin))
  }
}

GetFilterModel <- function(name, max.feats, n.core=3, 
                           origin.transform=TransformOriginSolidLiquid, ...){
  prep.fun <- GetPrepFun(origin.transform)
  pred.fun <- function(fit, d, i){ 
    X <- d$X.test[, fit$finalModel$xNames] %>% prep.fun
    list(
      prob=predict(fit, X, type='prob')[,1], 
      class=predict(fit, X, type='raw')
    )
  }
  list(
    name=name, predict=pred.fun, test=bin.test,
    train=function(d, idx, i){
      registerDoMC(n.core)
      # d$feat.scores contains scores for all 36k features
      top.feats <- d$feat.scores %>% arrange(score) %>% head(max.feats) %>% .$feature 
      X <- prep.fun(d$X.train[,top.feats])
      y <- d$y.train.bin
      train(X, y, ...)
    }
  )
}

GetXGBoostModel <- function(max.feats, n.core=3, tuneLength=8, k=5){
  name <- sprintf('xgb.%s', max.feats)
  GetFilterModel(
    name, max.feats, n.core=n.core,
    method='xgbTree', metric='Accuracy', 
    preProcess=c('zv'), tuneLength=tuneLength,
    trControl=trainControl(
      method='cv', number=k, classProbs=T, 
      returnData=F, savePredictions='final',
      allowParallel=T, verboseIter=F
    )
  )
}

model.xgb.15 <- GetXGBoostModel(15, n.core=3)

models <- list()
ec <- F
models$xgb.15 <- trainer$train(model.xgb.15, enable.cache=ec)

cv.res <- SummarizeTrainingResults(
  models, T, fold.summary=ResSummaryFun('roc'), model.summary=ResSummaryFun('roc'))

##### Simple Filter Models #####

# source('~/repos/portfolio/functional/ml/R/parallel.R')
# registerCores(3, log.file = '/tmp/musc_genomics_training.log')
registerDoMC(3)

RunXGBFilter <- function(
  X, y, max.feats=50, tune.length=10,
  outer.method='repeatedcv', outer.number=10, outer.repeats=1,
  inner.method='cv', inner.number=8){
  seed <- 1024
  GetXGBFilter <- function(limit){
    model.args <- list(
      method='xgbTree',
      metric='Accuracy', 
      preProcess=c('zv'),
      tuneLength=tune.length,
      trControl=trainControl(
        method=inner.method, number=inner.number, classProbs=T, 
        returnData=T, savePredictions='final',
        allowParallel=F, verboseIter=T
      )
    )
    GetLimitFilter(seed, limit, model.args)
  }
  set.seed(seed)
  sbfctrl <- sbfControl(
    functions=GetXGBFilter(max.feats), 
    method=outer.method, number=outer.number, repeats=outer.repeats, 
    saveDetails=T, verbose=T, allowParallel=T)
  set.seed(seed)
  sbf(X, y, sbfControl = sbfctrl)
}


set.seed(SEED)
d.samp <- d.prep %>% sample_frac(1)
# d.samp <- d.prep
X <- d.samp %>% select(-tumor_id, -response, -origin) %>%
  mutate(origin=ifelse(toupper(origin) == 'HAEMATOPOIETIC_AND_LYMPHOID_TISSUE', 1, 0))
y <- DichotomizeOutcome(d.samp[,'response'], threshold=RESPONSE_THRESH)
# table(x=X[,'origin'], y=y)
# fisher.test(X[,'origin'], y)$p.value

# Benchmarking runs
by.max.feats <- lapply(seq(5, 50, by=5), function(i){
  cat(sprintf('Running model for max feat %s\n', i))
  t.start <- proc.time()
  res.xgb <- RunXGBFilter(
    X, y, outer.method='cv', outer.number=10, outer.repeats=1,
    inner.method='cv', inner.number=5, max.feats=i, tune.length=8)
  time.diff <- proc.time() - t.start  
  list(max.feat=i, time=time.diff, res=res.xgb)
})
res.xgb.w.o <- res.xgb

data.frame(r1=res.xgb$resample$Accuracy, r2=res.xgb.w.o$resample$Accuracy) %>%
  melt(id.vars=NULL) %>%
  ggplot(aes(x=variable, y=value)) + geom_boxplot()

# 298.063 seconds for sample frac .3 (123 rows)
# 377 seconds for sample frac .3
# 311 seconds for sample frac 1, w/ domc(3)
# 629 seconds for sample frac 1, w/ domc(3)

# Actual run
# res.xgb <- RunXGBFilter(X, y)

