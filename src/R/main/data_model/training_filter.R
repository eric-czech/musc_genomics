FeatureScore <- function(x, y) {
  if (length(unique(x)) == 2){
    pv <- try(fisher.test(factor(x), y)$p.value, silent = TRUE)
  } else { 
    #pv <- try(anova(lm(x ~ y), test = "F")[1, "Pr(>F)"], silent = TRUE)
    pv <- try(t.test(x ~ y)$p.value, silent = TRUE)
  }
  if (any(class(pv) == "try-error") || is.na(pv) || is.nan(pv)) pv <- 1
  pv
}

GetLimitFilter <- function(seed, limit, model.args, verbose=T){
  .FeatureScore <- FeatureScore
  list(
    summary = defaultSummary,
    fit = function(x, y, ...){
      train.args <- model.args
      train.args$x <- x
      train.args$y <- y
      if (verbose) cat('(+ sbf) Training model within filter\n')
      set.seed(seed)
      r <- do.call('train', train.args)
      if (verbose) cat('(- sbf) Model training complete\n')
      r
    },
    pred = function(object, x) {
      if (verbose) cat('(+ sbf) Making predictions within filter\n')
      set.seed(seed)
      pred.class <- data.frame(pred=predict(object, newdata=x, type='raw'))
      set.seed(seed)
      pred.prob <- as.data.frame(predict(object, newdata=x, type = "prob"))
      if (verbose) cat('(- sbf) Predictions complete\n')
      cbind(pred.class, pred.prob)
    },
    score = function(x, y){
      ## should return a named logical vector
      if(!is.factor(y))
        stop('y must be factor')
      .FeatureScore(x, y)
    },
    filter = function(score, x, y) {
      if (verbose) cat(sprintf('(- sbf) Sorting scored feature vector of length %s\n', length(score)))
      score <- sort(score, decreasing=F)
      n.keep <- min(limit, length(score))
      if (verbose) cat('(- sbf) Feature sort complete\n')
      setNames(c(rep(T, n.keep), rep(F, length(score) - n.keep)), names(score))
    }
  )
}

# d <- twoClassSim(n=500, noiseVars=100)
# X <- d %>% select(-Class)
# y <- d$Class
# 
# model.args <- list(
#   method='xgbTree',
#   metric='Kappa', 
#   preProcess=c('zv', 'center', 'scale'),
#   tuneLength=3,
#   trControl=trainControl(
#     method='cv', number=10, classProbs=T, returnData=F,
#     savePredictions='final'
#   )
# )
# 
# registerDoMC(1)
# modelFilter <- GetLimitFilter(123, 5, model.args)
# sbfctrl <- sbfControl(functions=modelFilter, method='cv', number=5, saveDetails=T, verbose=T, allowParallel=F)
# m <- sbf(X, y, sbfControl = sbfctrl)

