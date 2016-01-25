#'-----------------------------------------------------------------------------
#' ML Model Training Script
#'
#' This module contains code for training a variety of regression models on 
#' prepared genomics datasets (post feature-selection).
#' 
#' @author eczech
#'-----------------------------------------------------------------------------
source('utils.R')
source('data_prep/prep.R')
source('data_model/training_lib.R')
lib('MASS')
lib('caret')
lib('doMC')
#registerDoMC(8)
SEED <- 1024


d <- GetPreparedData()

fields <- d$fields
c.ge <- fields$gene_expression
c.cn <- fields$copy_number
c.mu <- fields$mutations
d <- d$data

Prep <- function(d, c.numeric, c.binary, min.mutations=5){
  d %>% filter(!is.na(ic_50)) %>%                    # For now we're using COSMIC data only
    rename(response=ic_50) %>% select(-auc) %>%      # Select response field
    Filter(function(x)!all(is.na(x)), .) %>%         # Remove NA-only columns
    RemoveNA(row.threshold=.1) %>%                    # Remove rows with large % NA's
    RemoveRareMutations(c.binary, min.mutations) %>% # Remove cols for rare mutations
    mutate(response=ScaleVector(response))           # Scale response
  # nearZeroVar would make sense here as well but no features would be removed in past tests
}
d.prep <- Prep(d, c(c.ge, c.cn), c.mu)
c.mu <- GetFeatures(d.prep, 'mu')
c.cn <- GetFeatures(d.prep, 'cn')
c.ge <- GetFeatures(d.prep, 'ge')
if (any(is.na(d.prep)))
  stop('Dataset contains unexpected NA values')

# Remove features with low correlation with response
# dt <- ApplyFeatureFilter(d.prep, c(c.ge, c.cn), c.mu, 
#                              numeric.score.p=.0001, binary.score.p=.15)
# numeric.score.p=.025, binary.score.p=.075 --> bin=11, num=9978
# numeric.score.p=.01, binary.score.p=.1 --> bin=17, num=6034
# numeric.score.p=.005, binary.score.p=.15 --> bin=26, num=4116 
# numeric.score.p=.0001, binary.score.p=.15 --> bin=26, num=1007 
# length(GetFeatures(dt, 'mu')); length(GetFeatures(dt, 'cn')) + length(GetFeatures(dt, 'ge'))


# d.prep <- data.frame(
#   tumor_id=c('t1', 't2', 't3', 't4', 't5', 't6'),
#   ic_50=c(-2.3, -1.2, -.3, .9, 1.5, 2.8),
#   auc=rnorm(6),
#   'cn.A1BG'=c(2.1, 1.3, .4, -.7, -1.3, -3),
#   'cn.A1BG.AS1'=c(-1, -2, -5, -3, -4, -3),
#   'ge.A1BG'=rnorm(6),
#   'ge.A1BG.AS1'=c(-6, -3, -1, 1, 3, 3.5),
#   'mu.XPO4_L499F'=c(0, 0, 0, 1, 1, 1),
#   'mu.XPO4_X498.SPLICE'=c(1, 1, 0, 0, 0, 1)
# )
# d.prep <- rbind(d.prep, d.prep)
# d.prep <- rbind(d.prep, d.prep)
# d.prep <- rbind(d.prep, d.prep)
# c.mu <- GetFeatures(d.prep, 'mu')
# c.cn <- GetFeatures(d.prep, 'cn')
# c.ge <- GetFeatures(d.prep, 'ge')
# d.prep <- Prep(d.prep, c(c.ge, c.cn), c.mu)

set.seed(SEED)
d.prep.tr <- d.prep %>% sample_frac(.8, replace = F)
d.prep.ho <- d.prep %>% filter(!tumor_id %in% d.prep.tr$tumor_id)
trctrl <- function(...) trainControl(method = "cv", number = 3, savePredictions='final', ...) # 3
fsctrl <- function(...) sbfControl(method = "cv", number = 5, ...) # 5
preproc <- c('zv', 'center', 'scale')

fs.s <- GetFeatureSelector(.0001, .15)
fs.l <- GetFeatureSelector(.005, .3)
#fs.m <- GetFeatureSelector(.0025, .15)
#fs.l <- GetFeatureSelector(.005, .15)


# sbf(response ~ . - tumor_id, data=d.prep, method='lm', 
#     sbfControl=sbfControl(functions=sbfFunc, method='cv', number=3))

TrainModels <- function(X, y){
  
  
  browser()
  
  folds <- createFolds(Sonar$Class)
  
  # PLS 
  set.seed(SEED); registerDoMC(8)
  m.pls.1 <- sbf(
    X, y, method='pls', preProc = preproc, 
    tuneGrid=data.frame(ncomp=c(1,2,3)),
    trControl = trctrl(), sbfControl = fsctrl(functions = fs.l)
  )
  
  X.pp <- predict(preProcess(X, method=preproc), X)
  init <- glmnet(as.matrix(X.pp), y, 
                 family = 'gaussian', 
                 nlambda = 100, alpha = .01)
  lambda <- unique(init$lambda)
  lambda <- lambda[-c(1, length(lambda))]
  lambda <- lambda[1:length(lambda)]
  
  # GLMNET
  set.seed(SEED); registerDoMC(5)
  m.glmnet.1 <- sbf(
    X, y, method='glmnet', preProc = c(preproc), 
    tuneGrid=expand.grid(.alpha = seq(0, .12, by=.01), .lambda = lambda),
    #tuneGrid=expand.grid(.alpha = 10^seq(-1,-5,length=5), .lambda=10^seq(1,-5,length=100)),
    #tuneLength=25,
    trControl = trctrl(), 
    sbfControl = fsctrl(functions = fs.s),
    standardize=F
  )
  m.glmnet.1.custom <- m.glmnet.1
  save(m.glmnet.1, file='/home/eczech/data/musc_genomics/m.glmnet.1.Rdata')
  save(m.glmnet.2, file='/home/eczech/data/musc_genomics/m.glmnet.2.Rdata')
  
  # GLMNET w/ PCA
  set.seed(SEED)
  m.glmnet.2 <- train(
    X, y, method='glmnet', preProc = c(preproc, 'pca'),
    tuneGrid=expand.grid(.alpha = seq(0, .12, by=.01), .lambda = lambda),
    trControl = trctrl(preProcOptions = list(thresh = 0.3)), 
    type.gaussian="naive", standardize=F
  )
  
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
  
  list(pls1=m.pls.1, pls1=m.pls.2, glmnet1=m.glmnet.1, glmnet2=m.glmnet.2)
}

set.seed(SEED)
#d.samp <- d.prep.tr %>% sample_frac(.2)
#X <- d.samp[,sample(c(c.ge, c.cn, c.mu), replace=F, size = 100)]; y <- d.samp[,'response']
X <- d.prep.tr[,c(c.ge, c.cn, c.mu)]; y <- d.prep.tr[,'response']
m.full <- TrainModels(X, y)

##### Post Feature Selection Models #####
TrainSubsetModels <- function(d){
  # PLS
  set.seed(SEED)
  m.pls.2 <- train(
    response ~ ., data=d.tr.sbf, preProc = c('center', 'scale'), 
    tuneLength=3, trControl = ctrl
  )
  
  # SVM
  set.seed(SEED)
  svmFit <- train(
    response ~ .,
    data = GermanCreditTrain,
    method = "svmRadial",
    preProc = c("center", "scale"),
    tuneLength = 10,
    trControl = trainControl(method = "repeatedcv", repeats = 5))
}



