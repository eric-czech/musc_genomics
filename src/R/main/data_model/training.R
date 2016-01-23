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
lib('MASS')
lib('caret')
SEED <- 1024

d <- GetPreparedData()

fields <- d$fields
c.ge <- fields$gene_expression
c.cn <- fields$copy_number
c.mu <- fields$mutations
d <- d$data

Prep <- function(d, c.numeric, c.binary, min.mutations=3){
  d %>% filter(!is.na(ic_50)) %>%                    # For now, we're using COSMIC data
    rename(response=ic_50) %>% select(-auc) %>%      # Select response field
    Filter(function(x)!all(is.na(x)), .) %>%         # Remove NA-only columns
    RemoveRareMutations(c.binary, min.mutations) %>% # Remove cols for rare mutations
    mutate(response=scale(response))                 # Scale response
  # Apply nearZeroVar() as well?
}
d.prep <- Prep(d, c(c.ge, c.cn), c.mu)
c.mu <- GetFeatures(d.prep, 'mu')

d.prep <- data.frame(
  tumor_id=c('t1', 't2', 't3', 't4', 't5', 't6'),
  ic_50=c(-2.3, -1.2, -.3, .9, 1.5, 2.8),
  auc=rnorm(6),
  'cn.A1BG'=c(2.1, 1.3, .4, -.7, -1.3, -3),
  'cn.A1BG.AS1'=c(-1, -2, -5, -3, -4, -3),
  'ge.A1BG'=rnorm(6),
  'ge.A1BG.AS1'=c(-6, -3, -1, 1, 3, 3.5),
  'mu.XPO4_L499F'=c(0, 0, 0, 1, 1, 1),
  'mu.XPO4_X498.SPLICE'=c(1, 1, 0, 0, 0, 1)
)
d.prep <- rbind(d.prep, d.prep)
d.prep <- rbind(d.prep, d.prep)
d.prep <- rbind(d.prep, d.prep)
c.mu <- GetFeatures(d.prep, 'mu')
c.cn <- GetFeatures(d.prep, 'cn')
c.ge <- GetFeatures(d.prep, 'ge')
d.prep <- Prep(d.prep, c(c.ge, c.cn), c.mu)



set.seed(SEED)
d.prep.tr <- d.prep %>% sample_frac(.8, replace = F)
d.prep.ho <- d.prep %>% filter(!tumor_id %in% d.prep.tr$tumor_id)
trctrl <- function(...) trainControl(method = "cv", number = 3, ...)
fsctrl <- function(...) sbfControl(method = "cv", number = 3, ...)
form <- as.formula('response ~ . - tumor_id')

fs.s <- GetFeatureSelector(.05)
fs.m <- GetFeatureSelector(.1)
fs.l <- GetFeatureSelector(.25)

# sbf(response ~ . - tumor_id, data=d.prep, method='lm', 
#     sbfControl=sbfControl(functions=sbfFunc, method='cv', number=3))

TrainModels <- function(d){
  
  # PLS - No feature selection
  set.seed(SEED)
  m.pls.1 <- train(
    form, data=d, method='pls', preProc = c('center', 'scale'), 
    tuneGrid=data.frame(ncomp=2:4), trControl = trctrl()
  )
  
  # PLS - With feature selection
  set.seed(SEED)
  m.pls.2 <- sbf(
    form, data=d, method='pls', preProc = c('center', 'scale'), 
    tuneGrid=data.frame(ncomp=2:4), trControl = trctrl(), 
    sbfControl = fsctrl(functions = fs.s)
  )
  
  # GLMNET - No feature selection
  set.seed(SEED)
  m.glmnet.1 <- train(
    form, data=d, method='glmnet', preProc = c('center', 'scale'), 
    tuneGrid=expand.grid(.alpha = seq(.05, 1, length = 10), .lambda=seq(0,0.05,by=0.01)),
    trControl = trctrl()
  )
  
  # GLMNET - With feature selection
  set.seed(SEED)
  m.glmnet.2 <- sbf(
    form, data=d, method='glmnet', preProc = c('center', 'scale'), 
    tuneGrid=expand.grid(.alpha = seq(.05, 1, length = 10), .lambda=seq(0,0.05,by=0.01)),
    trControl = trctrl(), sbfControl = fsctrl(functions = fs.l)
  )
  
  list(pls1=m.pls.1, pls1=m.pls.2, glmnet1=m.glmnet.1, glmnet2=m.glmnet.2)
}

m.full <- TrainModels(d.prep.tr)

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



