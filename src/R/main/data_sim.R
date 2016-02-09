library(MASS)
library(glmnet)

set.seed(1)
N <- 100
r <- foreach(p=c(seq(10, 500, length.out=20), seq(250, 2000, by=250)), .combine=rbind)%do%{
  print(p)
  q <- ceiling(p * .1)
  X <- mvrnorm(n=N, rep(0, p), diag(rep(.8, p)) + .2)
  b <- c(rnorm(q), rep(0, p - q))
  y <- X %*% b + rnorm(n=N, sd = 1)  
  m <- cv.glmnet(x = X, y = y, alpha = 1)
  
  # Determine number of features with coef > 0 along 
  # regularization path where lambda = lambda.min in CV
  n.nonzero <- varImp(m$glmnet.fit, lambda=m$lambda.min) %>% filter(Overall > 0) %>% nrow
  data.frame(p=p, n=n.nonzero)
} 
r %>% ggplot(aes(x=p, y=n)) + geom_line() + theme_bw() +
  geom_abline(intercept=0, slope=.1, linetype='dashed')

trainer.i1 <- Trainer(
  cache.dir='/tmp/sim_data', 
  cache.project=paste0(RESPONSE_TYPE, '.all'), 
  seed=SEED
)
trainer.i1$generateFoldIndex(y, CreateFoldIndex)
fold.data.gen <- GetFoldDataGenerator(preproc, F, n.core=8, 
                                      sml.num.p=.0001, lrg.num.p=.01, sml.bin.p=.1, lrg.bin.p=.15)
trainer.i1$generateFoldData(X, y, fold.data.gen, GetDataSummarizer())

