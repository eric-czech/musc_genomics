# Data Modeling

## Feature Selection Benchmarking

For the genomics dataset (w/ COSMIC IC 50 response) having dimensions ~ 400 obs x 47k features, here are some recorded times necessary for removing features uncorrelated with response:

  dt <- ApplyFeatureFilter(d.prep, c(c.ge, c.cn), c.mu, 
                                numeric.score.p=.0001, binary.score.p=.15)
  # numeric.score.p=.025, binary.score.p=.075 --> bin=11, num=9978
  # numeric.score.p=.01, binary.score.p=.1 --> bin=17, num=6034
  # numeric.score.p=.005, binary.score.p=.15 --> bin=26, num=4116 
  # numeric.score.p=.0001, binary.score.p=.15 --> bin=26, num=1007 
  # length(GetFeatures(dt, 'mu')); length(GetFeatures(dt, 'cn')) + length(GetFeatures(dt, 'ge'))