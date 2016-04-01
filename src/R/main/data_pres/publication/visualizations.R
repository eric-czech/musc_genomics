
GetModelNameMap <- function(){
  c(
    'enet'='Elastic Net',
    'gbm'='Boosted Trees',
    'rf'='Random Forest',
    'knn'='Nearest Neighbors',
    'sda'='Shrinkage Discriminant Analysis',
    'svm'='Radial SVM',
    'xgb'='XGBoost Trees',
    'pls'='Partial Least Squares',
    'pam'='Predictive Analysis for Microarrays',
    'hdrda'='High-Dimensional Regularized DA',
    'scrda'='Shrunken Centroid Regularized DA'
  )
}

GetModelNameTransform <- function(){
  model.map <- GetModelNameMap()
  function(m) ifelse(m %in% names(model.map), model.map[m], m)
}

GeneratePerfProfileVis <- function(
    cv.res.all, response.type, 
    img.dir='~/repos/musc_genomics/src/R/main/data_pres/images/publication/genomics_conf',
    ignore.models=c('nnet')
  ){
  model.trans <- GetModelNameTransform()
  d.plt <- cv.res.all %>% 
    filter(n.feats > 0 & (!model %in% ignore.models)) %>%
    mutate(model=sapply(model, model.trans))
  
  p.acc <- d.plt %>%
    PlotFeatureCountProfile(metric='acc') +
    xlab('Number of Features Included') +
    ylab('Accuracy') +
    ggtitle(sprintf('Accuracy (%s)', toupper(response.type))) + 
    ggsave(sprintf('%s/%s/acc.png', img.dir, response.type))
  
  p.cacc <- d.plt %>%
    PlotFeatureCountProfile(metric='cacc') +
    xlab('Number of Features Included') +
    ylab('Accuracy Over Baseline') +
    ggtitle(sprintf('Accuracy Over Baseline (%s)', toupper(response.type))) + 
    ggsave(sprintf('%s/%s/cacc.png', img.dir, response.type))
  
  p.kappa <- d.plt %>%
    PlotFeatureCountProfile(metric='kappa') +
    xlab('Number of Features Included') +
    ylab('Kappa Statistic') +
    ggtitle(sprintf('Kappa (%s)', toupper(response.type))) + 
    ggsave(sprintf('%s/%s/kappa.png', img.dir, response.type))
  
  list(p.acc, p.cacc, p.kappa)
}
