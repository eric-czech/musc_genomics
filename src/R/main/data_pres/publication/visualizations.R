

GetModelLineTypeTransform <- function(){
  function(x){
    if (x %in% c('enet', 'sda', 'pls')) 1
    else if (x %in% c('rf', 'gbm', 'xgb')) 2
    else 3
  }
}

GeneratePerfProfileVis <- function(
    cv.res.all, response.type, 
    img.dir='~/repos/musc_genomics/src/R/main/data_pres/images/publication/genomics_conf',
    ignore.models=c('nnet'), scale=1.5, dpi=600
  ){
  
  d.plt <- cv.res.all %>% 
    filter(n.feats > 0 & (!model %in% ignore.models)) %>%
    mutate(linetype=sapply(model, GetModelLineTypeTransform())) %>%
    mutate(model=sapply(model, GetModelNameTransform()))
  
  model.levels <- unique(d.plt$model)
  ens.levels <- c('Ensemble (Average)', 'Ensemble (GLM)')
  first.levels <- model.levels[model.levels %in% ens.levels]
  rest.levels <- model.levels[!model.levels %in% ens.levels]
  ordered.levels <- c(first.levels, rest.levels)
  
  d.plt <- d.plt %>% mutate(model=factor(as.character(model), levels=ordered.levels))
  
  p.acc <- d.plt %>%
    PlotFeatureCountProfile(metric='acc') +
    xlab('Number of Features Included') +
    ylab('Accuracy') +
    ggtitle(sprintf('Accuracy (%s)', toupper(response.type))) + 
    ggsave(sprintf('%s/%s/acc.png', img.dir, response.type), scale=scale, dpi=dpi)
  
  p.cacc <- d.plt %>%
    PlotFeatureCountProfile(metric='cacc') +
    xlab('Number of Features Included') +
    ylab('Accuracy Over Baseline') +
    ggtitle(sprintf('Accuracy Over Baseline (%s)', toupper(response.type))) + 
    ggsave(sprintf('%s/%s/cacc.png', img.dir, response.type), scale=scale, dpi=dpi)
  
  p.kappa <- d.plt %>%
    PlotFeatureCountProfile(metric='kappa') +
    xlab('Number of Features Included') +
    ylab('Kappa Statistic') +
    ggtitle(sprintf('Kappa (%s)', toupper(response.type))) + 
    ggsave(sprintf('%s/%s/kappa.png', img.dir, response.type), scale=scale, dpi=dpi)
  
  list(p.acc, p.cacc, p.kappa)
}


GenerateErrorProfileVis <- function(cv.res.all, model.names=c('svm', 'enet')){
  cv.res.all %>% filter(n.feats==20 & model %in% model.names) %>% 
    mutate(model=toupper(model)) %>%
    ggplot(aes(x=kappa, fill=model)) + geom_density(alpha=.3) + 
    theme_bw() +
    xlab('Predictive Accuracy (Kappa)') + ylab('Density') +
    scale_fill_brewer(guide=guide_legend(title='Algorithm'), palette="Set1")
}

GenerateTopFeatVis <- function(top.feat.cv){
  top.feat.cv.names <- top.feat.cv %>% filter(frequency >= 18 & feature != 'origin') %>% .$feature %>% as.character 
  d.prep[,c(top.feat.cv.names, 'response')] %>%
    mutate(response=DichotomizeOutcome(response, RESPONSE_THRESH)) %>%
    mutate(response=ifelse(response == 'pos', 'Yes', 'No')) %>%
    melt(id.vars=c('response')) %>% 
    ggplot(aes(x=variable, y=value, color=response)) + geom_boxplot(position='dodge') +
    theme_bw() + scale_color_brewer(guide=guide_legend(title='Was Sensitive'), palette="Set1") +
    xlab('Feature Name') + ylab('Feature Value') 
}

