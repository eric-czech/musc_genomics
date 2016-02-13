source('utils.R')
source('data_model/training_lib.R')
lib('stringr')
lib('ggplot2')
lib('reshape2')

GetCVPerfSummary <- function(cv.res){
  cv.res$fold.summary %>% group_by(model, fold) %>%
    summarise_each(funs(head(., 1)), -model) %>%
    ungroup %>% group_by(model) %>% 
    summarise_each(funs(mean, sd, lo=quantile(., .25), hi=quantile(., .75)), -model) %>%
    ungroup
}

PlotFoldMetricByMargin <- function(cv.res, metric){
  d <- GetCVPerfSummary(cv.res) 
  metrics <- d %>% select(matches(sprintf('%s_margin_.*_mean$', metric))) %>% names %>% str_replace('_mean$', '')
  d <- d %>% select(model, 
                    one_of(paste0(metrics, '_mean')), 
                    one_of(paste0(metrics, '_lo')),
                    one_of(paste0(metrics, '_hi'))) %>% data.frame
  d <- foreach(m=metrics, .combine=cbind)%do%{
    r <- data.frame(
      lo=d[,paste0(m, '_lo')],
      hi=d[,paste0(m, '_hi')],
      mid=d[,paste0(m, '_mean')]
    )
    setNames(r, paste(m, c('lo', 'hi', 'mid'), sep='_'))
  } %>% cbind(data.frame(model=d$model))
  
  d %>% melt(id.vars='model') %>% #.$variable %>% as.character %>% str_split('_')
    mutate(variable=as.character(variable)) %>%
    mutate(margin=sapply(str_split(variable, '_'), function(x) as.numeric(x[3]))) %>%
    mutate(range=sapply(str_split(variable, '_'), function(x) x[4])) %>%
    select(-variable) %>% dcast(model + margin ~ range) %>%
    ggplot(aes(x=factor(margin), y=mid, ymax=hi, ymin=lo)) + 
    geom_pointrange() + facet_wrap(~model) + theme_bw()
}

PlotFoldMetric <- function(cv.res, metric){
  GetCVPerfSummary(cv.res) %>% 
    rename_(v_mean=paste0(metric, '_mean'), v_sd=paste0(metric, '_sd')) %>%
    arrange(v_mean) %>% mutate(model=factor(model, levels=model)) %>%
    ggplot(aes(x=model, y=v_mean, ymin=v_mean - v_sd, ymax=v_mean + v_sd, color=model)) +
    geom_pointrange() + coord_flip() + theme_bw() + ylab(metric) + 
    ggtitle(sprintf('CV %s Estimates', metric))
}

PlotPerFoldROC <- function(cv.res){
  cv.res$fold.summary %>% 
    select(model, fold, x, y, t) %>%
    inner_join(GetCVPerfSummary(cv.res), by='model') %>% ungroup %>%
    filter(str_detect(model, '.')) %>% mutate(
      model.label=sprintf('%s (auc=%s, acc=%s%%)', model, round(auc_mean, 2), round(acc_mean*100, 0)),
      fold=factor(fold)
    ) %>% ggplot(aes(x=x, y=y, color=fold)) + 
    geom_abline(alpha=.5) + geom_line() + theme_bw() + 
    facet_wrap(~model.label) + xlab('False Positive Rate') + ylab('True Positive Rate') +
    ggtitle('CV ROC')
}

PlotAllFoldROC <- function(cv.res){
  cv.res$model.summary %>%
    select(model, x, y, t) %>%
    inner_join(GetCVPerfSummary(cv.res), by='model') %>%
    filter(str_detect(model, '.')) %>% 
    mutate(model.label=paste0(model, ' (', round(auc_mean, 3), ')')) %>% 
    ggplot(aes(x=x, y=y, color=factor(model.label))) + 
    geom_abline(alpha=.5) + geom_line() + theme_bw() +
    xlab('False Positive Rate') + ylab('True Positive Rate') +
    ggtitle('CV ROC (combined)')
}

PlotHoldOutMetric <- function(ho.res, metric){
  ho.res$model.summary %>% group_by(model) %>% do({head(., 1)}) %>%
    select(one_of(metric), model) %>% melt(id.vars = 'model') %>%
    arrange(value) %>% mutate(model=factor(model, levels=model)) %>%
    ggplot(aes(x=model, y=value, color=model)) + ylab(metric) + 
    geom_point(stat='identity', size=3) + theme_bw() + coord_flip() + 
    ggtitle(sprintf('Holdout %s', metric))
}

PlotHoldOutROC <- function(ho.res){
  ho.res$model.summary %>%
    select(model, x, y, t) %>%
    inner_join(ho.res$model.summary %>% group_by(model) %>% summarise(auc=auc[1]), by='model') %>% 
    mutate(model.label=paste0(model, ' (', round(auc, 3), ')')) %>% 
    ggplot(aes(x=x, y=y, color=factor(model.label))) + 
    geom_line() + geom_abline(alpha=.5) + theme_bw() + 
    xlab('False Positive Rate') + ylab('True Positive Rate') +
    ggtitle('Holdout ROC')
    
}

PlotResponseDist <- function(response.type, response.data){
  data.frame(r=response.data) %>% ggplot(aes(x=r)) + 
    geom_histogram(binwidth=.15, alpha=.3) + 
    geom_density(aes(y=0.15*..count..)) + theme_bw() + 
    ggtitle(sprintf('%s Response Distribution', response.type))
}

