source('utils.R')
source('data_model/training_lib.R')
lib('stringr')
lib('ggplot2')
lib('reshape2')

GetCVScalarStats <- function(cv.res){
  cv.res$fold.summary %>% group_by(model, fold) %>%
    summarise_each(funs(head(., 1)), -model) %>% ungroup
}

GetCVPerfSummary <- function(cv.res){
  GetCVScalarStats(cv.res) %>% group_by(model) %>% 
    summarise_each(funs(mean, sd, lo=quantile(., .25), hi=quantile(., .75)), -model, -x, -y, -t) %>%
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
    geom_pointrange() + facet_wrap(~model) + theme_bw() + 
    xlab('Margin Size (on probability scale)') + 
    ylab('Accuracy Range (IQR)')
}

PlotFoldConfusion <- function(cv.res){
  cv.res$fold.summary %>% group_by(model, fold) %>% do({head(., 1)}) %>% ungroup %>%
    select(model, fold, len, starts_with('cm.')) %>% group_by(model) %>% do({
      d <- .
      pcts <- d %>% select(starts_with('cm.'))
      pcts %>% setNames(str_replace_all(names(pcts), 'cm\\.', ''))
    }) %>% melt(id.vars='model') %>% 
    ggplot(aes(x=model, y=value, color=model)) + 
    geom_boxplot(alpha=.5, outlier.size=0, position = 'dodge', width = 0.5) +
    geom_jitter(width = .3, alpha=.5) + 
    facet_wrap(~variable, scales='free') + theme_bw() + 
    theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
    scale_colour_discrete(guide = FALSE) +
    xlab('Model') + ylab('Count') + 
    ggtitle('Confusion Matrix Distributions')
}

PlotFoldMetric <- function(cv.res, metric){
  GetCVScalarStats(cv.res) %>% 
    rename_(value=metric) %>%
    #mutate(model=factor(model)) %>%
    #mutate(model=reorder(model, value, FUN = median)) %>%
    ggplot(aes(x=model, y=value, color=model)) +
    geom_boxplot() + coord_flip() + theme_bw() + ylab(metric) + 
    scale_colour_discrete(guide = FALSE) +
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

PlotPerFoldPR <- function(cv.res){
  cv.res$fold.summary %>% 
    select(model, fold, x, y, t) %>%
    inner_join(GetCVPerfSummary(cv.res), by='model') %>% ungroup %>%
    filter(str_detect(model, '.')) %>% mutate(
      model.label=sprintf('%s (sens=%s, kappa=%s)', model, round(sens_mean, 2), round(kappa_mean, 2)),
      fold=factor(fold)
    ) %>% ggplot(aes(x=x, y=y, color=fold)) + 
    geom_line() + theme_bw() + 
    facet_wrap(~model.label) + xlab('Recall') + ylab('Precision') +
    ggtitle('CV PR')
}

PlotAllFoldROC <- function(cv.res){
  cv.res$model.summary %>%
    select(model, x, y, t) %>%
    inner_join(GetCVPerfSummary(cv.res), by='model') %>%
    mutate(model.label=paste0(model, ' (', round(auc_mean, 3), ')')) %>% 
    ggplot(aes(x=x, y=y, color=factor(model.label))) + 
    geom_abline(alpha=.5) + geom_line() + theme_bw() +
    xlab('False Positive Rate') + ylab('True Positive Rate') +
    ggtitle('CV ROC (combined)')
}

PlotAllFoldPR <- function(cv.res){
  cv.res$model.summary %>%
    select(model, x, y, t) %>%
    inner_join(GetCVPerfSummary(cv.res), by='model') %>%
    mutate(model.label=sprintf('%s (sens=%s, kappa=%s)', model, round(sens_mean, 2), round(kappa_mean, 2))) %>%
    ggplot(aes(x=x, y=y, color=factor(model.label))) + 
    geom_line() + theme_bw() +
    facet_wrap(~model.label) + xlab('Recall') + ylab('Precision') +
    ggtitle('CV PR (combined)')
}


PlotHoldOutConfusion <- function(ho.res){
  ho.res$model.summary %>% group_by(model) %>% do({head(., 1)}) %>% ungroup %>%
    select(model, len, starts_with('cm.')) %>% group_by(model) %>% do({
      d <- .
      pcts <- d %>% select(starts_with('cm.'))
      pcts %>% setNames(str_replace_all(names(pcts), 'cm\\.', ''))
    }) %>% melt(id.vars='model') %>% 
    ggplot(aes(x=model, y=value, fill=model)) + geom_bar(stat='identity') +
    facet_wrap(~variable, scales='free') + theme_bw() + 
    xlab('Model') + ylab('Count') + 
    ggtitle('Confusion Matrix Distributions')
}

PlotHoldOutMetric <- function(ho.res, metric){
  ho.res$model.summary %>% group_by(model) %>% do({head(., 1)}) %>%
    select(one_of(metric), model) %>% melt(id.vars = 'model') %>%
    arrange(value) %>% mutate(model=factor(model, levels=model)) %>%
    ggplot(aes(x=model, y=value, color=model)) + ylab(metric) + 
    geom_point(stat='identity', size=3) + theme_bw() + coord_flip() + 
    ggtitle(sprintf('Holdout %s', metric))
}

PlotHoldOutPR <- function(ho.res){
  ho.res$model.summary %>%
    select(model, x, y, t) %>%
    inner_join(ho.res$model.summary %>% group_by(model) %>% summarise(sens=sens[1], kappa=kappa[1]), by='model') %>% 
    mutate(model.label=sprintf('%s (sens=%s, kappa=%s)', model, round(sens, 2), round(kappa, 2))) %>%
    ggplot(aes(x=x, y=y, color=factor(model.label))) + 
    geom_line() + theme_bw() + 
    facet_wrap(~model.label) + xlab('Recall') + ylab('Precision') +
    ggtitle('Holdout PR')
}

PlotHoldOutLift <- function(ho.res){
  ho.res$model.summary %>%
    select(model, x, y, t) %>%
    inner_join(ho.res$model.summary %>% group_by(model) %>% summarise(sens=sens[1], kappa=kappa[1]), by='model') %>% 
    mutate(model.label=sprintf('%s (sens=%s, kappa=%s)', model, round(sens, 2), round(kappa, 2))) %>%
    ggplot(aes(x=x, y=y, color=factor(model.label))) + 
    geom_line() + theme_bw() + 
    facet_wrap(~model.label) + xlab('% Cases Predicted') + ylab('Lift') +
    ggtitle('Holdout Lift')
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

PlotResponseDist <- function(response.type, response.data, response.cutoff=0){
  data.frame(r=response.data) %>% ggplot(aes(x=r)) + 
    geom_histogram(binwidth=.15, alpha=.3) + 
    geom_density(aes(y=.15 * ..count..)) + theme_bw() + 
    geom_vline(aes(xintercept=response.cutoff)) + 
    ggtitle(sprintf('%s Response Distribution', response.type)) 
}

