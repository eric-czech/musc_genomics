
GetCVPerfSummary <- function(cv.res){
  cv.res$fold.summary %>% group_by(model, fold) %>%
    summarise_each(funs(head(., 1)), one_of(c('auc', 'acc.max', 'acc.cut'))) %>%
    ungroup %>% group_by(model) %>% 
    summarise_each(funs(mean, sd), one_of(c('auc', 'acc.max', 'acc.cut'))) %>%
    ungroup
}

PlotFoldAUC <- function(cv.res){
  GetCVPerfSummary(cv.res) %>% arrange(auc_mean) %>%
    mutate(model=factor(model, levels=model)) %>%
    ggplot(aes(x=model, y=auc_mean, ymin=auc_mean - auc_sd, ymax=auc_mean + auc_sd, color=model)) +
    geom_pointrange() + coord_flip() + theme_bw() +
    ggtitle('CV AUC Estimates')
}

PlotPerFoldROC <- function(cv.res){
  cv.res$fold.summary %>% 
    select(model, fold, x, y, t) %>%
    inner_join(GetCVPerfSummary(cv.res), by='model') %>% ungroup %>%
    filter(str_detect(model, '.')) %>% mutate(
      model.label=sprintf('%s (auc=%s, acc=%s%%)', model, round(auc_mean, 2), round(acc.max_mean*100, 0)),
      fold=factor(fold)
    ) %>% ggplot(aes(x=x, y=y, color=fold)) + 
    geom_abline(alpha=.5) + geom_line() + theme_bw() + 
    facet_wrap(~model.label) + xlab('False Positive Rate') + ylab('True Positive Rate') +
    ggtitle('CV ROC')
}

PlotAllFoldROC <- function(cv.res){
  cv.res$model.summary %>%
    select(model, x, y, t) %>%
    inner_join(model.perf, by='model') %>%
    filter(str_detect(model, '.')) %>% 
    mutate(model.label=paste0(model, ' (', round(auc_mean, 3), ')')) %>% 
    ggplot(aes(x=x, y=y, color=factor(model.label))) + 
    geom_abline(alpha=.5) + geom_line() + theme_bw() +
    xlab('False Positive Rate') + ylab('True Positive Rate') +
    ggtitle('CV ROC (combined)')
}

PlotHoldOutMetric <- function(ho.res, metric){
  ho.res %>% group_by(model) %>% do({head(., 1)}) %>%
    select(one_of(metric), model) %>% melt(id.vars = 'model') %>%
    arrange(desc(value)) %>% mutate(model=factor(model, levels=model)) %>%
    ggplot(aes(x=model, y=value, color=model)) + 
    geom_point(stat='identity', size=3) + theme_bw()
}

