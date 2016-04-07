source('utils.R')
source('data_model/training_lib.R')
lib('stringr')
lib('ggplot2')
lib('reshape2')
lib('RColorBrewer')



GetModelNameMap <- function(){
  c(
    'enet'='Elastic Net',
    'enetwt'='Weighted Elastic Net',
    'mars'='MV Regression Splines',
    'gbm'='Boosted Trees',
    'rf'='Random Forest',
    'knn'='Nearest Neighbors',
    'sda'='Shrinkage Discriminant Analysis',
    'svm'='Radial SVM',
    'xgb'='XGBoost Trees',
    'pls'='Partial Least Squares',
    'pam'='Predictive Analysis for Microarrays',
    'hdrda'='High-Dimensional Regularized DA',
    'scrda'='Shrunken Centroid Regularized DA',
    'ensglm'='Ensemble (GLM)',
    'ensavg'='Ensemble (Average)'
  )
}

GetModelNameTransform <- function(){
  model.map <- GetModelNameMap()
  function(m) ifelse(m %in% names(model.map), model.map[m], m)
}

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

PlotFoldMetric <- function(cv.res, metric, order.by.score=T, extract.model.names=F){
  p <- GetCVScalarStats(cv.res) %>% 
    rename_(value=metric)
  if (extract.model.names){
    trans <- GetModelNameTransform()
    p <- p %>%
      mutate(model=trans(str_extract(model, '^.*?(?=\\.)')))
  }
  if (order.by.score){
    p <- p %>% mutate(model=factor(model)) %>%
      mutate(model=reorder(model, value, FUN = median))
  }
  p <- p %>%
    ggplot(aes(x=model, y=value, color=model)) +
    geom_boxplot() + coord_flip() + theme_bw() + 
    ylab(metric) + xlab('Model') +
    scale_colour_discrete(guide = FALSE) +
    ggtitle(sprintf('CV %s Estimates', metric))
  p
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
    theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
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


#' @title Aggregate CV results across models and feature subset sizes
#' @return data frame containing model name, number of features, and performance metrics
GetAggregateFilterCVRes <- function(models, metrics=c('kappa', 'cacc', 'acc')){
  cv.res.all <- foreach(m=names(models), .combine=rbind)%dopar%{
    loginfo('Processing model type %s', m)
    cv.res <- SummarizeTrainingResults(
      models[[m]], T, fold.summary=ResSummaryFun('roc'), model.summary=ResSummaryFun('roc')
    )
    GetCVScalarStats(cv.res) %>% select(one_of(metrics), fold, model) %>%
      rename(model.name=model) %>%
      mutate(
        n.feats=as.numeric(str_extract(model.name, '(?<=\\.).*?(?=\\.)')), 
        model=str_extract(model.name, '.*?(?=\\.)')
      )
  }
  n.fold <- length(models[[1]][[1]])
  is.valid <- cv.res.all %>% group_by(model, n.feats) %>% tally %>% .$n %>% unique == n.fold
  if (!is.valid) stop(sprintf(
    'CV result aggregation produced incorrect number of fold results per model (expecting %s)', n.fold
  ))
  cv.res.all
}


#' @title Plots performance profile over number of features included in an arbitrary model set
#' @param cv.res.agg results from \link{GetAggregateFilterCVRes}
#' @param metric performance metric to plot profile over
#' @return ggplot 
#' @seealso 
PlotFeatureCountProfile <- function(cv.res.agg, metric){
  
  n.feat.all <- sort(unique(cv.res.agg$n.feats))
  n.model <- length(unique(cv.res.agg$model))
  n.solid <- floor(n.model / 2)
  n.dash <- n.model - n.solid
  legend.guide <- guide_legend(title='Model Name')
  
  cv.res.agg %>% 
    rename_(value=metric) %>%
    group_by(model, n.feats) %>% mutate(mean_value=mean(value)) %>% ungroup %>%
    mutate(i.feats=factor(n.feats, ordered=T)) %>% 
    ggplot + 
    geom_line(aes(x=as.integer(i.feats), y=mean_value, color=model, linetype=model), lwd=.75, alpha=.5) +
    geom_smooth(aes(x=as.integer(i.feats), y=mean_value), color='black', method='loess', level=1-1E-9) +
    scale_x_continuous(labels=as.character(n.feat.all), breaks=1:length(n.feat.all)) +
    scale_color_manual(
      values = c(brewer.pal(n.solid, "Set1"), brewer.pal(n.dash, "Set1")),
      guide=legend.guide
    ) + scale_linetype_manual(
      values = c(rep("solid", n.solid), rep("dotted", n.dash)),
      guide=legend.guide
    ) +
    theme_bw() + theme(
      panel.grid.minor=element_blank(), 
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(size=.1)
    )
    #theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())
}

PlotHoldoutGainFinal <- function(d.ho.lift){
  max.x <- max(d.ho.lift$x)
  d.ho.lift %>%
    ggplot(aes(x=x, y=y)) + geom_line(color='#336633', size=1, alpha=.5) +
    theme_bw() + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=.5)) +
    xlab('Number of Predictions Used') +
    ylab('Number of Actually Sensitive Samples') +
    geom_abline(intercept=0, slope=n.pos/n.pred, color='#CC3333', size=1, alpha=.5) +
    geom_abline(intercept=0, slope=1, color='#0066CC', size=1, alpha=.5) + 
    scale_x_continuous(breaks=seq(0, max.x, by=4)) +
    coord_cartesian(xlim=c(0, n.pred), ylim=c(0, n.pos), expand = FALSE) 
}

PlotHoldoutSensitivityFinal <- function(d.ho.lift){
  max.x <- max(d.ho.lift$x)
  d.ho.lift %>%
    ggplot(aes(x=x, y=sensitivity)) + geom_line() +
    scale_x_continuous(breaks = seq(0, max.x, by=4)) +
    scale_y_continuous(labels=scales::percent) + 
    theme_bw() + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=.5)) +
    xlab('Number of Predictions Used') +
    ylab('% of Predictions Actually Sensitive') 
}

