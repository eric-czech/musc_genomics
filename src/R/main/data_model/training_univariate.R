#'-----------------------------------------------------------------------------
#' Univariate Feature Analysis
#'
#' This module contains code for scoring gene-based features based on univariate
#' tests (like t-tests) as well as determining the accuracy of classification
#' rules based on single features.
#' 
#' @author eczech
#'-----------------------------------------------------------------------------
source('utils.R')
source('data_model/training_filter.R')
source('data_model/filtering_lib.R')

##### Data Prep #####

# Choose dataset to use for modeling (must pick one of the following)
EnableCosmic()
#EnableCtd()

# 
d.prep <- GetTrainingData(TRAIN_CACHE, RESPONSE_TYPE, RESPONSE_SELECTOR, min.mutations=3)

X <- d.prep %>% select(-tumor_id, -response) %>% mutate(origin=TransformOriginSolidLiquid(origin))
y <- DichotomizeOutcome(d.prep[,'response'], threshold = RESPONSE_THRESH)

##### Univariate Feature Scoring #####

# Compute scores for all features
options(stringsAsFactors=F)
registerDoMC(3)
scores <- foreach(xn=names(X), .combine=rbind)%dopar%{
  x <- X[,xn]
  score <- FeatureScore(x, y, t.test.alt='two.sided')
  if (length(unique(x)) == 2) {
    direction <- 'na'
  } else {
    mu <- mean(x[y == 'pos']) - mean(x[y == 'neg'])
    direction <- ifelse(mu > 0, 'greater', 'less')
  }
  data.frame(feature=xn, score=score, direction=direction)
}


# Rank features across the whole set as well as by type (e.g. Copy Number vs Gene Expression)
gene.scores <- scores %>%
  mutate(type=str_extract(feature, '.*?(?=\\.)')) %>% 
  mutate(gene=str_replace(feature, paste0(type, '.'), '')) %>% 
  mutate(overall.rank=dense_rank(score)) %>% 
  arrange(type, overall.rank) %>% group_by(type) %>% 
  mutate(type.rank=row_number()) %>% ungroup %>% 
  arrange(direction, overall.rank) %>% group_by(direction) %>% 
  mutate(direction.rank=row_number()) %>% ungroup %>% 
  arrange(overall.rank)

# Compute average ranks for genes with both copy number and gene expression
gene.avgs <- gene.scores %>% 
  filter(type %in% c('ge', 'cn')) %>% group_by(gene) %>% 
  dplyr::summarise(
    avg.rank=mean(type.rank), 
    direction=ifelse(length(unique(direction))==1, direction[1], 'opposing'), 
    n=n()) %>% 
  ungroup %>% filter(n > 1) %>% arrange(avg.rank) %>% 
  select(-n)

export.path <- '~/repos/musc_genomics/src/R/main/data_export/univariate_analysis'
gene.scores %>% write.csv(file=file.path(export.path, 'gene.scores.csv'), row.names=F)
gene.avgs %>% write.csv(file=file.path(export.path, 'gene.averages.csv'), row.names=F)



##### Feature Plots #####

data.frame(x=X[,'cn.SMAD4'], y=y) %>% 
  ggplot(aes(x=y, y=x)) + geom_boxplot() + theme_bw()


top.avg.feats <- gene.avgs %>% filter(direction=='greater') %>% head(9) %>% .$gene
# top.avg.feats <- gene.scores %>% head(9) %>% .$gene
# gene.avgs %>% filter(gene %in% top.avg.feats)

top.avg.feats <- foreach(gene=top.avg.feats, .combine=rbind) %do%{
  rbind(
    data.frame(gene=gene, value=X[,paste0('cn.', gene)], type='CopyNum', y=y),
    data.frame(gene=gene, value=X[,paste0('ge.', gene)], type='Expression', y=y)
  )
}
top.avg.feats %>%
  ggplot(aes(x=type, y=value, color=y)) + geom_boxplot(position='dodge') +
  facet_wrap(~gene, scales='free') + 
  #theme_bw() + ggtitle('Top Features Across CN and GE') +
  theme_bw() + ggtitle('Top Standalone Features') +
  scale_color_discrete(guide=guide_legend(title='Sensitivity')) +
  ggsave(file.path(export.path, 'top.genes.plot.png'))

##### Univariate Classification Rules #####

#' @title Compute pos + neg class accuracy using a single feature
#' @return data frame with accuracy for each break in feature
#' tested as classification rule
get.univariate.acc <- function(x, y){
  d <- data.frame(x=x, y=y)
  breaks <- seq(range(x)[1], range(x)[2], len=100)
  breaks <- breaks[-length(breaks)]
  breaks <- breaks[-1]
  foreach(i=breaks, .combine=rbind)%do%{
    y.lo <- subset(d, x <= i)$y
    y.hi <- subset(d, x > i)$y
    acc <- (sum(y.lo == 'neg') + sum(y.hi == 'pos')) / length(y)
    data.frame(brk=i, acc=acc, x.max=max(x), x.min=min(x))
  }
}

# Plot accuracy profile for univariate rule
#univ.feat <- 'ge.YAP1'
univ.feat <- 'ge.SMAD4'

acc <- get.univariate.acc(X[,univ.feat], y)
acc %>% ggplot(aes(x=brk, y=acc)) + geom_line() + 
  theme_bw() + ggtitle(univ.feat)

