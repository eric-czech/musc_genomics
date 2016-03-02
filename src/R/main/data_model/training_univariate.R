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

##### Data Prep #####

# Choose dataset to use for modeling (must pick one of the following)
EnableCosmic()
#EnableCtd()

# 
d.prep <- GetTrainingData(TRAIN_CACHE, RESPONSE_TYPE, RESPONSE_SELECTOR, min.mutations=3)

X <- d.prep %>% select(-tumor_id, -response)
y <- DichotomizeOutcome(d.prep[,'response'], threshold = RESPONSE_THRESH)

##### Univariate Feature Scoring #####

#' @title Compute score for single classification feature
#' @description Computes score for feature as either p-value 
#' for Fisher's exact test if the feature is binary or as
#' the p-value from a two-sided t-test otherwise
#' @param x feature to get score for 
#' @param y response factor (e.g. 'pos' or 'neg')
GetFeatureScore <- function(x, y) {
  if (length(unique(x)) == 2){
    pv <- try(fisher.test(factor(x), y)$p.value, silent = TRUE)
  } else {
    pv <- try(t.test(x ~ y)$p.value, silent = TRUE)
  }
  if (any(class(pv) == "try-error") || is.na(pv) || is.nan(pv)) pv <- NA
  pv
}

# Compute scores for all features
scores <- sapply(X, function(x){ GetFeatureScore(x, y) } )

# Rank features across the whole set as well as by type (e.g. Copy Number vs Gene Expression)
gene.scores <- data.frame(score=scores) %>% add_rownames(var='feature') %>% 
  mutate(type=str_extract(feature, '.*?(?=\\.)')) %>% 
  mutate(gene=str_replace(feature, paste0(type, '.'), '')) %>% 
  mutate(overall.rank=dense_rank(score)) %>% 
  arrange(type, overall.rank) %>%
  group_by(type) %>% mutate(type.rank=row_number()) %>%
  ungroup %>% arrange(overall.rank)

# Compute average ranks for genes with both copy number and gene expression
gene.avgs <- gene.scores %>% 
  filter(type %in% c('ge', 'cn')) %>%
  group_by(gene) %>% summarise(avg.rank=mean(type.rank), n=n()) %>% 
  ungroup %>% filter(n > 1) %>% arrange(avg.rank) %>% 
  select(-n) %>% head(100)

export.path <- '~/repos/musc_genomics/src/R/main/data_export/univariate_analysis'
gene.scores %>% write.csv(file=file.path(export.path, 'gene.scores.csv'), row.names=F)
gene.avgs %>% write.csv(file=file.path(export.path, 'gene.averages.csv'), row.names=F)

##### Single Feature Plots #####

data.frame(x=X[,'ge.KRAS'], y=y) %>% 
  ggplot(aes(x=y, y=x)) + geom_boxplot() + theme_bw()

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
univ.feat <- 'ge.KRAS'

acc <- get.univariate.acc(X[,univ.feat], y)
acc %>% ggplot(aes(x=brk, y=acc)) + geom_line() + 
  theme_bw() + ggtitle(univ.feat)

