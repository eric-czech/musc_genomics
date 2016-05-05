#'-----------------------------------------------------------------------------
#' Feature Analysis Script
#'
#' This script conatins logic to analyze selected feature values for various
#' cancer types (e.g. small cell lung cancer)
#' 
#' @author eczech
#'-----------------------------------------------------------------------------

source('utils.R')
source('data_model/training_lib.R')
source('data_model/training_models.R')
source('data_model/training_filter.R')
source('data_model/filtering_lib.R')
source('~/repos/portfolio/functional/ml/R/trainer.R')
source('data_model/filtering_models.R')
source('data_model/filtering_settings.R')
library(plotly)

## Choose dataset to use for modeling (must pick one of the following)
EnableCosmic()
#EnableCtd()

RES_CACHE_DIRNAME <- 'filter_result_data.ga'
RESULT_CACHE <- Cache(dir=file.path(CACHE_DIR, RES_CACHE_DIRNAME), project=RESPONSE_TYPE)

##### Data Loading #####

# Load training data for this response type
d.prep <- GetTrainingData(TRAIN_CACHE, RESPONSE_TYPE, RESPONSE_SELECTOR, min.mutations=3)

# Load a dataset equivalent to the above but with missing response labels (ultimately to be predicted)
d.unk <- GetPredictionData(TRAIN_CACHE, RESPONSE_TYPE, PREDICTION_SELECTOR, names(d.prep)) %>% select(-tissue)

# Add tissue type for truly unknown cell lines
subtype.path <- file.path(CACHE_DIR, 'raw_data/lung_cancer_cell_line_subtypes.csv')
d.subtype <- read.csv(subtype.path, stringsAsFactors=F) %>%
  select(tumor_id, tissue=subtype, prob=ensglm.prob) %>% 
  mutate(tissue=toupper(tissue))
d.predict <- d.unk %>% left_join(d.subtype, by='tumor_id') %>%
  select(tumor_id, origin, tissue, prob, response, everything())

# Make sure every single record in subtype file joined successfully 
# (it was a manually generated dataset so it should be perfectly aligned with known tumor ids)
if (sum(!is.na(d.predict$tissue)) != nrow(d.subtype))
  stop('Some rows in subtype data frame failed to join successfully')


# Fetch top features found in training
top.feats <- GetTopFeatures(RESULT_CACHE, GetOptimalFeatCt())
top.feats <- top.feats %>% filter(selected.in.uk )

# d.prep %>%
#   filter(origin == 'LUNG' & tissue == 'LUNG_SMALL_CELL_CARCINOMA') %>% 
#   select(one_of(c(top.feats$feature, 'tumor_id')), -origin) %>% 
#   melt(id.vars=c('tumor_id')) %>% 
#   ggplot(aes(x=variable, y=tumor_id, fill=value)) + geom_tile()

feat.map <- d.predict %>%
  select(one_of(top.feats$feature), tumor_id, origin, tissue, prob) %>% 
  mutate(pred_type=ifelse(prob > .8, 'HI', 'LO')) %>% select(-prob) %>%
  filter(origin == 'LUNG') %>%
  filter(tissue == 'LUNG_SMALL_CELL_CARCINOMA') %>% 
  select(-origin, -tissue) %>%
  melt(id.vars=c('tumor_id', 'pred_type')) 

# feat.map %>%
#   ggplot(aes(x=variable, y=tumor_id, fill=value)) + geom_tile() +
#   facet_wrap(~pred_type, scales='free_y', ncol=1)
    

prep.features <- function(d, pred.type){
  r <- d %>% filter(pred_type == pred.type)
  r <- reshape2::dcast(r, tumor_id ~ variable, value.var='value') 
  rownames(r) <- as.character(r$tumor_id)
  r <- r %>% select(-tumor_id)
  r[order(predict(pca, r)[,'PC1']),]
}


d.hi <- prep.features(feat.map, 'HI')
d.lo <- prep.features(feat.map, 'LO')
high <- as.matrix(d.hi)
low <- as.matrix(d.lo)
subplot(
  plot_ly(z = high, x = colnames(d.hi), y = rownames(d.hi), type = "heatmap", colors='YlOrRd', showscale = F, zmin=-3, zmax=6, width=1200),
  plot_ly(z = low, x = colnames(d.lo), y = rownames(d.lo), type = "heatmap", colors='YlOrRd', showscale = T, zmin=-3, zmax=6, colorbar=list(title='Feature Value'), width=1200),
  margin = 0.05, 
  nrows=2
) %>% layout(
  showlegend = F, 
  title='Feature Values for Hi/Lo Predictions', 
  xaxis=list(title='Cell Lines w/ High Predicted Sensitivity'),
  yaxis=list(title=''),
  xaxis2=list(title='Cell Lines w/ Low Predicted Sensitivity'),
  yaxis2=list(title='')
)


