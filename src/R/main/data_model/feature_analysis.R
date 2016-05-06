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

d.prep[,c(top.feats$feature, 'tumor_id')] %>%
  filter(origin == 'LUNG' & tissue == 'LUNG_SMALL_CELL_CARCINOMA') %>%
  select(-origin) %>%
  melt(id.vars=c('tumor_id')) %>%
  ggplot(aes(x=variable, y=tumor_id, fill=value)) + geom_tile()


get.feat.map <- function(d){
  d[,c(top.feats$feature, 'tumor_id', 'origin', 'tissue', 'prob')] %>%
    mutate(pred_type=ifelse(prob <= .8, ifelse(prob < 0, 'TRAIN', 'LO'), 'HI')) %>% select(-prob) %>%
    filter(origin == 'LUNG') %>%
    filter(tissue == 'LUNG_SMALL_CELL_CARCINOMA') %>% 
    select(-origin, -tissue) 
}

reorder.pca <- function(d){
  d.order <- d %>% select(-tumor_id, -pred_type)
  d[order(predict(prcomp(d.order), d.order)[,'PC1']),]
}

feat.map <- rbind(get.feat.map(d.predict), get.feat.map(d.prep %>% mutate(prob=-1)))


d.hi <- feat.map %>% filter(pred_type=='HI') %>% reorder.pca
p.hi <- plot_ly(
  x=d.hi %>% select(-pred_type, -tumor_id) %>% names,
  y=d.hi$tumor_id, 
  z=d.hi %>% select(-pred_type, -tumor_id) %>% as.matrix,
  type='heatmap', showscale = T, zmin=-3, zmax=6, 
  colorbar=list(title='Feature Value'), width=1200,
  colors='YlOrRd'
)

d.lo <- feat.map %>% filter(pred_type=='LO') %>% reorder.pca
p.lo <- plot_ly(
  x=d.lo %>% select(-pred_type, -tumor_id) %>% names,
  y=d.lo$tumor_id, 
  z=d.lo %>% select(-pred_type, -tumor_id) %>% as.matrix, 
  type='heatmap', showscale = F, zmin=-3, zmax=6
)

d.tr <- feat.map %>% filter(pred_type=='TRAIN') %>% reorder.pca
d.tr.res <- d.prep[,c('tumor_id', 'response')] %>%
  mutate(response = as.character(DichotomizeOutcome(response, RESPONSE_THRESH)))
d.tr <- d.tr %>% left_join(d.tr.res, by='tumor_id') %>%
  mutate(tumor_id=sprintf('%s (%s)', tumor_id, response)) %>% select(-response)

p.tr <- plot_ly(
  x=d.tr %>% select(-pred_type, -tumor_id) %>% names,
  y=d.tr$tumor_id, 
  z=d.tr %>% select(-pred_type, -tumor_id) %>% as.matrix, 
  type='heatmap', showscale = F, zmin=-3, zmax=6
)


subplot(p.hi, p.lo, p.tr, margin = 0.05, nrows=3) %>% layout(
  showlegend = F, 
  title='Feature Values for Training and Predicted Cell Lines (by Hi/Lo prediction)', 
  xaxis=list(title='Cell Lines w/ High Predicted Sensitivity'),
  yaxis=list(title=''),
  xaxis2=list(title='Cell Lines w/ Low Predicted Sensitivity'),
  yaxis2=list(title=''),
  xaxis3=list(title='Training Cell Lines'),
  yaxis3=list(title='')
)



# Defunct Now:

# 
# # feat.map %>%
# #   ggplot(aes(x=variable, y=tumor_id, fill=value)) + geom_tile() +
# #   facet_wrap(~pred_type, scales='free_y', ncol=1)
#     
# 
# prep.features <- function(d, pred.type){
#   r <- d %>% filter(pred_type == pred.type)
#   r <- reshape2::dcast(r, tumor_id ~ variable, value.var='value') 
#   rownames(r) <- as.character(r$tumor_id)
#   r <- r %>% select(-tumor_id)
#   r[order(predict(prcomp(r), r)[,'PC1']),]
# }
# 
# ### Single Prob Type (training data)
# d.hi <- prep.features(feat.map, 'HI')
# plot_ly(
#   z = unname(d.hi), x = colnames(d.hi), y = rownames(d.hi), 
#   type = "heatmap", #colors='YlOrRd', 
#   showscale = F, zmin=-3, zmax=6, width=1200
# )
# 
# ### Binary Prob Type (prediction data)
# d.hi <- prep.features(feat.map, 'HI')
# d.lo <- prep.features(feat.map, 'LO')
# high <- as.matrix(d.hi)
# low <- as.matrix(d.lo)
# subplot(
#   plot_ly(z = high, x = colnames(d.hi), y = rownames(d.hi), type = "heatmap", colors='YlOrRd', showscale = F, zmin=-3, zmax=6, width=1200),
#   plot_ly(z = low, x = colnames(d.lo), y = rownames(d.lo), type = "heatmap", colors='YlOrRd', showscale = T, zmin=-3, zmax=6, colorbar=list(title='Feature Value'), width=1200),
#   margin = 0.05, 
#   nrows=2
# ) %>% layout(
#   showlegend = F, 
#   title='Feature Values for Hi/Lo Predictions', 
#   xaxis=list(title='Cell Lines w/ High Predicted Sensitivity'),
#   yaxis=list(title=''),
#   xaxis2=list(title='Cell Lines w/ Low Predicted Sensitivity'),
#   yaxis2=list(title='')
# )
# 
# 
