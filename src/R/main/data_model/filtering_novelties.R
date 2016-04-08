
#' @title Collect information about the origin of tumors in never-before-seen samples
NoveltyAttributesByOrigin <- function(d.predict, d.tr.uk, d.ho.uk, origin.trans){
  
  # Get a unique list of tumor origins seen in training
  uk.tr.origin <- unique(d.tr.uk$X[,'origin'])
  uk.tr.origin.common <- origin.trans$convert(uk.tr.origin) %>% .[. > 0] %>% names
  
  # For each origin seen in new samples, return flags indicating whether or 
  # not that origin is a "common" origin (occurs frequently) and whether or
  # not it was seen in training
  sapply(d.ho.uk$X[,'origin'], function(x) c(
    is.common=(x %in% uk.tr.origin.common),
    is.known=(x %in% uk.tr.origin)
  )) %>% t %>% as.data.frame %>% add_rownames('origin') %>%
    mutate(tumor_id = d.predict[,'tumor_id']) # Attach tumor id's for convenience
}

#' @title Build novelty detection model using a two-class SVM
NoveltyAttributesTwoClassSVM <- function(X.tr.nov, X.ho.nov, seed, N.sim.factor=3, n.core=3){
  
  # Create simulated data meant to represent data similar to that which would be expected
  # in never-before-seen samples that deviate significantly from what was seen in training.
  N.sim <- nrow(X.tr.nov)*N.sim.factor
  X.tr.sim <- apply(X.tr.nov, 2, function(x) runif(n=N.sim, min=min(x), max=max(x))) %>% 
    data.frame %>% setNames(names(X.tr.nov))
  
  # Combine the actual training data with the simulated data and give the records from
  # each a label that will be used as the supervised learning target
  d.sim <- rbind(X.tr.nov %>% mutate(InSample='Yes'), X.tr.sim %>% mutate(InSample='No')) %>% 
    mutate(InSample=factor(InSample, levels=c('Yes', 'No')))
  rownames(d.sim) <- NULL
  
  # Train a binary classification SVM model on the above data which can then be used to
  # detect anomalous samples in unlabled datasets
  registerDoMC(n.core)
  set.seed(seed)
  m.svm <- train(
    InSample ~ ., data=d.sim, method='svmRadial', tuneLength=10,
    trControl=trainControl(method='cv', number=10, classProbs=T, verboseIter=T)
  )
  
  # Make predictions on which samples seen in unlabeled dataset are anomalous,
  # with respect to the original training data
  p <- predict(m.svm, X.ho.nov, type='prob')[,1] # Produces probability predictions
  
  # Return both the model and the predictions
  list(novelty.model=m.svm, novelty.predictions=p)
}

#' @title Build novelty detection model using a one-class SVM
NoveltyAttributesOneClassSVM <- function(X.tr.nov, X.ho.nov, seed){
  require(e1071)
  set.seed(seed)
  
  # Train one class SVM on training data alone
  m.svm <- svm(X.tr.nov %>% as.matrix, y=NULL, type='one-classification', nu=0.10, scale=FALSE, kernel="radial")
  
  # Make anomaly predictions for unlabeled holdout data
  p <- predict(m.svm, X.ho.nov %>% as.matrix) # Produces boolean predictions
  
  # Return both the model and the predictions
  list(novelty.model=m.svm, novelty.predictions=p)
}

NoveltyAttributesByDist <- function(X.tr.nov, X.ho.nov){
  add.names <- function(d, p) {
    r <- as.matrix(d)
    rownames(r) <- sprintf('%s.%s', p, 1:nrow(d))
    r
  }
  
  # Compute pairwise distances between all training and unknown samples
  nov.dist <- dist(rbind(add.names(X.tr.nov, 'tr'), add.names(X.ho.nov, 'ho'))) %>% as.matrix %>% as.data.frame
  
  # Extract row/col names corresponding to each data type
  n.tr.dist <- nov.dist %>% select(starts_with('tr.')) %>% names
  n.ho.dist <- nov.dist %>% select(starts_with('ho.')) %>% names
  
  # Compute average distance between a single unknown holdout sample and all training samples,
  # for each holdout sample (returning a named vector of mean distances)
  ho.dist <- sapply(n.ho.dist, function(x) mean(nov.dist[n.tr.dist, x]) ) %>% setNames(n.ho.dist)
  
  ho.dist
}


PlotUnlabeledHoldout <- function(X.ho.nov, p.svm.two.class){
  require(GGally)
  ggparcoord(
    X.ho.nov %>% 
      mutate(
        p=cut(p.svm.two.class, breaks=3, labels=c('Least Likely', 'Somewhat Likely', 'Most Likely')), 
        a=1-p.svm.two.class
      ),
    columns=1:ncol(X.ho.nov),
    groupColumn=ncol(X.ho.nov)+1,
    alphaLines='a'
  ) + scale_color_discrete(guide=guide_legend(title='Pr(Sample similar to Training)')) +
    scale_alpha_continuous(guide=F) +
    theme_bw() + theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
    xlab('Feature Name') + ylab('Feature Value') 
}