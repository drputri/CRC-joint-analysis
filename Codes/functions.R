############################################################################################
## FUNCTIONS FOR ELASTIC NET 
## FEATURE SELECTION & EVALUATION
############################################################################################
## Library
library(sampling)
library(mixOmics)
library(glmnet)
library(CMA)
library(caret)



## Stratified samples
StratSampleSelection <- function(x, strat.by = NULL, sample.id = "row.names", fold = 3) {
  require(sampling)
  if(class(x) != "data.frame") {
    x <- as.data.frame(x)
  }
  
  if(sample.id != "row.names") {
    rownames(x) <- NULL
  }
  
  group <- unique(x[[strat.by]])
  n.group <- length(unique(x[[strat.by]]))
  n <- nrow(x)
  size <- n/fold
  size.int <- floor(size)     # size for testing 
  
  # Calculate the sample proportion by group
  tabx <- table(x[[strat.by]])
  prop <- tabx/sum(tabx)
  mySample <- round(prop*(n-size.int), digits = 0)
  
  strat.sampled <- strata(x, stratanames = strat.by, size = mySample, 
                          method = "srswr")
  
  if(sample.id != "row.names") {
    id <- as.factor(x[[sample.id]][strat.sampled$ID_unit])
  } else {
    id <- as.factor(rownames(x)[strat.sampled$ID_unit])
  }
  

  return(id)
}



## Function for feature selection
select_feature <- function(X, y, fold = 3, group.name = NULL, alpha = 1, seed = 1) {
  set.seed(seed)
  
  if(class(X) != "matrix"){
    X <- as.matrix.data.frame(X)
  }
  
  if(is.null(group.name)){
    cv_enet <- cv.glmnet(type.measure = "class", x = X, y = as.vector.factor(y),
                         alpha = alpha, family='binomial', nfolds = fold)
    md <- glmnet(X, as.vector.factor(y), 
                 family = "binomial", alpha = alpha, lambda = cv_enet$lambda.1se)
  } else {
    cv_enet <- cv.glmnet(type.measure = "class", x = X, y = as.vector.factor(y[[group.name]]),
                         alpha = alpha, family='binomial', nfolds = fold)
    md <- glmnet(X, as.vector.factor(y[[group.name]]), 
                 family = "binomial", alpha = alpha, lambda = cv_enet$lambda.1se)
  }
  
  
  
  ## Get the gene's name
  co <- as.data.frame.matrix(coef(md))
  inds <-  which(co != 0)
  co <-  rownames(co)[inds]
  co <- co[which(co != "(Intercept)")]
  
  return(co)
}



## Function for feature selection using LassoCMA function
select_feat_CMA <- function(X, y, fold = 3, seed = 1) {
  set.seed(seed)
  
  if(class(X) != "matrix"){
    X <- as.matrix.data.frame(X)
  }
  
  CV <- GenerateLearningsets(y = y, method = "CV", fold = fold, strat = TRUE)
  varsel <- GeneSelection(X = X, y = y, learningsets = CV, method = "lasso", trace = FALSE)
  
  var <- NULL
  for(i in 1:fold) {
    temp <- as.data.frame(toplist(varsel, iter = i, k = 50, show = FALSE))
    temp <- temp[temp$importance !=0, ]
    
    if(is.null(var)) {
      var <- temp
    } else {
      var <- rbind(var, temp)
    }
  }
  
  return(var)
}



## Function for feature selection: keeping uArray or uBiome in the model
select_feature_penalty <- function(X, y, fold = 3, group.name, seed = 1, name.penalty) {
  set.seed(seed)
  
  penalty.fctr <- ifelse(grepl(name.penalty, colnames(df.new), ignore.case = TRUE), 0, 1)
  
  cv_enet <- cv.glmnet(type.measure = "class", x = as.matrix.data.frame(X), y = as.vector.factor(y[[group.name]]),
                       family='binomial', nfolds = fold, penalty.factor = penalty.fctr)
  md <- glmnet(as.matrix.data.frame(X), as.vector.factor(y[[group.name]]), family = "binomial", lambda = cv_enet$lambda.1se)
  
  
  ## Get the gene's name
  co <- as.data.frame.matrix(coef(md))
  inds <-  which(co != 0)
  co <-  rownames(co)[inds]
  co <- co[which(co != "(Intercept)")]
  
  return(co)
}







## Function for model evaluation: CMA
summary.stat.cma <- function(class.lasso) {
  require(mixOmics)
  require(caret)
  
  df.eval <- NULL
  for(i in 1:length(class.lasso)) {
    pred <- data.frame(y = class.lasso[[i]]@y,
                       yhat = class.lasso[[i]]@yhat)
    conf <- caret::confusionMatrix(as.factor(pred$yhat), as.factor(pred$y))
    temp <- data.frame(BER = get.BER(t(conf$table)),
                       mce = 1-conf$overall[["Accuracy"]],
                       kappa = conf$overall[["Kappa"]],
                       specificity = conf$byClass[["Specificity"]],
                       sensitivity = conf$byClass[["Sensitivity"]],
                       ppv = conf$byClass[["Pos Pred Value"]],
                       npv = conf$byClass[["Neg Pred Value"]])
    if(is.null(df.eval)) {
      df.eval <- temp
    } else {
      df.eval <- rbind(df.eval, temp)
    }
  }
  df.eval <- colMeans(df.eval)
}

enet_eval_CMA <- function(k, X, y, df.feat, fold = 3, iti = 1000, df.feat.all = TRUE, feat.name) {
  require(magrittr)
  require(mixOmics)
  require(caret)
  
  out <- NULL
  
  # Look for the top K genera
  df.feat <- as.data.frame(df.feat)
  if(df.feat.all) {
    df.feat <- df.feat[df.feat$ind == 1, ]
  }
  df.feat <- df.feat[order(-df.feat$freq),] 
  
  feat.to.select <- df.feat[[feat.name]][1:k]
  
  X <- X[, colnames(X) %in% feat.to.select]
  
  for(i.k in 1:iti) {
    set.seed(i.k)
    
    CV <- GenerateLearningsets(y = y, method ="CV",
                               fold = fold, strat = TRUE)
    tune.lasso <- CMA::tune(X = as.matrix.data.frame(X), y = group, learningsets = CV, 
                            classifier = LassoCMA, trace = FALSE)
    cl.lasso <- classification(X = as.matrix.data.frame(X), y = group, learningsets = CV, 
                               tuneres = tune.lasso, classifier = LassoCMA, trace = FALSE)
    temp <- summary.stat.cma(cl.lasso) 
    
    if(is.null(out)) {
      out <- temp
    } else {
      out <- rbind(out, temp)
    }
  }
  
  rownames(out) <- NULL
  
  return(out)
}



## Function for model evaluation
enet_eval <- function(k, X, y, df.feat, fold = 3, group.name, sample.id = "row.names",
                      iti = 1000, df.feat.all = TRUE, feat.name, alpha = 1) {
  require(glmnet)
  require(limma)
  require(magrittr)
  require(mixOmics)
  require(caret)
  
  
  ## Split data set between train & test, stratified by sample proportion
  test <- as.data.frame(StratSampleSelection(x = y, strat.by = group.name, sample.id = "row.names", fold = fold))
  colnames(test) <- c("sampleID")
  
  # uArray matrix
  df.train <- merge(test, X, by.x = "sampleID", by.y = "row.names", all.x = TRUE)
  df.train <- subset(df.train, select = -c(sampleID))
  df.test <- X[!(rownames(X) %in% as.character(test$sampleID)), ]
  
  # Response / group
  if(sample.id == "row.names") {
    y.train <- merge(test, y, by.x = "sampleID", by.y = "row.names", all.x = TRUE)
    y.train <- as.vector.factor(y.train[[group.name]])
    y.test <- as.vector(y[[group.name]][!(rownames(y) %in% as.character(test$sampleID))])
  } else {
    y.train <- merge(test, y, by = "sampleID", all.x = TRUE)
    y.train <- as.vector.factor(y.train[[group.name]])
    y.test <- as.vector.factor(y[[group.name]][!(y$sampleID %in% as.character(test$sampleID))])
  }
  
  
  # Top feature
  df.feat <- as.data.frame(df.feat)
  if(df.feat.all) {
    df.feat <- df.feat[df.feat$ind == 1, ]
  }
  df.feat <- df.feat[order(-df.feat$freq),] 
  
  out <- list()
  
  feat.to.select <- df.feat[[feat.name]][1:k]
  
  # X2 <- as.matrix.data.frame(X[, colnames(X) %in% feat.to.select])
  # y2 <- as.vector.factor(y[[group.name]])
  # y2 <- ifelse(y2 == "Adenoma", 0, 1)
  
  df.train2 <- as.matrix.data.frame(df.train[, colnames(df.train) %in% feat.to.select])
  df.test2 <- as.matrix.data.frame(df.test[, colnames(df.test) %in% feat.to.select])
  
  out <- vector("list", iti)
  for(i.k in 1:iti) {
    set.seed(i.k)
    cv.enet <- cv.glmnet(type.measure = "class", x = df.train2, y = y.train,
                         alpha = alpha, family='binomial', nfolds = fold, parallel = TRUE)
    pred_enet = as.data.frame(predict(cv.enet, newx = df.test2, s = "lambda.1se", type = "response",
                        alpha = alpha))
    colnames(pred_enet) <- "pred"
    pred_enet$pred <- as.factor(ifelse(pred_enet$pred > 0.5, "CRC", "Adenoma"))
    conf <- caret::confusionMatrix(pred_enet$pred, as.factor(y.test))
    
    df.eval <- data.frame(BER = get.BER(t(conf$table)),
                          mce = 1-conf$overall[["Accuracy"]],
                          kappa = conf$overall[["Kappa"]],
                          specificity = conf$byClass[["Specificity"]],
                          sensitivity = conf$byClass[["Sensitivity"]],
                          ppv = conf$byClass[["Pos Pred Value"]],
                          npv = conf$byClass[["Neg Pred Value"]])
    subjectPred.temp <- data.frame(pred_enet, 
                                   obs = y.test,
                                   sampleID = rownames(pred_enet))
    
    out[[i.k]] <- list(df.eval, subjectPred.temp)
    names(out[[i.k]]) <- c("evaluation", "subject prediction")
    names(out)[[i.k]] <- i.k
  }
  return(out)
}



