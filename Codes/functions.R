############################################################################################
## FUNCTIONS FEATURE SELECTION & EVALUATION                                                #
############################################################################################
## Library
library(sampling)
library(mixOmics)
library(glmnet)
library(CMA)
library(caret)



## Function for feature selection
select_feature <- function(X, y, fold=3, group.name=NULL, alpha=1, seed=1) {
  set.seed(seed)
  
  if(class(X) != "matrix"){
    X <- as.matrix.data.frame(X)
  }
  
  if(is.null(group.name)){
    cv <- cv.glmnet(type.measure="class", x=X, y=as.vector.factor(y),
                    alpha=alpha, family='binomial', nfolds=fold)
    md <- glmnet(X, as.vector.factor(y), 
                 family="binomial", alpha=alpha, lambda=cv$lambda.min)
  } else {
    cv <- cv.glmnet(type.measure="class", x=X, y=as.vector.factor(y[[group.name]]),
                    alpha = alpha, family='binomial', nfolds = fold)
    md <- glmnet(X, as.vector.factor(y[[group.name]]), 
                 family="binomial", alpha=alpha, lambda=cv$lambda.min)
  }
  
  
  
  ## Get the gene's name
  co <- as.data.frame.matrix(coef(md))
  inds <-  which(co != 0)
  co <-  rownames(co)[inds]
  co <- co[which(co != "(Intercept)")]
  
  return(co)
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







## Function for model evaluation
perf_eval <- function(k, X, y, df.feat, fold=3, group.name,
                      iti=1000, df.feat.all=TRUE, feat.name, alpha=1) {
  require(glmnet)
  require(limma)
  require(magrittr)
  require(mixOmics)
  require(caret)
  
  if(class(X) != "matrix"){
    X <- as.matrix.data.frame(X)
  }
  
  # sample 2/3 of the subjects at random
  id_keep <- sample(x=rownames(X),
                    size=floor(2*nrow(X)/3),
                    replace=FALSE)
  
  # uArray matrix
  df.train <- X[rownames(X) %in% id_keep, ]
  df.test <- X[!(rownames(X) %in% id_keep), ]
  
  # Response / group
  y.train <- y[names(y) %in% id_keep]
  y.test <- y[!(names(y) %in% id_keep)]
  
  # Top feature
  df.feat <- as.data.frame(df.feat)
  
  if(df.feat.all) {
    df.feat <- df.feat[df.feat$Ind==1, ]
  }
  
  df.feat <- df.feat[order(-df.feat$Freq),] 
  
  out <- list()
  
  feat.to.select <- df.feat$Feat
  
  out <- vector("list", iti)
  
  for(f in 2:length(feat.to.select)) {
    
    f.red <- droplevels(feat.to.select[1:f])
    x.train <- as.matrix.data.frame(df.train[, colnames(df.train) %in% f.red])
    x.test <- as.matrix.data.frame(df.test[, colnames(df.test) %in% f.red])
    
    for (i in 1:iti) {
      
      set.seed(i)
      
      m0 <- cv.glmnet(type.measure="class", x=x.train, y=y.train,
                      alpha=alpha, family='binomial', nfolds=fold, parallel=TRUE)
      pred <- as.data.frame(predict(m0, newx=x.test, s="lambda.min", type = "response",
                                        alpha = alpha))
        
    }
    
  }
  
  
  
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



