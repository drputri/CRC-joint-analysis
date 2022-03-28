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
perf_eval <- function(k, X, y, df.feat, fold=3, group.name, alpha=1) {
  
  require(glmnet)
  require(limma)
  require(magrittr)
  require(mixOmics)
  require(caret)
  
  if(class(X) != "matrix"){
    X <- as.matrix.data.frame(X)
  }
  
  # Top feature
  df.feat <- as.data.frame(df.feat)
  
  df.feat <- df.feat[order(-df.feat$Freq),] 
  feat.to.select <- droplevels(df.feat$Feat[1:k])
  
  # sample 2/3 of the subjects at random
  id_keep <- sample(x=rownames(X),
                    size=floor(2*nrow(X)/3),
                    replace=FALSE)
  
  out <- vector("list", k-1)
  
  i=1
  
  for(f in 2:length(feat.to.select)) {
    
    f.red <- droplevels(feat.to.select[1:f])
    
    # uArray matrix
    df.train <- X[rownames(X) %in% id_keep, colnames(X) %in% f.red]
    df.test <- X[!(rownames(X) %in% id_keep), colnames(X) %in% f.red]
      
    # Response / group
    y.train <- y[names(y) %in% id_keep]
    y.test <- y[!(names(y) %in% id_keep)]
      
      
    # Ridge regression
    m0 <- cv.glmnet(type.measure="class", x=df.train, y=y.train,
                    alpha=alpha, family=binomial, nfolds=fold)
    m1 <- glmnet(x=as.matrix(df.train),
                 y=y.train,
                 family="binomial",
                 alpha=alpha,
                 lambda=m0$lambda.min)
    pred <- predict(m1, newx=as.matrix(df.test), type="class") %>%
      data.frame() %>%
      rename(pred=s0) %>%
      mutate(pred=factor(pred, levels=c("Adenoma", "CRC")))
      
    # Calculate the confusion matrix
      conf <- caret::confusionMatrix(pred$pred, as.factor(y.test))
      
      df.eval <- data.frame(BER = get.BER(t(conf$table)),
                            mce = 1-conf$overall[["Accuracy"]],
                            kappa = conf$overall[["Kappa"]],
                            specificity = conf$byClass[["Specificity"]],
                            sensitivity = conf$byClass[["Sensitivity"]],
                            ppv = conf$byClass[["Pos Pred Value"]],
                            npv = conf$byClass[["Neg Pred Value"]])
    
    # Store the information about the test set
      subjectPred.temp <- data.frame(sampleID=rownames(pred),
                                     pred, 
                                     obs=y.test
                                     )
      
      out[[i]] <- list(df.eval, subjectPred.temp)
      names(out[[i]]) <- c("evaluation", "subject prediction")
      names(out)[[i]] <- f
      
      i <- i+1
    
  }
  
  return(out)
}



