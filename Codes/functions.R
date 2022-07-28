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
perf_eval <- function(k, X, y, rank.by, order.desc = TRUE, 
                      df.feat, fold = 3, group.name, alpha = 1,
                      MCCV.type = 0) {
  
  require(glmnet)
  require(limma)
  require(magrittr)
  require(mixOmics)
  require(caret)
  
  # Combine the RNA and group
  dt <- merge(y, X, by="row.names")
  rownames(dt) <- dt$Row.names 
  dt <- subset(dt, select=-c(Row.names))
  
  # Top feature
  df.feat <- as.data.frame(df.feat)
  
  if(order.desc) {
    df.feat <- df.feat[order(-df.feat[[rank.by]]),] 
  } else {
    df.feat <- df.feat[order(df.feat[[rank.by]]),] 
  }
  
  feat.to.select <- droplevels(df.feat$Feat[1:k])
  
  
  # Samples on the different MCCV type
  if(MCCV.type==0) {
    
    # sample 2/3 of the subjects at random
    id_keep <- sample(x=rownames(X),
                      size=floor(2*nrow(X)/3),
                      replace=FALSE)
    
  } else if(MCCV.type==1) {
    
    # MCCV type 1: 11-11 samples on test set
    ad <- names(group)[group=="Adenoma"]
    crc <- names(group)[group=="CRC"]
    id_keep <- c(sample(x=ad, size=30, replace=FALSE),
                 sample(x=crc, size=10, replace=FALSE))
    
  } else if (MCCV.type==2) {
    
    # MCCV type 2: 20-10 samples on test set
    ad <- names(group)[group=="Adenoma"]
    crc <- names(group)[group=="CRC"]
    id_keep <- c(sample(x=ad, size=21, replace=FALSE),
                 sample(x=crc, size=11, replace=FALSE))
    
  } else {
    
    # MCCV type 2: 15-6 samples on test set
    ad <- names(group)[group=="Adenoma"]
    crc <- names(group)[group=="CRC"]
    id_keep <- c(sample(x=ad, size=26, replace=FALSE),
                 sample(x=crc, size=15, replace=FALSE))
    
  }
  
  
  
  out <- vector("list", k-1)
  
  i=1
  
  for(f in 2:length(feat.to.select)) {
    
    f.red <- droplevels(feat.to.select[1:f])
    
    # uArray matrix
    df.train <- dt[rownames(dt) %in% id_keep, colnames(dt) %in% c("x", as.character(f.red))]
    df.test <- dt[!rownames(dt) %in% id_keep, colnames(dt) %in% c("x", as.character(f.red))]
      
    # # Response / group
    # y.train <- y[names(y) %in% id_keep]
    # y.test <- y[!(names(y) %in% id_keep)]
      
      
    # # Ridge regression
    # m0 <- cv.glmnet(type.measure="class", x=df.train, y=y.train,
    #                 alpha=alpha, family=binomial, nfolds=fold)
    # m1 <- glmnet(x=as.matrix(df.train),
    #              y=y.train,
    #              family="binomial",
    #              alpha=alpha,
    #              lambda=m0$lambda.min)
    
    
    # Logistic regression
    m0 <- glm(x ~ ., data=df.train, family="binomial")
    probs <- predict(m0, newdata=df.test[, -1], type="response")
    pred <- ifelse(probs > 0.5, "CRC", "Adenoma")
    
    # pred <- predict(m1, newx=as.matrix(df.test), type="class") %>%
    #   data.frame() %>%
    #   rename(pred=s0) %>%
    #   mutate(pred=factor(pred, levels=c("Adenoma", "CRC")))
      
    # Calculate the confusion matrix
      conf <- caret::confusionMatrix(as.factor(pred), df.test$x)
      
      df.eval <- data.frame(BER = conf$byClass[["Balanced Accuracy"]],
                            mce = 1-conf$overall[["Accuracy"]],
                            kappa = conf$overall[["Kappa"]],
                            specificity = conf$byClass[["Specificity"]],
                            sensitivity = conf$byClass[["Sensitivity"]],
                            ppv = conf$byClass[["Pos Pred Value"]],
                            npv = conf$byClass[["Neg Pred Value"]])
    
    # Store the information about the test set
      subjectPred.temp <- data.frame(sampleID=names(pred),
                                     pred, 
                                     obs=df.test$x
                                     )
      
      out[[i]] <- list(df.eval, subjectPred.temp)
      names(out[[i]]) <- c("evaluation", "subject prediction")
      names(out)[[i]] <- f
      
      i <- i+1
    
  }
  
  return(out)
}



