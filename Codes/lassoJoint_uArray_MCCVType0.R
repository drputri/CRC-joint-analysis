#####################################################
## STATISTICAL ANALYSIS                            ##
## LASSO: joint                                    ##
## Prepared by: Dea Putri (dputri@its.jnj.com)     ##
#####################################################



## Load data
load("Data/rna_crcad.Rda")
load("Data/otu_431_clr.Rda")
load("Data/group_crcad.Rda")


# Get the selected features of the different MCCV type
feat0 <- get(load("Data/uArray_selected feat_MCCVType0.Rda"))
feat1 <- get(load("Data/uArray_selected feat_MCCVType1.Rda"))
feat2 <- get(load("Data/uArray_selected feat_MCCVType2.Rda"))
feat3 <- get(load("Data/uArray_selected feat_MCCVType3.Rda"))

genus0 <- get(load("Data/uBiome_selected feat_MCCVType0.Rda"))
genus1 <- get(load("Data/uBiome_selected feat_MCCVType1.Rda"))
genus2 <- get(load("Data/uBiome_selected feat_MCCVType2.Rda"))
genus3 <- get(load("Data/uBiome_selected feat_MCCVType3.Rda"))


# Order the selected features based on the frequency of selection
feat <- feat0[order(-feat0$Freq), ]
genus <- genus[order(-genus$Freq), ]


# Get top 6 uArray & 8 microbiome features
get.feat <- droplevels(feat$Feat[1:10])
get.genus <- droplevels(genus$Feat[1:8])



## Subset uArray and uBiome data
rna.10 <- rna.crcad[, colnames(rna.crcad) %in% as.character(get.feat)]
rna.10 <- cbind(rna.10, group)
otu.8 <- otu.clr[, colnames(otu.clr) %in% as.character(get.genus)]
otu.8 <- cbind(otu.8, group)

save(rna.10, file = "Data/df_lasso_10uarray_MCCVType0.Rda")
save(otu.8, file = "Data/df_lasso_8uBiome.Rda")



## Lasso data frame for non-penalized uArray (439 features)
df.uArray <- merge(rna.10, otu.clr, by = "row.names")
rownames(df.uArray) <- df.uArray$Row.names
df.uArray <- subset(df.uArray, select = -c(Row.names))
save(df.uArray, file = "Data/df_lasso_uArray_MCCVType0.Rda")



## Lasso data frame for non-penalized uBiome (9,704 features)
df.uBiome <- merge(rna.crcad, otu.8, by = "row.names")
rownames(df.uBiome) <- df.uBiome$Row.names
df.uBiome <- subset(df.uBiome, select = -c(Row.names))
save(df.uBiome, file = "Data/df_lasso_uBiome.Rda")



## Lasso data frame for all uBiome and uArray data (10,127 features)
df.all <- merge(rna.crcad, otu.clr, by = "row.names")
rownames(df.all) <- df.all$Row.names
df.all <- subset(df.all, select = -c(Row.names))
df.all <- cbind(df.all, group)
save(df.all, file = "Data/df_lasso_all.Rda")







#----------------------------------------------------------------------------------------
# PULLING 10 MICROARRAY + 431 MICROBIOME: no penalisation on microarray 
#----------------------------------------------------------------------------------------
load("Data/df_lasso_uArray_MCCVType0.Rda")

library(parallel)
library(ggplot2)


# Choose the MCCV type
MCCV.type <- 0


# Create a cluster of cores
cl <- makeCluster(getOption("cl.cores", 4))
clusterExport(cl=cl, 
              varlist=c("df.uArray", "MCCV.type"))


# Main algorithm
system.time({
  out <- parLapply(cl=cl,
                   X=1:500,
                   fun=function(i) {
                     require(glmnet)
                     require(magrittr)
                     require(mixOmics)
                     require(dplyr)
                     
                     # Set seed for the sample randomization
                     set.seed(i)
                     
                     
                     # Samples on the different MCCV type
                     if(MCCV.type==0) {
                       
                       # sample 2/3 of the subjects at random
                       id_keep <- sample(x=rownames(df.uArray),
                                         size=floor(2*nrow(df.uArray)/3),
                                         replace=FALSE)
                       
                     } else if(MCCV.type==1) {
                       
                       # MCCV type 1: 11-11 samples on test set
                       ad <- rownames(df.uArray)[df.uArray$group == "Adenoma"]
                       crc <- rownames(df.uArray)[df.uArray$group == "CRC"]
                       id_keep <- c(sample(x=ad, size=30, replace=FALSE),
                                    sample(x=crc, size=10, replace=FALSE))
                       
                     } else if (MCCV.type==2) {
                       
                       # MCCV type 2: 20-10 samples on test set
                       ad <- rownames(df.uArray)[df.uArray$group == "Adenoma"]
                       crc <- rownames(df.uArray)[df.uArray$group == "CRC"]
                       id_keep <- c(sample(x=ad, size=21, replace=FALSE),
                                    sample(x=crc, size=11, replace=FALSE))
                       
                     } else {
                       
                       # MCCV type 2: 15-6 samples on test set
                       ad <- rownames(df.uArray)[df.uArray$group == "Adenoma"]
                       crc <- rownames(df.uArray)[df.uArray$group == "CRC"]
                       id_keep <- c(sample(x=ad, size=11, replace=FALSE),
                                    sample(x=crc, size=11, replace=FALSE))
                       
                     }
                     
                     # Sample test and train data
                     df.train <- df.uArray[rownames(df.uArray) %in% id_keep, ]
                     df.test <- df.uArray[!(rownames(df.uArray) %in% id_keep), ]
                     
                     x <- model.matrix(group~., df.train)[,-1]
                     y <- as.factor(factor(df.train$group, levels = c("Adenoma", "CRC")))
                     
                     
                     #Unpenalized microarray data
                     penalty.fctr <- ifelse(grepl("SVs", colnames(x), ignore.case = TRUE), 1, 0)
                     
                     
                     # LASSO
                     cv.lasso <- NULL
                     while(is(cv.lasso, 'try-error') || is.null(cv.lasso)) {
                       cv.lasso <- try(cv.glmnet(x, y, alpha=1, family="binomial", type.measure="class",
                                                 penalty.factor=penalty.fctr, standardize=TRUE,
                                                 nfolds=5, parallel=TRUE), silent=FALSE)
                     }
                     
                     m1 <- glmnet(x, y, alpha=1, family="binomial",
                                  penalty.factor=penalty.fctr, standardize=TRUE,
                                  lambda=cv.lasso$lambda.min)
                     
                     x.test <- model.matrix(group ~., df.test)[,-1]
                     pred.out <- as.data.frame(predict(m1, newx=x.test, s="lambda.min", 
                                                       type="class", alpha=1)) %>%
                       rename(pred="s1") %>%
                       mutate(pred=factor(pred, levels=c("Adenoma", "CRC")))
                     
                     # Model performance
                     conf <- caret::confusionMatrix(pred.out$pred, as.factor(df.test$group))
                     
                     df.eval <- data.frame(BER = get.BER(t(conf$table)),
                                           mce = 1-conf$overall[["Accuracy"]],
                                           kappa = conf$overall[["Kappa"]],
                                           specificity = conf$byClass[["Specificity"]],
                                           sensitivity = conf$byClass[["Sensitivity"]],
                                           ppv = conf$byClass[["Pos Pred Value"]],
                                           npv = conf$byClass[["Neg Pred Value"]])
                     
                     # Test subject 
                     subjectPred.temp <- data.frame(pred.out, 
                                                    obs = df.test$group,
                                                    sampleID = rownames(pred.out))
                     
                     coeff <- as.data.frame.matrix(coef(m1))
                     coeff$feat <- rownames(coeff)
                     coeff <- coeff[coeff$s0 != 0, ]
                     
                     out <- list(df.eval, subjectPred.temp, coeff)
                     names(out) <- c("evaluation", "subject prediction", "coefficients")
                     
                     
                     return(out)
                   }
  )
})


# Stop the cluster
stopCluster(cl)
gc() # Running time: 16.460s

save(out, file = "Data/jointLasso_uArray_MCCVType0.Rda")



## Calculate the data summary
load("Data/jointLasso_uArray_MCCVType0.Rda")
library(tibble)
library(tidyr)
library(magrittr)
library(dplyr)

tmp <- lapply(out, function(a) { #out[[1]]
  as.data.frame(a$evaluation)
})

perf <- do.call("rbind", tmp) %>%
  rownames_to_column(var="i")

save(perf, file="Data/jointLasso_uArray_perf_MCCVType0.Rda")


# Density plot
ggplot(perf, aes(x=BER)) +
  geom_density(adjust=2)



## Lasso coefficients
tmp <- lapply(out, function(a) { #out[[1]]
  as.data.frame(a$coefficients)
})

coef <- do.call("rbind", tmp) 


# Remove those that is always 0
coef <- coef[abs(coef$s0) > 0, ]
coef <- coef[coef$feat != "(Intercept)", ]
rownames(coef) <- NULL

coef <- coef %>%
  select(-s0) %>%
  table %>%
  data.frame() %>%
  set_colnames(c("Feat", "Freq")) 



## Get the genes and genus name
# uArray
load("Data/esetRna_new.Rda") 

fdata <- as(featureData(esetRna), "data.frame")
temp <- subset(fdata, select = c(SYMBOL))

coef$Feat <- gsub("|`", "", coef$Feat)
coef <- merge(coef, temp, by.x="Feat", by.y="row.names", all.x=TRUE) %>%
  rename(name="SYMBOL")


# uBiome
load("Data/uBiome_taxtable.Rda")

coef <- merge(coef, taxTable, by.x="Feat", by.y="row.names", all.x=TRUE) %>%
  mutate(name=ifelse(is.na(name), as.character(Genus), name)) %>%
  select(-Genus) %>%
  mutate(prop=Freq/500)

save(coef, file = "Data/joint_lasso_uArray_coef.Rda")


give.n <- function(x, upper_limit = max(coef.lasso$s0, na.rm = TRUE)*1.15){
  return(data.frame(y = as.numeric(.95*upper_limit),
                    label = paste('n=', 
                                  format(length(x), big.mark = ",", decimal.mark = ".", scientific = FALSE))))
}



## Microbiome frequency of selection
d <- coef[grepl("SVs", coef$Feat), ]
d <- d[order(-d$Freq), ]
d$features <- as.factor(factor(d$features, levels = unique(d$features[order(-d$freq)]), ordered = TRUE))
d.20 <- d[1:20, ]
d.20$ind <- ifelse(d.20$Freq >= 250, "1", "0")


ggplot(data=d.20, aes(x=reorder(name, -Freq), y=Freq)) +
  geom_col(aes(fill=ind)) +
  labs(x="", y="Selection frequency") +
  scale_fill_manual(name = "> 500",
                    labels = c("No", "Yes"),
                    values = c("1" = "tomato3", "0" = "grey54")) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none", 
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(angle=90, hjust=1, size = 12))







#----------------------------------------------------------------------------------------
# PULLING 9695 MICROARRAY + 8 MICROBIOME: no penalisation on microbiome 
#----------------------------------------------------------------------------------------
load("Data/df_lasso_uBiome.Rda")

library(parallel)
library(ggplot2)


# Create a cluster of cores
cl <- makeCluster(getOption("cl.cores", 4))
clusterExport(cl=cl, 
              varlist=c("df.uBiome"))


# Main algorithm
system.time({
  out <- parLapply(cl=cl,
                   X=1:500,
                   fun=function(i) {
                     require(glmnet)
                     require(magrittr)
                     require(mixOmics)
                     require(dplyr)
                     
                     cv.lasso <- NULL
                     while(is(cv.lasso, 'try-error') || is.null(cv.lasso)) {
                       
                     # sample 2/3 of the subjects at random
                     id_keep <- sample(x=rownames(df.uBiome),
                                       size=floor(2*nrow(df.uBiome)/3),
                                       replace=FALSE)
                     
                     # Sample test and train data
                     df.train <- df.uBiome[rownames(df.uBiome) %in% id_keep, ]
                     df.test <- df.uBiome[!(rownames(df.uBiome) %in% id_keep), ]
                     
                     x <- df.train[, names(df.train)!="group"] %>%
                       as.matrix.data.frame()
                     
                     #Unpenalized microbiome data
                     penalty.fctr <- ifelse(grepl("SVs", colnames(x), ignore.case = TRUE), 0, 1)
                     
                     
                     # LASSO
                     cv.lasso <- try(cv.glmnet(x=x, 
                                               y=df.train$group, alpha=1, family="binomial", 
                                               type.measure="class",
                                               penalty.factor=penalty.fctr, standardize=TRUE,
                                               nfolds=5, parallel=TRUE), silent=FALSE)
                     }
                     
                     m1 <- glmnet(x, df.train$group, alpha=1, family="binomial",
                                  penalty.factor=penalty.fctr, standardize=TRUE,
                                  lambda=cv.lasso$lambda.min)
                     
                     x.test <- model.matrix(group ~., df.test)[,-1]
                     pred.out <- as.data.frame(predict(m1, newx=x.test, s="lambda.min", 
                                                       type="class", alpha=1)) %>%
                       rename(pred="s1") %>%
                       mutate(pred=factor(pred, levels=c("Adenoma", "CRC")))
                     
                     # Model performance
                     conf <- caret::confusionMatrix(pred.out$pred, as.factor(df.test$group))
                     
                     df.eval <- data.frame(BER = get.BER(t(conf$table)),
                                           mce = 1-conf$overall[["Accuracy"]],
                                           kappa = conf$overall[["Kappa"]],
                                           specificity = conf$byClass[["Specificity"]],
                                           sensitivity = conf$byClass[["Sensitivity"]],
                                           ppv = conf$byClass[["Pos Pred Value"]],
                                           npv = conf$byClass[["Neg Pred Value"]])
                     
                     # Test subject 
                     subjectPred.temp <- data.frame(pred.out, 
                                                    obs = df.test$group,
                                                    sampleID = rownames(pred.out))
                     
                     coeff <- as.data.frame.matrix(coef(m1))
                     coeff$feat <- rownames(coeff)
                     coeff <- coeff[coeff$s0 != 0, ]
                     
                     out <- list(df.eval, subjectPred.temp, coeff)
                     names(out) <- c("evaluation", "subject prediction", "coefficients")
                     
                     
                     return(out)
                   }
  )
})


# Stop the cluster
stopCluster(cl)
gc() # Running time: 2,290.911s

save(out, file = "Data/joint_lasso_uBiome.Rda")



## Calculate the data summary
load("Data/joint_lasso_uBiome.Rda")
library(tibble)
library(tidyr)
library(magrittr)
library(dplyr)

tmp <- lapply(out, function(a) { #out[[1]]
  as.data.frame(a$evaluation)
})

perf <- do.call("rbind", tmp) %>%
  rownames_to_column(var="i")

save(perf, file="Data/joint_lasso_uBiome_perf.Rda")


# Density plot
ggplot(perf, aes(x=BER)) +
  geom_density(adjust=2)



## Lasso coefficients
tmp <- lapply(out, function(a) { #out[[1]]
  as.data.frame(a$coefficients)
})

coef <- do.call("rbind", tmp) 


# Remove those that is always 0
coef <- coef[abs(coef$s0) > 0, ]
coef <- coef[coef$feat != "(Intercept)", ]
rownames(coef) <- NULL

coef <- coef %>%
  select(-s0) %>%
  table %>%
  data.frame() %>%
  set_colnames(c("Feat", "Freq")) 



## Get the genes and genus name
# uArray
library(Biobase)
load("Data/esetRna_new.Rda") 

fdata <- as(featureData(esetRna), "data.frame")
temp <- subset(fdata, select = c(SYMBOL))

coef$Feat <- gsub("|`", "", coef$Feat)
coef <- merge(coef, temp, by.x="Feat", by.y="row.names", all.x=TRUE) %>%
  rename(name="SYMBOL")




# uBiome
load("Data/uBiome_taxtable.Rda")

coef <- merge(coef, taxTable, by.x="Feat", by.y="row.names", all.x=TRUE) %>%
  mutate(name=ifelse(is.na(name), as.character(Genus), name)) %>%
  select(-Genus) %>%
  mutate(prop=Freq/500)

save(coef, file = "Data/joint_lasso_uBiome_coef.Rda")


give.n <- function(x, upper_limit = max(coef.lasso$s0, na.rm = TRUE)*1.15){
  return(data.frame(y = as.numeric(.95*upper_limit),
                    label = paste('n=', 
                                  format(length(x), big.mark = ",", decimal.mark = ".", scientific = FALSE))))
}



## Microarray frequency of selection
d <- coef[!grepl("SVs", coef$Feat), ]
d <- d[order(-d$Freq), ]
d$features <- as.factor(factor(d$features, levels = unique(d$features[order(-d$freq)]), ordered = TRUE))
d.20 <- d[1:20, ]
d.20$ind <- ifelse(d.20$Freq >= 250, "1", "0")


ggplot(data=d.20, aes(x=reorder(name, -Freq), y=Freq)) +
  geom_col(aes(fill=ind)) +
  labs(x="", y="Selection frequency") +
  scale_fill_manual(name = "> 500",
                    labels = c("No", "Yes"),
                    values = c("1" = "tomato3", "0" = "grey54")) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none", 
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(angle=90, hjust=1, size = 12))







#----------------------------------------------------------------------------------------
# LASSO 9,695 MICROARRAY + 431 MICROBIOME, NO PENALISATION
#----------------------------------------------------------------------------------------
load("Data/df_lasso_all.Rda")

library(parallel)
library(ggplot2)


# Create a cluster of cores
cl <- makeCluster(getOption("cl.cores", 4))
clusterExport(cl=cl, 
              varlist=c("df.all"))


# Main algorithm
system.time({
  out <- parLapply(cl=cl,
                   X=1:500,
                   fun=function(i) {
                     require(glmnet)
                     require(magrittr)
                     require(mixOmics)
                     require(dplyr)
                     
                     cv.lasso <- NULL
                     while(is(cv.lasso, 'try-error') || is.null(cv.lasso)) {
                       
                       # sample 2/3 of the subjects at random
                       id_keep <- sample(x=rownames(df.all),
                                         size=floor(2*nrow(df.all)/3),
                                         replace=FALSE)
                       
                       # Sample test and train data
                       df.train <- df.all[rownames(df.all) %in% id_keep, ]
                       df.test <- df.all[!(rownames(df.all) %in% id_keep), ]
                       
                       x <- df.train[, names(df.train)!="group"] %>%
                         as.matrix.data.frame()
                       
                       
                       # LASSO
                       cv.lasso <- try(cv.glmnet(x=x, 
                                                 y=df.train$group, alpha=1, family="binomial", 
                                                 type.measure="class",
                                                 standardize=TRUE,
                                                 nfolds=5, parallel=TRUE), silent=FALSE)
                     }
                     
                     m1 <- glmnet(x, df.train$group, alpha=1, family="binomial",
                                  standardize=TRUE,
                                  lambda=cv.lasso$lambda.min)
                     
                     # Predict test set
                     x.test <- df.test[, names(df.test)!="group"] %>%
                       as.matrix.data.frame()
                     pred.out <- as.data.frame(predict(m1, newx=x.test, s="lambda.min", 
                                                       type="class", alpha=1)) %>%
                       rename(pred="s1") %>%
                       mutate(pred=factor(pred, levels=c("Adenoma", "CRC")))
                     
                     # Model performance
                     conf <- caret::confusionMatrix(pred.out$pred, as.factor(df.test$group))
                     
                     df.eval <- data.frame(BER = get.BER(t(conf$table)),
                                           mce = 1-conf$overall[["Accuracy"]],
                                           kappa = conf$overall[["Kappa"]],
                                           specificity = conf$byClass[["Specificity"]],
                                           sensitivity = conf$byClass[["Sensitivity"]],
                                           ppv = conf$byClass[["Pos Pred Value"]],
                                           npv = conf$byClass[["Neg Pred Value"]])
                     
                     # Test subject 
                     subjectPred.temp <- data.frame(pred.out, 
                                                    obs = df.test$group,
                                                    sampleID = rownames(pred.out))
                     
                     coeff <- as.data.frame.matrix(coef(m1))
                     coeff$feat <- rownames(coeff)
                     coeff <- coeff[coeff$s0 != 0, ]
                     
                     out <- list(df.eval, subjectPred.temp, coeff)
                     names(out) <- c("evaluation", "subject prediction", "coefficients")
                     
                     
                     return(out)
                   }
  )
})


# Stop the cluster
stopCluster(cl)
gc() # Running time: 86.290s

save(out, file = "Data/joint_lasso_all.Rda")



## Calculate the data summary
load("Data/joint_lasso_all.Rda")
library(tibble)
library(tidyr)
library(magrittr)
library(dplyr)

tmp <- lapply(out, function(a) { #out[[1]]
  as.data.frame(a$evaluation)
})

perf <- do.call("rbind", tmp) %>%
  rownames_to_column(var="i")

save(perf, file="Data/joint_lasso_all_perf.Rda")


# Density plot
ggplot(perf, aes(x=BER)) +
  geom_density(adjust=2)



## Lasso coefficients
tmp <- lapply(out, function(a) { #out[[1]]
  as.data.frame(a$coefficients)
})

coef <- do.call("rbind", tmp) 


# Remove those that is always 0
coef <- coef[abs(coef$s0) > 0, ]
coef <- coef[coef$feat != "(Intercept)", ]
rownames(coef) <- NULL

coef <- coef %>%
  select(-s0) %>%
  table %>%
  data.frame() %>%
  set_colnames(c("Feat", "Freq")) 



## Get the genes and genus name
# uArray
library(Biobase)
load("Data/esetRna_new.Rda") 

fdata <- as(featureData(esetRna), "data.frame")
temp <- subset(fdata, select = c(SYMBOL))

coef$Feat <- gsub("|`", "", coef$Feat)
coef <- merge(coef, temp, by.x="Feat", by.y="row.names", all.x=TRUE) %>%
  rename(name="SYMBOL")




# uBiome
load("Data/uBiome_taxtable.Rda")

coef <- merge(coef, taxTable, by.x="Feat", by.y="row.names", all.x=TRUE) %>%
  mutate(name=ifelse(is.na(name), as.character(Genus), name)) %>%
  select(-Genus) %>%
  mutate(prop=Freq/500)

save(coef, file = "Data/joint_lasso_all_coef.Rda")


give.n <- function(x, upper_limit = max(coef.lasso$s0, na.rm = TRUE)*1.15){
  return(data.frame(y = as.numeric(.95*upper_limit),
                    label = paste('n=', 
                                  format(length(x), big.mark = ",", decimal.mark = ".", scientific = FALSE))))
}



## Microarray frequency of selection
d <- coef[!grepl("SVs", coef$Feat), ]
d <- d[order(-d$Freq), ]
d$features <- as.factor(factor(d$features, levels = unique(d$features[order(-d$freq)]), ordered = TRUE))
d.20 <- d[1:20, ]
d.20$ind <- ifelse(d.20$Freq >= 250, "1", "0")


ggplot(data=d.20, aes(x=reorder(name, -Freq), y=Freq)) +
  geom_col(aes(fill=ind)) +
  labs(x="", y="Selection frequency") +
  scale_fill_manual(name = "> 500",
                    labels = c("No", "Yes"),
                    values = c("1" = "tomato3", "0" = "grey54")) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none", 
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(angle=90, hjust=1, size = 12))

# Microbiome frequency of selection
d <- coef[grepl("SVs", coef$Feat), ]
d <- d[order(-d$Freq), ]
d$ind <- ifelse(d$Freq >= 250, "1", "0")

ggplot(data=d, aes(x=reorder(name, -Freq), y=Freq)) +
  geom_col(aes(fill=ind)) +
  labs(x="", y="Selection frequency") +
  scale_fill_manual(name = "> 500",
                    labels = c("No", "Yes"),
                    values = c("1" = "tomato3", "0" = "grey54")) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none", 
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(angle=90, hjust=1, size = 12))







#----------------------------------------------------------------------------------------
# GRAPH COMPARING 3 LASSO 
#----------------------------------------------------------------------------------------
# LASSO 1: X2|X1
# LASSO 2: X1|X2
# LASSO 3: X1,X2




### BER
eval1 <- get(load("Data/joint_lasso_uArray_perf.Rda"))
eval2 <- get(load("Data/joint_lasso_uBiome_perf.Rda"))
eval3 <- get(load("Data/joint_lasso_all_perf.Rda"))

eval1$model <- "1"
eval2$model <- "2"
eval3$model <- "3"


# Combine all evaluation data
eval.all <- rbind(eval1, eval2)
eval.all <- rbind(eval.all, eval3)

med <- eval.all %>%
  select(i, BER, model) %>%
  group_by(model) %>%
  summarise(grp.mean=mean(BER),
            grp.std=sd(BER))

# Plot
ggplot(eval.all, aes(x=BER, color=model)) +
  geom_density() +
  geom_vline(data = med, aes(xintercept = grp.mean, color = model), linetype = "dashed") +
  annotate("text", x = med$grp.mean+0.02, y = c(15, 12, 15), label = paste0("Mean: ", round(med$grp.mean, 2)), 
           size = 4.5, color = c("orange", "steelblue4", "orangered2"))+
  scale_color_manual(values = c("orange", "steelblue4", "orangered2"),
                     labels = c(bquote(eta~"("~X[2]~"|"~tilde(X)[1]~")"), 
                                bquote(eta~"("~X[1]~"|"~tilde(X)[2]~")"), 
                                bquote(eta~"("~X[1]~","~X[2]~")"))) +
  labs(color = "Model") +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        #legend.position = "none", 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14))+
  theme_bw()
