#####################################################
## STATISTICAL ANALYSIS                            ##
## LASSO: joint                                    ##
## Prepared by: Dea Putri (dputri@its.jnj.com)     ##
#####################################################

setwd("/mnt/exports/shared/home/dputri/PhD/CRC microbiome/drp_crcmicrobiome/Joint analysis")
load("Data/rna_crcad.Rda")
load("Data/group_crcad.Rda")
load("Data/otu_431_clr.Rda")

feat <- get(load("Data/selected feat.Rda"))
genus <- get(load("Data/selected genus.Rda"))

feat <- feat[order(-feat$freq), ]
genus <- genus[order(-genus$freq), ]


# Get top 10 uArray & 4 microbiome features
get.feat <- feat$features[1:14]
get.genus <- genus$genus[1:4]



## Subset uArray and uBiome data
rna.14 <- rna.crcad[, colnames(rna.crcad) %in% as.character(get.feat)]
rna.14 <- cbind(rna.14, group)
otu.clr.4 <- otu.clr[, colnames(otu.clr) %in% as.character(get.genus)]
otu.clr.4 <- cbind(otu.clr.4, group)

save(rna.14, file = "Data/df_lasso_14uarray.Rda")
save(otu.clr.4, file = "Data/df_lasso_4uBiome.Rda")



## Lasso data frame for non-penalized uArray (445 features)
df.uArray <- merge(rna.14, otu.clr, by = "row.names")
rownames(df.uArray) <- df.uArray$Row.names
df.uArray <- subset(df.uArray, select = -c(Row.names))
save(df.uArray, file = "Data/df_lasso_uArray.Rda")



## Lasso data frame for non-penalized uBiome (9,699 features)
df.uBiome <- merge(rna.crcad, otu.clr.4, by = "row.names")
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
# PULLING 13 MICROARRAY + 431 MICROBIOME: no penalisation on microarray 
#----------------------------------------------------------------------------------------
load("Data/df_lasso_uArray.Rda")

fit.lasso <- function(df.comb, fold = 3, iter = 1000, standardized = TRUE) {
  require(glmnet)
  require(limma)
  require(magrittr)
  require(mixOmics)
  require(caret)
  
  out <- vector("list", iter)
  
  for(i.k in 1:iter) {
    set.seed(i.k)
    
    smpl <- createDataPartition(y = df.comb$group, p=0.7, list = FALSE)
    
    df.train <- df.comb[smpl, ]
    df.test <- df.comb[-smpl, ]
    
    x <- model.matrix(group~., df.train)[,-1]
    y <- as.factor(factor(df.train$group, levels = c("Adenoma", "CRC")))
    
    
    #Unpenalized microarray data
    penalty.fctr <- ifelse(grepl("SVs", colnames(x), ignore.case = TRUE), 1, 0)
    
    cv.lasso <- NULL
    while(is(cv.lasso, 'try-error') || is.null(cv.lasso)) {
      cv.lasso <- try(cv.glmnet(x, y, alpha = 1, family = "binomial", type.measure = "class",
                                penalty.factor = penalty.fctr, standardize = standardized,
                                nfolds = fold, parallel = TRUE), silent = FALSE)
    }
    
    mod.combined <- glmnet(x, y, alpha = 1, family = "binomial",
                           penalty.factor = penalty.fctr, standardize = standardized,
                           lambda = cv.lasso$lambda.min)
    
    x.test <- model.matrix(group ~., df.test)[,-1]
    pred.out = as.data.frame(predict(mod.combined, newx = x.test, s = "lambda.min", type = "class",
                                     alpha = 1))
    colnames(pred.out) <- "pred"
    conf <- caret::confusionMatrix(pred.out$pred, as.factor(df.test$group))
    
    df.eval <- data.frame(BER = get.BER(t(conf$table)),
                          mce = 1-conf$overall[["Accuracy"]],
                          kappa = conf$overall[["Kappa"]],
                          specificity = conf$byClass[["Specificity"]],
                          sensitivity = conf$byClass[["Sensitivity"]],
                          ppv = conf$byClass[["Pos Pred Value"]],
                          npv = conf$byClass[["Neg Pred Value"]])
    
    subjectPred.temp <- data.frame(pred.out, 
                                   obs = df.test$group,
                                   sampleID = rownames(pred.out))
    
    coeff <- as.data.frame.matrix(coef(mod.combined))
    coeff$feat <- rownames(coeff)
    coeff <- coeff[coeff$s0 != 0, ]
    
    out[[i.k]] <- list(df.eval, subjectPred.temp, coeff)
    names(out[[i.k]]) <- c("evaluation", "subject prediction", "coefficients")
    names(out)[[i.k]] <- i.k
  }
  return(out)
}


# Not-standardized
set.seed(1)
t2 = Sys.time()
joint.lasso.1 <- fit.lasso(df.comb = df.uArray, iter = 500, standardized = FALSE)
t <- Sys.time()
print((t - t2))
save(joint.lasso.1, file = "Data/joint_lasso_1.Rda")      # time = 1.07 mins





#### PERFORMANCE ASSESSMENT
## Calculate the evaluation statistics
load("Data/joint_lasso_1.Rda")


# Evaluation
eval.joint.1 <- NULL
for(i in 1:length(joint.lasso.1)) {
  temp <- joint.lasso.1[[i]]$evaluation
  temp$iter <- names(joint.lasso.1)[i]
  eval.joint.1 <- rbind(eval.joint.1, temp)
}

ggplot(eval.joint.1, aes(x = BER)) +
  geom_density() +
  geom_vline(aes(xintercept = mean(eval.joint.1$BER)), color="red", linetype="dashed", size=1) +
  annotate("text", x = mean(eval.joint.1$BER)+0.02, y = c(15), label = paste0("Mean: ", round(mean(eval.joint.1$BER), 2)), 
           size = 4.5, color = c("red"))+
  labs(x = "BER",
       y = "Count") +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        #legend.position = "none", 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14))





## Lasso coefficients
coef.lasso <- NULL
for(i in 1:length(joint.lasso.1)) {
  temp <- joint.lasso.1[[i]]$coefficients
  temp$iter <- names(joint.lasso.1)[i]
  coef.lasso <- rbind(coef.lasso, temp)
}


# Remove those that is always 0
coef.lasso <- coef.lasso[abs(coef.lasso$s0) > 0, ]
coef.lasso <- coef.lasso[coef.lasso$feat != "(Intercept)", ]



## Get the feature and genus name
# uArray
load("Data/esetRna_new.Rda") 

fdata <- as(featureData(esetRna), "data.frame")
temp <- subset(fdata, select = c(SYMBOL))

uBiome <- colnames(df.uArray)[grepl("SVs", colnames(df.uArray), ignore.case = TRUE)]

temp.uarray <- coef.lasso[!coef.lasso$feat %in% uBiome, ]
temp.uarray$feat <- gsub("|`", "", temp.uarray$feat)
temp.uarray <- merge(temp.uarray, temp, by.x = "feat", by.y = "row.names", all.x = TRUE)
names(temp.uarray)[names(temp.uarray) == "SYMBOL"] <- "name"

save(temp.uarray, file = "Data/uarraySelect_lasso_1.Rda")


# uBiome
load("Data/Obj_Genus.rda")

taxtable <- as.data.frame(tax_table(Obj_Genus))
temp <- subset(taxtable, select = c(Genus))
taxTable <- temp
save(taxTable, file = "Data/uBiome_taxtable.Rda")

temp.ubiome <- coef.lasso[coef.lasso$feat %in% uBiome, ]
temp.ubiome <- merge(temp.ubiome, temp, by.x = "feat", by.y = "row.names", all.x = TRUE)
temp.ubiome$Genus <- as.factor(factor(temp.ubiome$Genus))
names(temp.ubiome)[names(temp.ubiome) == "Genus"] <- "name"

save(temp.ubiome, file = "Data/ubiomeSelect_lasso_1.Rda")

coef.lasso <- rbind(temp.uarray, temp.ubiome)
coef.lasso$name <- as.factor(as.character(coef.lasso$name))
order.x <- factor(unique(coef.lasso$name), levels = unique(coef.lasso$name))

give.n <- function(x, upper_limit = max(coef.lasso$s0, na.rm = TRUE)*1.15){
  return(data.frame(y = as.numeric(.95*upper_limit),
                    label = paste('n=', 
                                  format(length(x), big.mark = ",", decimal.mark = ".", scientific = FALSE))))
}



## Microbiome frequency of selection
d <- as.data.frame(table(temp.ubiome$name))
colnames(d) <- c("features", "freq")
d <- d[order(-d$freq), ]
d$features <- as.factor(factor(d$features, levels = unique(d$features[order(-d$freq)]), ordered = TRUE))
d.20 <- d[1:20, ]
d.20$ind <- ifelse(d.20$freq >= 250, "1", "0")


ggplot(data=d.20, aes(x = features, y = freq)) +
  geom_col(aes(fill = ind)) +
  labs(x = "", y = "Selection frequency") +
  scale_fill_manual(name = "> 500",
                    labels = c("No", "Yes"),
                    values = c("1" = "tomato3", "0" = "grey54")) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none", 
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(angle=90, hjust=1, size = 12))


# Microbiome features
temp.ubiome$name <- as.factor(factor(temp.ubiome$name, levels = levels(d$features)))
coef.ubiome <- temp.ubiome[temp.ubiome$name %in% d.20$features, ]

give.n <- function(x, upper_limit = max(coef.ubiome$s0, na.rm = TRUE)*1.15){
  return(data.frame(y = as.numeric(.95*upper_limit),
                    label = paste('n=', 
                                  format(length(x), big.mark = ",", decimal.mark = ".", scientific = FALSE))))
}

ggplot(coef.ubiome, aes(x=name, y=s0)) + 
  geom_boxplot() + 
  #geom_jitter(shape = 16, position = position_jitter(0.2)) +
  stat_summary(fun.data = give.n, geom = "text", col="red", hjust = .5, vjust = .9, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylab("Lasso coefficients") +
  xlab("Features") +
  #scale_x_discrete(limits = order.x) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, size = 12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14))



## Top 5 uArray genes profiles on CRC and adenoma patients
d <- as.data.frame(table(coef.lasso$name))
d <- d[order(-d$Freq), ]


# Top 5 uArray genes
get.genes <- as.character(d$Var1[1:5])

load("Data/esetRna_new.Rda") 

fdata <- as(featureData(esetRna), "data.frame")
temp <- as.data.frame(subset(fdata, select = c(SYMBOL)))
temp$genes <- rownames(temp)
temp <- as.data.frame(temp[temp$SYMBOL %in% get.genes, ])

load("Data/rna_crcad.Rda")
load("Data/group_crcad.Rda")

group <- as.data.frame(group)
group$sampleID <- rownames(group)
rna <- rna.crcad[, colnames(rna.crcad) %in% temp$genes]
rna$sampleID <- rownames(rna)

rna.long <- melt(rna, 
                 id.vars = c("sampleID"),
                 variable.name = "genes",
                 value.name = "exprs")
rna.long <- merge(rna.long, group, by = "sampleID", all.x = TRUE)
rna.long <- merge(rna.long, temp, by = "genes", all.x = TRUE)


# Plot
ggplot(rna.long, aes(x = group, y = exprs)) +
  geom_boxplot() +
  facet_wrap(~ SYMBOL) +
  ylab("Expression") +
  xlab("Group") +
  theme(panel.background = element_blank(),
        panel.border = element_blank())







#----------------------------------------------------------------------------------------
# PULLING 9695 MICROARRAY + 4 MICROBIOME: no penalisation on microbiome 
#----------------------------------------------------------------------------------------
load("Data/df_lasso_uBiome.Rda")

fit.lasso <- function(df.comb, fold = 3, iter = 1000, standardized = TRUE) {
  require(glmnet)
  require(limma)
  require(magrittr)
  require(mixOmics)
  require(caret)
  
  out <- vector("list", iter)
  
  for(i.k in 1:iter) {
    set.seed(i.k)
    
    smpl <- df.comb$group %>%
      createDataPartition(p = 0.7, list = FALSE)
    
    df.train <- df.comb[smpl, ]
    df.test <- df.comb[-smpl, ]
    
    x <- model.matrix(group~., df.train)[,-1]
    y <- as.factor(factor(df.train$group, levels = c("Adenoma", "CRC")))
    
    #Unpenalized microbiome data
    penalty.fctr <- ifelse(grepl("SVs", colnames(x), ignore.case = TRUE), 0, 1)
    
    cv.lasso <- NULL
    while(is(cv.lasso, 'try-error') || is.null(cv.lasso)) {
      cv.lasso <- try(cv.glmnet(x, y, alpha = 1, family = "binomial", type.measure = "class",
                                penalty.factor = penalty.fctr, standardize = standardized,
                                nfolds = fold, parallel = TRUE), silent = FALSE)
    }
    
    mod.combined <- glmnet(x, y, alpha = 1, family = "binomial",
                           penalty.factor = penalty.fctr, standardize = standardized,
                           lambda = cv.lasso$lambda.min)
    
    x.test <- model.matrix(group ~., df.test)[,-1]
    pred.out = as.data.frame(predict(mod.combined, newx = x.test, s = "lambda.min", type = "class",
                                     alpha = 1))
    colnames(pred.out) <- "pred"
    conf <- caret::confusionMatrix(pred.out$pred, as.factor(df.test$group))
    
    df.eval <- data.frame(BER = get.BER(t(conf$table)),
                          mce = 1-conf$overall[["Accuracy"]],
                          kappa = conf$overall[["Kappa"]],
                          specificity = conf$byClass[["Specificity"]],
                          sensitivity = conf$byClass[["Sensitivity"]],
                          ppv = conf$byClass[["Pos Pred Value"]],
                          npv = conf$byClass[["Neg Pred Value"]])
    
    subjectPred.temp <- data.frame(pred.out, 
                                   obs = df.test$group,
                                   sampleID = rownames(pred.out))
    
    coeff <- as.data.frame.matrix(coef(mod.combined))
    coeff$feat <- rownames(coeff)
    coeff <- coeff[coeff$s0 != 0, ]
    
    out[[i.k]] <- list(df.eval, subjectPred.temp, coeff)
    names(out[[i.k]]) <- c("evaluation", "subject prediction", "coefficients")
    names(out)[[i.k]] <- i.k
    print(i.k)
  }
  return(out)
}


# Not-standardized
set.seed(1)
t2 = Sys.time()
joint.lasso.2 <- fit.lasso(df.comb = df.uBiome, iter = 500, standardized = FALSE)
t <- Sys.time()
print((t - t2))
save(joint.lasso.2, file = "Data/joint_lasso_2.Rda")      # time = 3.056 hours





#### PERFORMANCE ASSESSMENT
## Calculate the evaluation statistics
load("Data/joint_lasso_2.Rda")


# Evaluation
eval.joint.2 <- NULL
for(i in 1:length(joint.lasso.2)) {
  temp <- joint.lasso.2[[i]]$evaluation
  temp$iter <- names(joint.lasso.2)[i]
  eval.joint.2 <- rbind(eval.joint.2, temp)
}

ggplot(eval.joint.2, aes(x = BER)) +
  geom_density() +
  geom_vline(aes(xintercept = mean(eval.joint.2$BER)), color="red", linetype="dashed", size=1) +
  annotate("text", x = mean(eval.joint.2$BER)+0.02, y = c(15), label = paste0("Mean: ", round(mean(eval.joint.2$BER), 2)), 
           size = 4.5, color = c("red"))+
  labs(x = "BER",
       y = "Count") +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        #legend.position = "none", 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14))



## Lasso coefficients
coef.lasso <- NULL
for(i in 1:length(joint.lasso.2)) {
  temp <- joint.lasso.2[[i]]$coefficients
  
  if(nrow(temp) == 0) {
    next
  } else {
    temp$iter <- names(joint.lasso.2)[i]
    coef.lasso <- rbind(coef.lasso, temp)
  }
  
  print(i)
}


# Remove those that is always 0
coef.lasso <- coef.lasso[abs(coef.lasso$s0) > 0, ]
coef.lasso <- coef.lasso[coef.lasso$feat != "(Intercept)", ]



## Get the feature and genus name
ubiome <- unique(coef.lasso$feat[grepl("svs", coef.lasso$feat, ignore.case = TRUE)])

# uArray
load("Data/uArray_featTable.Rda") 

temp.uarray <- coef.lasso[!coef.lasso$feat %in% ubiome, ]
temp.uarray$feat <- gsub("|`", "", temp.uarray$feat)
temp.uarray <- merge(temp.uarray, temp, by.x = "feat", by.y = "row.names", all.x = TRUE)
names(temp.uarray)[names(temp.uarray) == "SYMBOL"] <- "name"

save(temp.uarray, file = "Data/uarraySelect_lasso_2.Rda")


# uBiome
load("Data/uBiome_taxtable.Rda")
taxtable <- subset(taxtable, select = -c(Kingdom, Phylum, Class, Order, Family))

temp.ubiome <- coef.lasso[coef.lasso$feat %in% ubiome, ]
temp.ubiome <- merge(temp.ubiome, taxtable, by.x = "feat", by.y = "row.names", all.x = TRUE)
temp.ubiome$Genus <- as.factor(factor(temp.ubiome$Genus))
names(temp.ubiome)[names(temp.ubiome) == "Genus"] <- "name"

save(temp.ubiome, file = "Data/ubiomeSelect_lasso_2.Rda")

coef.lasso <- rbind(temp.uarray, temp.ubiome)
coef.lasso$name <- as.factor(as.character(coef.lasso$name))
order.x <- factor(unique(coef.lasso$name), levels = unique(coef.lasso$name))

give.n <- function(x, upper_limit = max(coef.lasso$s0, na.rm = TRUE)*1.15){
  return(data.frame(y = as.numeric(.95*upper_limit),
                    label = paste('n=', 
                                  format(length(x), big.mark = ",", decimal.mark = ".", scientific = FALSE))))
}



## Microarray frequency of selection
d <- as.data.frame(table(temp.uarray$name))
colnames(d) <- c("features", "freq")
d <- d[order(-d$freq), ]
d$features <- as.factor(factor(d$features, levels = unique(d$features[order(-d$freq)]), ordered = TRUE))
d$prop <- d$freq/500
d.20 <- d[1:20, ]
d.20$ind <- ifelse(d.20$prop >= .5, "1", "0")

ggplot(data=d.20, aes(x = features, y = prop)) +
  geom_col(aes(fill = ind)) +
  labs(x = "", y = "Selection proportion") +
  scale_fill_manual(name = "> 0.5",
                    labels = c("No", "Yes"),
                    values = c("1" = "tomato3", "0" = "grey54")) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none", 
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(angle=90, hjust=1, size = 12))


# Microarray features
temp.uarray.20 <- temp.uarray[temp.uarray$name %in% d.20$features, ]
temp.uarray.20$name <- as.factor(factor(temp.uarray.20$name, levels = levels(d.20$features)))
coef.uarray <- temp.uarray.20[temp.uarray.20$name %in% d.20$features, ]

give.n <- function(x, upper_limit = max(coef.uarray$s0, na.rm = TRUE)*1.15){
  return(data.frame(y = as.numeric(.95*upper_limit),
                    label = paste('n=', 
                                  format(length(x), big.mark = ",", decimal.mark = ".", scientific = FALSE))))
}

ggplot(coef.uarray, aes(x=name, y=s0)) + 
  geom_boxplot() + 
  #geom_jitter(shape = 16, position = position_jitter(0.2)) +
  stat_summary(fun.data = give.n, geom = "text", col="red", hjust = .5, vjust = .9, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  ylab("Lasso coefficients") +
  xlab("Features") +
  #scale_x_discrete(limits = order.x) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, size = 12),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14))



## Top 5 uArray genes profiles on CRC and adenoma patients
d <- as.data.frame(table(coef.lasso$name))
d <- d[order(-d$Freq), ]


# Top 5 uArray genes
get.genes <- as.character(d$Var1[1:5])

load("Data/esetRna_new.Rda") 

fdata <- as(featureData(esetRna), "data.frame")
temp <- as.data.frame(subset(fdata, select = c(SYMBOL)))
temp$genes <- rownames(temp)
temp <- as.data.frame(temp[temp$SYMBOL %in% get.genes, ])

load("Data/rna_crcad.Rda")
load("Data/group_crcad.Rda")

group <- as.data.frame(group)
group$sampleID <- rownames(group)
rna <- rna.crcad[, colnames(rna.crcad) %in% temp$genes]
rna$sampleID <- rownames(rna)

rna.long <- melt(rna, 
                 id.vars = c("sampleID"),
                 variable.name = "genes",
                 value.name = "exprs")
rna.long <- merge(rna.long, group, by = "sampleID", all.x = TRUE)
rna.long <- merge(rna.long, temp, by = "genes", all.x = TRUE)


# Plot
ggplot(rna.long, aes(x = group, y = exprs)) +
  geom_boxplot() +
  facet_wrap(~ SYMBOL) +
  ylab("Expression") +
  xlab("Group") +
  theme(panel.background = element_blank(),
        panel.border = element_blank())







#----------------------------------------------------------------------------------------
# LASSO 9,695 MICROARRAY + 431 MICROBIOME, NO PENALISATION
#----------------------------------------------------------------------------------------
load("Data/df_lasso_all.Rda")

fit.lasso <- function(df.comb, fold = 3, iter = 1000, standardized = FALSE) {
  require(glmnet)
  require(limma)
  require(magrittr)
  require(mixOmics)
  require(caret)
  
  out <- vector("list", iter)
  
  for(i.k in 1:iter) {
    set.seed(i.k)
    smpl <- df.comb$group %>%
      createResample(times = 1, list = FALSE)
    
    df.train <- df.comb[smpl, ]
    df.test <- df.comb[-smpl, ]
    
    x <- model.matrix(group~., df.train)[,-1]
    y <- as.factor(factor(df.train$group, levels = c("Adenoma", "CRC")))
    
    cv.lasso <- try(cv.glmnet(x, y, alpha = 1, family = "binomial", type.measure = "class",
                              nfolds = fold, parallel = TRUE, standardize = standardized), silent = FALSE)
    
    if(is(cv.lasso, 'try-error')) {
      while(is(cv.lasso, 'try-error')) {
        cv.lasso <- try(cv.glmnet(x, y, alpha = 1, family = "binomial", type.measure = "class",
                                  nfolds = fold, parallel = TRUE, standardize = standardized), silent = FALSE)
      }
    }
    
    mod.combined <- glmnet(x, y, alpha = 1, family = "binomial", standardize = standardized,
                           lambda = cv.lasso$lambda.min)
    
    x.test <- model.matrix(group ~., df.test)[,-1]
    
    # Evaluation
    pred.out = as.data.frame(predict(mod.combined, newx = x.test, s = "lambda.min", type = "class",
                                     alpha = 1))
    colnames(pred.out) <- "pred"
    conf <- caret::confusionMatrix(pred.out$pred, as.factor(df.test$group))
      
    df.eval <- data.frame(BER = get.BER(t(conf$table)),
                       mce = 1-conf$overall[["Accuracy"]],
                       kappa = conf$overall[["Kappa"]],
                       specificity = conf$byClass[["Specificity"]],
                       sensitivity = conf$byClass[["Sensitivity"]],
                       ppv = conf$byClass[["Pos Pred Value"]],
                       npv = conf$byClass[["Neg Pred Value"]])
    
    # Subject prediction
    subjectPred.temp <- data.frame(pred.out, 
                                   obs = df.test$group,
                                   sampleID = rownames(pred.out))
    
    # Coefficients
    coeff <- as.data.frame.matrix(coef(mod.combined))
    coeff$feat <- rownames(coeff)
    coeff <- coeff[coeff$s0 != 0, ]
    
    out[[i.k]] <- list(df.eval, subjectPred.temp, coeff)
    names(out[[i.k]]) <- c("evaluation", "subject prediction", "coefficients")
    names(out)[[i.k]] <- i.k
    print(i.k)
  }
  return(out)
}


# Not-standardized
set.seed(1)
t2 = Sys.time()
joint.lasso.3 <- fit.lasso(df = df.all, iter = 500, standardized = FALSE)
t <- Sys.time()
print((t - t2))
save(joint.lasso.3, file = "Data/joint_lasso_3.Rda")      # time = 3.05 hours





#### PERFORMANCE ASSESSMENT
## Calculate the evaluation statistics
# Evaluation
load("Data/joint_lasso_3.Rda")
eval.joint.3 <- NULL
for(i in 1:length(joint.lasso.3)) {
  temp <- joint.lasso.3[[i]]$evaluation
  temp$iter <- names(joint.lasso.3)[i]
  eval.joint.3 <- rbind(eval.joint.3, temp)
}

ggplot(eval.joint.3, aes(x = BER)) +
  geom_density() +
  geom_vline(aes(xintercept = mean(eval.joint.3$BER)), color="red", linetype="dashed", size=1) +
  annotate("text", x = mean(eval.joint.3$BER)+0.02, y = c(15), label = paste0("Mean: ", round(mean(eval.joint.3$BER), 2)), 
           size = 4.5, color = c("red"))+
  labs(x = "BER",
       y = "Count") +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        #legend.position = "none", 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14))



## Lasso coefficients
coef.lasso <- NULL
for(i in 1:length(joint.lasso.3)) {
  temp <- joint.lasso.3[[i]]$coefficients
  
  if(nrow(temp) == 0) {
    next
  } else {
    temp$iter <- names(joint.lasso.3)[i]
    coef.lasso <- rbind(coef.lasso, temp)
  }
  
  print(i)
}


# Remove those that is always 0
coef.lasso <- coef.lasso[abs(coef.lasso$s0) > 0, ]
coef.lasso <- coef.lasso[coef.lasso$feat != "(Intercept)", ]



## Get the feature and genus name
ubiome <- unique(coef.lasso$feat[grepl("svs", coef.lasso$feat, ignore.case = TRUE)])


# uArray
load("Data/uArray_featTable.Rda") 
fdata <- subset(fdata, select = SYMBOL)

temp.uarray <- coef.lasso[!coef.lasso$feat %in% ubiome, ]
temp.uarray$feat <- gsub("|`", "", temp.uarray$feat)
temp.uarray <- merge(temp.uarray, fdata, by.x = "feat", by.y = "row.names", all.x = TRUE)
names(temp.uarray)[names(temp.uarray) == "SYMBOL"] <- "name"

save(temp.uarray, file = "Data/uarraySelect_lasso_3.Rda")


# uBiome
load("Data/uBiome_taxtable.Rda")
temp.ubiome <- coef.lasso[coef.lasso$feat %in% ubiome, ]
temp.ubiome <- merge(temp.ubiome, taxTable, by.x = "feat", by.y = "row.names", all.x = TRUE)
temp.ubiome$Genus <- as.factor(factor(temp.ubiome$Genus))
names(temp.ubiome)[names(temp.ubiome) == "Genus"] <- "name"

save(temp.ubiome, file = "Data/ubiomeSelect_lasso_3.Rda")

coef.lasso <- rbind(temp.uarray, temp.ubiome)
coef.lasso$name <- as.factor(as.character(coef.lasso$name))
order.x <- factor(unique(coef.lasso$name), levels = unique(coef.lasso$name))

give.n <- function(x, upper_limit = max(coef.lasso$s0, na.rm = TRUE)*1.15){
  return(data.frame(y = as.numeric(.95*upper_limit),
                    label = paste('n=', 
                                  format(length(x), big.mark = ",", decimal.mark = ".", scientific = FALSE))))
}



## Microarray frequency of selection
d <- as.data.frame(table(temp.uarray$name))
colnames(d) <- c("features", "freq")
d <- d[order(-d$freq), ]
d$features <- as.factor(factor(d$features, levels = unique(d$features[order(-d$freq)]), ordered = TRUE))
d$prop <- d$freq/500
d.20 <- d[1:20, ]
d.20$ind <- ifelse(d.20$prop >= .5, "1", "0")

ggplot(data=d.20, aes(x = features, y = prop)) +
  geom_col(aes(fill = ind)) +
  labs(x = "", y = "Selection proportion") +
  scale_fill_manual(name = "> 0.5",
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
load("Data/joint_lasso_1.Rda")
load("Data/joint_lasso_2.Rda")
load("Data/joint_lasso_3.Rda")

eval.joint.1 <- NULL
for(i in 1:length(joint.lasso.1)) {
  temp <- joint.lasso.1[[i]]$evaluation
  temp$iter <- names(joint.lasso.1)[i]
  eval.joint.1 <- rbind(eval.joint.1, temp)
}
eval.joint.1$model <- "1"

eval.joint.2 <- NULL
for(i in 1:length(joint.lasso.2)) {
  temp <- joint.lasso.2[[i]]$evaluation
  temp$iter <- names(joint.lasso.2)[i]
  eval.joint.2 <- rbind(eval.joint.2, temp)
}
eval.joint.2$model <- "2"

eval.joint.3 <- NULL
for(i in 1:length(joint.lasso.3)) {
  temp <- joint.lasso.3[[i]]$evaluation
  temp$iter <- names(joint.lasso.3)[i]
  eval.joint.3 <- rbind(eval.joint.3, temp)
}
eval.joint.3$model <- "3"

eval.all <- rbind(eval.joint.1, eval.joint.2)
eval.all <- rbind(eval.all, eval.joint.3)

ber <- subset(eval.all, select = c(BER, model))
med <- ddply(ber, .(model), summarise, grp.mean = mean(BER), grp.std = sd(BER))

ggplot(ber, aes(x = BER, color = model)) +
  geom_density() +
  geom_vline(data = med, aes(xintercept = grp.mean, color = model), linetype = "dashed") +
  annotate("text", x = med$grp.mean+0.02, y = c(15, 12, 15), label = paste0("Mean: ", round(med$grp.mean, 2)), 
           size = 4.5, color = c("orange", "steelblue4", "orangered2"))+
  scale_color_manual(values = c("orange", "steelblue4", "orangered2"),
                     labels = c(bquote(eta~"("~X[2]~"|"~tilde(X)[1]~")"), 
                                bquote(eta~"("~X[1]~"|"~tilde(X)[2]~")"), 
                                bquote(eta~"("~X[1]~","~X[2]~")"))) +
  labs(x = "MCE",
       y = "Count",
       color = "Model") +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        #legend.position = "none", 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14))+
  theme_bw()
