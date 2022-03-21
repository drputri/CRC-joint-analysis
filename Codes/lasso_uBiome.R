#####################################################
## STATISTICAL ANALYSIS                            ##
## LASSO: uBiome                                   ##
## Prepared by: Dea Putri (dputri@its.jnj.com)     ##
#####################################################

setwd("/home/dputri/02. PhD/CRC microbiome/drp_crcmicrobiome/Joint analysis/")







#----------------------------------------------------------------------------------------
# FEATURE SELECTION: 431 genera -- clr transformed
#----------------------------------------------------------------------------------------
library(compositions)

load("Data/uBiome_OTU_GMPR.Rda")


# Order the subject ID in otu and group
otu <- otu[order(row.names(otu)), ]


# Remove unclassified genera
load("Data/Obj_Genus.rda")

taxtable <- as.data.frame(tax_table(Obj_Genus))
rmv <- rownames(taxtable)[taxtable$Genus == "unclassified"]

otu <- otu[, !colnames(otu) %in% rmv]
otu.t <- data.frame(t(otu))
otu.clr <- data.frame(t(clr(otu.t)))
rownames(otu.clr) <- sub(".", "", rownames(otu.clr))   # removing the X in the first rownames 

save(otu.clr, file = "Data/otu_431_clr.Rda")


# Using CMA
load("Data/otu_431_clr.Rda")
load("Data/group_crcad.Rda")
source("Codes/enet_functions.R")

compute.local <- function() {
  out.all <- vector("list", 1000)
  for(i in 1:1000) {
    out.all[[i]] <- select_feat_CMA(X = otu.clr, y = group, fold = 3, seed = i)
  }
  return(out.all)
}



## Create the cluster
library(doSNOW)
library(parallel)
library(snow)
library(foreach)

message("...Start to run glmnet...")
t2 = Sys.time()

cl <- makeCluster(4, type = "SOCK")
clusterEvalQ(cl, library(data.table))
clusterEvalQ(cl, library(CMA))

registerDoSNOW(cl)

genus <- compute.local()
stopCluster(cl)

registerDoSEQ()
t <- Sys.time()
print((t - t2))
save(genus, file = "Data/feat_uBiome_lasso_CMA.Rda")      # 1.942011 mins







#----------------------------------------------------------------------------------------
# FEATURES SELECTED: 431 genera - clr transformed
#----------------------------------------------------------------------------------------
## Graph: features selected
load("Data/feat_uBiome_lasso_CMA.Rda")
load("Data/otu_431_clr.Rda")
load("Data/uBiome_taxtable.Rda")

genus <- do.call("rbind", genus)
genus <- ddply(genus, .(index), summarize, "mean" = mean(importance), "freq" = length(importance))
genus$genus <- colnames(otu.clr)[genus$index]
genus <- subset(genus, select = -index)

d <- merge(genus, taxTable, by.x = "genus", by.y = "row.names", all.x = TRUE)
d$freq <- d$freq/3
d$ind <- ifelse(d$prop >= .5, "1", "0")

save(d, file = "Data/selected genus.Rda")


d$prop <- d$freq/1000
d <- d[order(-d$prop), ]

ggplot(data=d[1:20, ], aes(x = reorder(Genus, -prop), y = prop)) +
  geom_col(aes(fill = ind)) +
  scale_fill_manual(name = "> 500",
                    labels = c("No", "Yes"),
                    values = c("1" = "tomato3", "0" = "grey54")) +
  labs(x = "Genus", y = "Frequency", title = "Number of repeats: 1,000") +
  ylim(0, 1) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none", 
        axis.text.x = element_text(angle=90, hjust=1, size = 8))







#----------------------------------------------------------------------------------------
# FIXED GENERA: 3-fold cross validation, 431 genera - clr transformed
#----------------------------------------------------------------------------------------
load("Data/group_crcad.Rda")
load("Data/otu_431_clr.Rda")
load("Data/selected genus.Rda")
source("Codes/enet_functions.R")



## Cluster
library(doSNOW)
library(parallel)
library(snow)
library(foreach)
library(magrittr)


# Input 
k <- c(seq(2, nrow(d[d$ind == 1, ])))

compute.local <- function() {
  out.all <- vector("list", length(k))
  for(i in seq_along(k)) {
    out.all[[i]] <- enet_eval_CMA(k = k[i], X = otu.clr, y = group, df.feat = d, feat.name = "genus")
    names(out.all)[[i]] <- paste0(k[i])
  }
  return(out.all)
}



## Create the cluster
message("...Start to run glmnet...")
t2 = Sys.time()

cl <- makeCluster(4, type = "SOCK")
clusterExport(cl, c("StratSampleSelection"), envir=globalenv())
clusterEvalQ(cl, library(data.table))
clusterEvalQ(cl, library(magrittr))
clusterEvalQ(cl, library(glmnet))
clusterEvalQ(cl, library(sampling))
clusterEvalQ(cl, library(mixOmics))
clusterEvalQ(cl, library(caret))

registerDoSNOW(cl)

set.seed(1)
out.list <- compute.local()
stopCluster(cl)

registerDoSEQ()
t <- Sys.time()
print((t - t2))
save(out.list, file = "Data/uBiome_eval_lasso_431_clr.Rda")       #32.67 mins








#----------------------------------------------------------------------------------------
# PERFORMANCE ASSESSMENT: 431 genera - clr transformed
#----------------------------------------------------------------------------------------
## Calculate the data summary
load("Data/uBiome_eval_lasso_431_clr.Rda")

eval <- NULL
for(i in 1:length(out.list)) {
  x <- as.data.frame(out.list[[i]])
  x$k <- names(out.list)[i]
  
  if(is.null(eval)) {
    eval <- x
  } else {
    eval <- rbind(eval, x)
  }
}

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

eval$k <- as.factor(eval$k)
sorted_labels <- paste(sort(as.integer(levels(eval$k))))
eval$k <- factor(eval$k, levels = sorted_labels)

ber <- subset(eval, select = c(BER, k))
med <- ddply(ber, .(k), summarise, grp.mean = mean(BER), grp.sd = sd(BER))

ggplot(ber, aes(x = BER, color = k)) +
  geom_density() +
  geom_vline(data = med, aes(xintercept = grp.mean, color = k), linetype = "dashed") +
  annotate("text", x = med$grp.mean+0.02, y = c(15, 12, 13.5), label = paste0("Mean: ", round(med$grp.mean, 2)), 
           size = 4.8, color = c("orange", "steelblue4", "orangered2"))+
  scale_color_manual(values = c("orange", "steelblue4", "orangered2")) +
  labs(x = "MCE",
       y = "Count") +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        #legend.position = "none", 
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14))







#----------------------------------------------------------------------------------------
# FEATURE SELECTION: proportion of zero < 70% in one of the group
#----------------------------------------------------------------------------------------
load("Data/uBiome_OTU_GMPR.Rda")
load("Data/group_crcad.Rda")
load("Data/prop_zero_wide.Rda")
source("Codes/enet_functions.R")

temp <- prop_zero
temp$prop.AdenomaG <- ifelse(temp$prop.AdenomaG < .7, 1, 0)
temp$prop.CRC <- ifelse(temp$prop.CRC < .7, 1, 0)
temp$keep <- rowSums(temp[, c(1:2)])
genusToKeep <- temp$genus[which(temp$keep > 1)]


# Order the subject ID in otu and group
otu <- otu[order(row.names(otu)), ]


# Remove unclassified genera
load("Data/Obj_Genus.rda")

taxtable <- as.data.frame(tax_table(Obj_Genus))
rmv <- rownames(taxtable)[taxtable$Genus == "unclassified"]

otu <- otu[, !colnames(otu) %in% rmv]


# Remove genera that the proportion of zero > 70% in the two groups
otu <- otu[, which(names(otu) %in% genusToKeep)]

save(otu, file = "Data/otu_61.Rda")


compute.local <- function() {
  out.all <- vector("list", 1000)
  for(i in 1:1000) {
    out.all[[i]] <- select_feature(X = otu, y = group, fold = 3, alpha = 1, seed = i)
  }
  return(out.all)
}



## Create the cluster
library(doSNOW)
library(parallel)
library(snow)
library(foreach)

message("...Start to run glmnet...")
t2 = Sys.time()

cl <- makeCluster(4, type = "SOCK")
clusterEvalQ(cl, library(data.table))
clusterEvalQ(cl, library(glmnet))

registerDoSNOW(cl)

genus <- compute.local()
stopCluster(cl)

registerDoSEQ()
t <- Sys.time()
print((t - t2))
save(genus, file = "Data/feat_uBiome_lasso_61.Rda")







#----------------------------------------------------------------------------------------
# FEATURES SELECTED: 61 genera
#----------------------------------------------------------------------------------------
load("Data/feat_uBiome_lasso_61.Rda")
load("Data/Obj_Genus.rda")

taxtable <- as.data.frame(tax_table(Obj_Genus))
temp <- subset(taxtable, select = c(Family, Genus))

genus <- as.data.frame(unlist(genus))
colnames(genus) <- "genus"

d <- as.data.frame(table(genus$genus))
colnames(d) <- c("genus", "freq")
d <- merge(d, temp, by.x = "genus", by.y = "row.names", all.x = TRUE)
d$ind <- ifelse(d$freq >= 500, "1", "0")

save(d, file = "Data/selected genus_61.Rda")

ggplot(data=d, aes(x = reorder(Genus, -freq), y = freq)) +
  geom_col(aes(fill = ind)) +
  scale_fill_manual(name = "> 500",
                    labels = c("No", "Yes"),
                    values = c("1" = "tomato3", "0" = "grey54")) +
  labs(x = "Genus", y = "Frequency", title = "Genera selected") +
  ylim(0, 1000) +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none", 
        axis.text.x = element_text(angle=90, hjust=1, size = 8))







#----------------------------------------------------------------------------------------
# FIXED GENERA: 3-fold cross validation, 61 genera
#----------------------------------------------------------------------------------------
load("Data/group_crcad.Rda")
load("Data/otu_61.Rda")
load("Data/selected genus_61.Rda")
source("Codes/enet_functions.R")



## Cluster
library(doSNOW)
library(parallel)
library(snow)
library(foreach)
library(magrittr)


# Input 
k <- c(seq(2, nrow(d[d$ind == 1, ])))
X = otu
y = group
y <- as.data.frame(y)
colnames(y) <- "group"


compute.local <- function() {
  out.all <- vector("list", length(k))
  for(i in seq_along(k)) {
    out.all[[i]] <- enet_eval(k = k[i], X = X, y = y, group.name = "group", df.feat = d, 
                              df.feat.all = TRUE, feat.name = "genus")
    names(out.all)[[i]] <- paste0(k[i])
  }
  return(out.all)
}



## Create the cluster
message("...Start to run glmnet...")
t2 = Sys.time()

cl <- makeCluster(4, type = "SOCK")
clusterExport(cl, c("StratSampleSelection"), envir=globalenv())
clusterEvalQ(cl, library(data.table))
clusterEvalQ(cl, library(magrittr))
clusterEvalQ(cl, library(glmnet))
clusterEvalQ(cl, library(sampling))
clusterEvalQ(cl, library(mixOmics))
clusterEvalQ(cl, library(caret))

registerDoSNOW(cl)

set.seed(1)
out.list <- compute.local()
stopCluster(cl)

registerDoSEQ()
t <- Sys.time()
print((t - t2))
save(out.list, file = "Data/uBiome_eval_lasso_61.Rda")







#----------------------------------------------------------------------------------------
# PERFORMANCE ASSESSMENT: 61 genera
#----------------------------------------------------------------------------------------
## Calculate the data summary
load("Data/uBiome_eval_lasso_61.Rda")

eval <- NULL
for(i in 1:length(out.list)) {
  temp <- out.list[[i]]
  for(j in 1:length(temp)) {
    x <- data.frame(temp[[j]]$evaluation)
    x$k <- names(out.list)[i]
    eval <- rbind(eval, x)
  }
}

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

eval$k <- as.factor(eval$k)
sorted_labels <- paste(sort(as.integer(levels(eval$k))))
eval$k <- factor(eval$k, levels = sorted_labels)

library(reshape2)
eval_long <- melt(eval, id.vars = c("k"))
eval_long <- data_summary(eval_long, varname = "value", 
                          groupnames = c("k", "variable"))

pd <- position_dodge(0.1)
ggplot(eval_long, aes(x=k, y=value)) + 
  facet_wrap(~ variable) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.1, position=pd) +
  geom_line(position = pd) +
  geom_point(position = pd) +
  theme(panel.background = element_blank(),
        panel.border = element_blank()#,
        #legend.position = "none", 
        #axis.text.x = element_text(angle=90, hjust=1, size = 8)
  )

