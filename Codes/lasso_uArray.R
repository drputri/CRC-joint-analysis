#####################################################
## STATISTICAL ANALYSIS                            ##
## LASSO: uArray                                   ##
## Prepared by: Dea Putri (dputri@its.jnj.com)     ##
#####################################################

setwd("/home/dputri/02. PhD/CRC microbiome/drp_crcmicrobiome/Joint analysis/")
load("Data/rna_crcad.Rda")
load("Data/group_crcad.Rda")
source("Codes/enet_functions.R")


# Using CMA
compute.local <- function() {
  out.all <- vector("list", 1000)
  for(i in 1:1000) {
    out.all[[i]] <- select_feat_CMA(X = rna.crcad, y = group, fold = 3, seed = i)
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

feat <- compute.local()
stopCluster(cl)

registerDoSEQ()
t <- Sys.time()
print((t - t2))
save(feat, file = "Data/feat_uArray_lasso.Rda")







#----------------------------------------------------------------------------------------
# FEATURES SELECTED
#----------------------------------------------------------------------------------------
library(Biobase)


## Graph: features selected
load("Data/feat_uArray_lasso.Rda")
load("Data/esetRna_new.Rda") 

fdata <- as(featureData(esetRna), "data.frame")
temp <- subset(fdata, select = c(SYMBOL))

feat <- as.data.frame(unlist(feat))
colnames(feat) <- "features"

d <- as.data.frame(table(feat$features))
colnames(d) <- c("features", "freq")
d <- merge(d, temp, by.x = "features", by.y = "row.names", all.x = TRUE)
d$ind <- ifelse(d$freq >= 500, "1", "0")

save(d, file = "Data/selected feat.Rda")

ggplot(data=d[d$freq > 50, ], aes(x = reorder(SYMBOL, -prop), y = prop)) +
  geom_col(aes(fill = ind)) +
  scale_fill_manual(name = "> 500",
                    labels = c("No", "Yes"),
                    values = c("1" = "tomato3", "0" = "grey54")) +
  labs(x = "Features", y = "Frequency", title = "Number of repeats: 1,000") +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "none", 
        axis.text.x = element_text(angle=90, hjust=1, size = 8))







#----------------------------------------------------------------------------------------
# FIXED GENES: 3-fold cross validation
#----------------------------------------------------------------------------------------
load("Data/selected feat.Rda")

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
    out.all[[i]] <- enet_eval_CMA(k = k[i], X = rna.crcad, y = group, df.feat = d, feat.name = "features")
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
clusterEvalQ(cl, library(CMA))
clusterEvalQ(cl, library(mixOmics))
clusterEvalQ(cl, library(caret))

registerDoSNOW(cl)

set.seed(1)
out.list <- compute.local()
stopCluster(cl)

registerDoSEQ()
t <- Sys.time()
print((t - t2))
save(out.list, file = "Data/uArray_eval_lasso.Rda")      #2.426232 hours







#----------------------------------------------------------------------------------------
# PERFORMANCE ASSESSMENT
#----------------------------------------------------------------------------------------
## Calculate the data summary
load("Data/uArray_eval_lasso.Rda")

eval <- NULL
for(i in 1:length(out.list)) {
  temp <- as.data.frame(out.list[[i]])
  temp$k <- names(out.list)[i]
  eval <- rbind(eval, temp)
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
ggplot(eval_long[eval_long$variable == "BER",], aes(x=k, y=value)) + 
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.1, position=pd) +
  #geom_jitter(data = eval[, c("BER", "k")], aes(x = k, y = BER), color = "lightblue", alpha = 0.5, width = 0.1)
  geom_line(position = pd) +
  geom_point(position = pd) +
  ylim(0, .2) +
  labs(y = "BER") +
  theme(panel.background = element_blank(),
        panel.border = element_blank()#,
        #legend.position = "none", 
        #axis.text.x = element_text(angle=90, hjust=1, size = 8)
  )


# Other performance metrices (MCE, SEN, SPE)
library(gridExtra)
p1 <- ggplot(eval_long[eval_long$variable == 'mce',], aes(x=k, y=value)) + 
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.1, position=pd) +
  #geom_jitter(data = eval[, c("BER", "k")], aes(x = k, y = BER), color = "lightblue", alpha = 0.5, width = 0.1)
  geom_line(position = pd) +
  geom_point(position = pd) +
  ylim(0, .2) +
  labs(title='MCE') +
  theme(panel.background = element_blank(),
        panel.border = element_blank()#,
        #legend.position = "none", 
        #axis.text.x = element_text(angle=90, hjust=1, size = 8)
  )

p2 <- ggplot(eval_long[eval_long$variable == 'sensitivity',], aes(x=k, y=value)) + 
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.1, position=pd) +
  geom_line(position = pd) +
  geom_point(position = pd) +
  ylim(0.95, 1.01) +
  labs(title='Sensitivity (SEN)') +
  theme(panel.background = element_blank(),
        panel.border = element_blank()#,
        #legend.position = "none", 
        #axis.text.x = element_text(angle=90, hjust=1, size = 8)
  )

p3 <- ggplot(eval_long[eval_long$variable == 'specificity',], aes(x=k, y=value)) + 
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.1, position=pd) +
  geom_line(position = pd) +
  geom_point(position = pd) +
  ylim(0.8, 1) +
  labs(title='Specificity (SPE)') +
  theme(panel.background = element_blank(),
        panel.border = element_blank()#,
        #legend.position = "none", 
        #axis.text.x = element_text(angle=90, hjust=1, size = 8)
  )
grid.arrange(p1, p2, p3, ncol=3)
