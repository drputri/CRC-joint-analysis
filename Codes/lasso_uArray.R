#####################################################
## STATISTICAL ANALYSIS                            ##
## LASSO: uArray                                   ##
## Prepared by: Dea Putri (dputri@its.jnj.com)     ##
#####################################################

rna <- get(load("Data/rna_crcad.Rda"))
load("Data/group_crcad.Rda")
source("Codes/functions.R")


## Load libraries
library(parallel)
library(ggplot2)
library(glmnet)
library(magrittr)
library(dplyr)
library(tidyr)






#----------------------------------------------------------------------------------------
# FEATURES SELECTION
#----------------------------------------------------------------------------------------
# Choose the MCCV type
MCCV.type <- 2


# Create a cluster of cores
cl <- makeCluster(getOption("cl.cores", 4))
clusterExport(cl=cl, 
              varlist=c("rna", "group", "select_feature", "MCCV.type"))


# Main algorithm
system.time({
  out <- parLapply(cl=cl,
                   X=1:500,
                   fun=function(i) {
                     require(glmnet)
                     
                     
                     # Set seed for the sample randomization
                     set.seed(i)
                     
                     
                     # Samples on the different MCCV type
                     if(MCCV.type==0) {
                       
                       # sample 2/3 of the subjects at random
                       id_keep <- sample(x=rownames(rna),
                                         size=floor(2*nrow(rna)/3),
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
                       id_keep <- c(sample(x=ad, size=11, replace=FALSE),
                                    sample(x=crc, size=11, replace=FALSE))
                       
                     }
                     
                     
                     # subset data
                     x <- rna[rownames(rna) %in% id_keep, ]
                     
                     # x <- rna #--old version
                     
                     # LASSO glmnet
                     # Set seed for lasso model
                     feat_select <- try(select_feature(X=x, y=group[names(group) %in% id_keep], 
                                                       fold=3, seed=i+1200), silent = TRUE)
                     # feat_select <- select_feature(X=x, y=group, fold=5, seed=i) #--old version
                     
                     return(feat_select)
                   }
  )
})


# Stop the cluster
stopCluster(cl)
gc() # Running time: 125.270s

tmp <- out


# Gather the selected feature
load("Data/esetRna_new.Rda") 

feat_select <- do.call("c", out) %>%
  table %>%
  data.frame() %>%
  set_colnames(c("Feat", "Freq")) 

fdata <- as(featureData(esetRna), "data.frame")
temp <- subset(fdata, select=c(SYMBOL))

feat_select <- merge(feat_select, temp, by.x="Feat", by.y="row.names", all.x = TRUE) %>%
  arrange(desc(Freq)) %>%
  mutate(Ind=ifelse(Freq>=250, "1", "0"))

save(feat_select, file="Data/uArray_selected feat_MCCVType2.Rda")


# Plot
ggplot(data=feat_select[1:25, ], aes(x=reorder(SYMBOL, -Freq), y=Freq)) +
  geom_col(aes(fill=Ind)) +
  scale_fill_manual(name="> 500",
                    labels=c("No", "Yes"),
                    values=c("1" = "tomato3", "0" = "grey54")) +
  labs(x="Features", y="Frequency", title="Number of repeats: 500") +
  lims(y=c(0, 500)) +
  theme(panel.background=element_blank(),
        panel.border=element_blank(),
        legend.position="none", 
        axis.text.x=element_text(angle=90, hjust=1, size = 8))






#----------------------------------------------------------------------------------------
# FIXED GENES: 3-fold cross validation + performance assessment
#----------------------------------------------------------------------------------------
## Load libraries
library(parallel)
library(ggplot2)
library(glmnet)
library(magrittr)
library(dplyr)
library(tidyr)


# Load data set and functions
rna <- get(load("Data/rna_crcad.Rda"))
load("Data/group_crcad.Rda")
source("Codes/functions.R")


# Setup the selected feature data
feat0 <- get(load("Data/uArray_selected feat_MCCVType0.Rda"))
feat1 <- get(load("Data/uArray_selected feat_MCCVType1.Rda"))
feat2 <- get(load("Data/uArray_selected feat_MCCVType2.Rda"))
feat3 <- get(load("Data/uArray_selected feat_MCCVType3.Rda"))


# Find the most common uarray features
feat0$rank <- seq(1:nrow(feat0))
feat1$rank <- seq(1:nrow(feat1))
feat2$rank <- seq(1:nrow(feat2))
feat3$rank <- seq(1:nrow(feat3))

feat.all <- merge(feat0, feat1, by = c("Feat", "SYMBOL")) %>%
  rename("Freq.feat0" = Freq.x,
         "Ind.feat0" = Ind.x,
         "rank.feat0" = rank.x,
         "Freq.feat1" = Freq.y,
         "Ind.feat1" = Ind.y,
         "rank.feat1" = rank.y,)
feat.all <- merge(feat.all, feat2, by = c("Feat", "SYMBOL")) %>%
  rename("Freq.feat2" = Freq,
         "Ind.feat2" = Ind,
         "rank.feat2" = rank)
feat.all <- merge(feat.all, feat3, by = c("Feat", "SYMBOL")) %>%
  rename("Freq.feat3" = Freq,
         "Ind.feat3" = Ind,
         "rank.feat3" = rank) %>%
  mutate(sum.rank = rank.feat0 + rank.feat1 + rank.feat2 + rank.feat3) %>%
  arrange(sum.rank)


# Create a cluster of cores
cl <- makeCluster(getOption("cl.cores", 4))
clusterExport(cl=cl, 
              varlist=c("rna", "group", 
                        "feat.all", "perf_eval"))


# Main algorithm
system.time({
  out <- parLapply(cl = cl,
                   X = 1:500,
                   fun = function(i) {
                     
                     require(glmnet)
                     require(limma)
                     require(magrittr)
                     require(mixOmics)
                     require(caret)
                     require(dplyr)
                     require(magrittr)
                     
                     
                     # Performance measure
                     # Change the MCCV.type to obtain the different scenario
                     set.seed(i)
                     model.perf <- perf_eval(k = 10, X = rna, y = group, 
                                             df.feat = feat.all, rank.by = "sum.rank",
                                             MCCV.type = 2, fold = 3)
                     
                     return(model.perf)
                     
                   }
  )
})


# Stop the cluster
stopCluster(cl)
gc() # Running time: 33.270s

# Change the file name based on the MCCV type run
save(out, file = "Data/uArray_model performance_MCCV2.Rda")



## Calculate the data summary
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(magrittr)


# Original MCCV
load("Data/uArray_model performance_MCCV0.Rda")

tmp <- lapply(out, function(a) { #out[[1]]
  x <- lapply(a, function(b) { #out[[1]]$k
    as.data.frame(b$evaluation)
  })
  
  y <- do.call("rbind", x) %>%
    rownames_to_column(var="k")
})

eval <- do.call("rbind", tmp) %>%
  rownames_to_column(var="i") %>%
  pivot_longer(!c(k, i), names_to="stat", values_to="value") %>%
  group_by(k, stat) %>%
  summarise(avg = mean(value, na.rm = TRUE),
            se = sd(value, na.rm = TRUE)/sqrt(500)) %>%
  mutate(ymin = avg - 1.96*se,
         ymax = avg + 1.96*se) %>%
  mutate(ymin = ifelse(ymin < 0, 0, ymin),
         ymax = ifelse(ymax > 1, 1, ymax),
         MCCV.type= "0",
         k = factor(k, levels = c(2:10)))


# MCCV type 1
load("Data/uArray_model performance_MCCV1.Rda")

tmp <- lapply(out, function(a) { #out[[1]]
  x <- lapply(a, function(b) { #out[[1]]$k
    as.data.frame(b$evaluation)
  })
  
  y <- do.call("rbind", x) %>%
    rownames_to_column(var="k")
})

tmp <- do.call("rbind", tmp) %>%
  rownames_to_column(var="i") %>%
  pivot_longer(!c(k, i), names_to="stat", values_to="value") %>%
  group_by(k, stat) %>%
  summarise(avg = mean(value, na.rm = TRUE),
            se = sd(value, na.rm = TRUE)/sqrt(500)) %>%
  mutate(ymin = avg - 1.96*se,
         ymax = avg + 1.96*se) %>%
  mutate(ymin = ifelse(ymin < 0, 0, ymin),
         ymax = ifelse(ymax > 1, 1, ymax),
         MCCV.type= "1",
         k = factor(k, levels = c(2:10)))

eval <- rbind(eval, tmp)


# MCCV type 2
load("Data/uArray_model performance_MCCV2.Rda")

tmp <- lapply(out, function(a) { #out[[1]]
  x <- lapply(a, function(b) { #out[[1]]$k
    as.data.frame(b$evaluation)
  })
  
  y <- do.call("rbind", x) %>%
    rownames_to_column(var="k")
})

tmp <- do.call("rbind", tmp) %>%
  rownames_to_column(var="i") %>%
  pivot_longer(!c(k, i), names_to="stat", values_to="value") %>%
  group_by(k, stat) %>%
  summarise(avg = mean(value, na.rm = TRUE),
            se = sd(value, na.rm = TRUE)/sqrt(500)) %>%
  mutate(ymin = avg - 1.96*se,
         ymax = avg + 1.96*se) %>%
  mutate(ymin = ifelse(ymin < 0, 0, ymin),
         ymax = ifelse(ymax > 1, 1, ymax),
         MCCV.type= "2",
         k = factor(k, levels = c(2:10)))

eval <- rbind(eval, tmp)


# MCCV type 3
load("Data/uArray_model performance_MCCV3.Rda")

tmp <- lapply(out, function(a) { #out[[1]]
  x <- lapply(a, function(b) { #out[[1]]$k
    as.data.frame(b$evaluation)
  })
  
  y <- do.call("rbind", x) %>%
    rownames_to_column(var="k")
})

tmp <- do.call("rbind", tmp) %>%
  rownames_to_column(var="i") %>%
  pivot_longer(!c(k, i), names_to="stat", values_to="value") %>%
  group_by(k, stat) %>%
  summarise(avg = mean(value, na.rm = TRUE),
            se = sd(value, na.rm = TRUE)/sqrt(500)) %>%
  mutate(ymin = avg - 1.96*se,
         ymax = avg + 1.96*se) %>%
  mutate(ymin = ifelse(ymin < 0, 0, ymin),
         ymax = ifelse(ymax > 1, 1, ymax),
         MCCV.type= "3",
         k = factor(k, levels = c(2:10)))

eval <- rbind(eval, tmp)

save(eval, file = "Data/lassoUArray_eval_differentMCCV.Rda")

  
# Plot
pd <- position_dodge(0.1)
p0 <- ggplot(eval[eval$stat=="BER", ], aes(x=k, y=avg, color=MCCV.type, group=MCCV.type)) + 
  #geom_jitter(data = eval[, c("BER", "k")], aes(x = k, y = BER), color = "lightblue", alpha = 0.5, width = 0.1)
  geom_line(position=pd) +
  geom_point(position=pd) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax), width=.1, position=pd) +
  # ylim(0.8, 1) +
  labs(y = "BER") +
  theme(panel.background = element_blank(),
        panel.border = element_blank()#,
        #legend.position = "none", 
        #axis.text.x = element_text(angle=90, hjust=1, size = 8)
  )


# Other performance metrices (MCE, SEN, SPE)
library(gridExtra)
p1 <- ggplot(eval[eval$stat=="mce", ], aes(x=k, y=avg, color=MCCV.type, group=MCCV.type)) + 
  #geom_jitter(data = eval[, c("BER", "k")], aes(x = k, y = BER), color = "lightblue", alpha = 0.5, width = 0.1)
  geom_line(position=pd) +
  geom_point(position=pd) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax), width=.1, position=pd) +
  # ylim(0.8, 1) +
  labs(y = "MCE") +
  theme(panel.background = element_blank(),
        panel.border = element_blank()#,
        #legend.position = "none", 
        #axis.text.x = element_text(angle=90, hjust=1, size = 8)
  )

p2 <- ggplot(eval[eval$stat=="sensitivity", ], aes(x=k, y=avg, color=MCCV.type, group=MCCV.type)) + 
  #geom_jitter(data = eval[, c("BER", "k")], aes(x = k, y = BER), color = "lightblue", alpha = 0.5, width = 0.1)
  geom_line(position=pd) +
  geom_point(position=pd) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax), width=.1, position=pd) +
  # ylim(0.8, 1) +
  labs(y = "Sensitivity") +
  theme(panel.background = element_blank(),
        panel.border = element_blank()#,
        #legend.position = "none", 
        #axis.text.x = element_text(angle=90, hjust=1, size = 8)
  )

p3 <- ggplot(eval[eval$stat=="specificity", ], aes(x=k, y=avg, color=MCCV.type, group=MCCV.type)) + 
  #geom_jitter(data = eval[, c("BER", "k")], aes(x = k, y = BER), color = "lightblue", alpha = 0.5, width = 0.1)
  geom_line(position=pd) +
  geom_point(position=pd) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax), width=.1, position=pd) +
  # ylim(0.8, 1) +
  labs(y = "Specificity") +
  theme(panel.background = element_blank(),
        panel.border = element_blank()#,
        #legend.position = "none", 
        #axis.text.x = element_text(angle=90, hjust=1, size = 8)
  )

p4 <- ggplot(eval[eval$stat=="npv", ], aes(x=k, y=avg, color=MCCV.type, group=MCCV.type)) + 
  #geom_jitter(data = eval[, c("BER", "k")], aes(x = k, y = BER), color = "lightblue", alpha = 0.5, width = 0.1)
  geom_line(position=pd) +
  geom_point(position=pd) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax), width=.1, position=pd) +
  # ylim(0.8, 1) +
  labs(y = "NPV") +
  theme(panel.background = element_blank(),
        panel.border = element_blank()#,
        #legend.position = "none", 
        #axis.text.x = element_text(angle=90, hjust=1, size = 8)
  )

p5 <- ggplot(eval[eval$stat=="ppv", ], aes(x=k, y=avg, color=MCCV.type, group=MCCV.type)) + 
  #geom_jitter(data = eval[, c("BER", "k")], aes(x = k, y = BER), color = "lightblue", alpha = 0.5, width = 0.1)
  geom_line(position=pd) +
  geom_point(position=pd) +
  geom_errorbar(aes(ymin=ymin, ymax=ymax), width=.1, position=pd) +
  # ylim(0.8, 1) +
  labs(y = "PPV") +
  theme(panel.background = element_blank(),
        panel.border = element_blank()#,
        #legend.position = "none", 
        #axis.text.x = element_text(angle=90, hjust=1, size = 8)
  )

grid.arrange(p0, p1, p2, p3, p4, p5, ncol=3)
