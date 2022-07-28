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







#----------------------------------------------------------------------------------------
# FEATURES SELECTION
#----------------------------------------------------------------------------------------
load("Data/otu_431_clr.Rda")
load("Data/group_crcad.Rda")
source("Codes/functions.R")


## Load libraries
library(parallel)
library(ggplot2)
library(glmnet)
library(magrittr)
library(dplyr)


# Choose the MCCV type
MCCV.type <- 3


# Create a cluster of cores
cl <- makeCluster(getOption("cl.cores", 4))
clusterExport(cl=cl, 
              varlist=c("otu.clr", "group", "select_feature", "MCCV.type"))


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
                       id_keep <- sample(x=rownames(otu.clr),
                                         size=floor(2*nrow(otu.clr)/3),
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
                     
                     
                     # Subset data
                     x <- otu.clr[rownames(otu.clr) %in% id_keep, ]
                     
                     # x <- otu.clr #--old version
                     
                     # LASSO glmnet
                     feat_select <- try(select_feature(X=x, y=group[names(group) %in% id_keep], 
                                                       fold=3, seed=i+1200), silent = TRUE)
                     # feat_select <- select_feature(X=x, y=group, fold=5, seed=i) #--old version
                     
                     return(feat_select)
                   }
  )
})


# Stop the cluster
stopCluster(cl)
gc() # Running time: 14.615


# Gather the selected feature
load("Data/uBiome_taxtable.Rda") 

feat_select <- do.call("c", out) %>%
  table %>%
  data.frame() %>%
  set_colnames(c("Feat", "Freq")) 

feat_select <- merge(feat_select, taxTable, by.x="Feat", by.y="row.names", all.x = TRUE) %>%
  arrange(desc(Freq)) %>%
  mutate(Ind=ifelse(Freq>=250, "1", "0"))

save(feat_select, file="Data/uBiome_selected feat_MCCVType3.Rda")


# Plot
ggplot(data=feat_select[1:25, ], aes(x=reorder(Genus, -Freq), y=Freq)) +
  geom_col(aes(fill=Ind)) +
  scale_fill_manual(name="> 500",
                    labels=c("No", "Yes"),
                    values=c("1" = "tomato3", "0" = "grey54")) +
  labs(x="Features", y="Frequency", title="Number of repeats: 1,000") +
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
load("Data/otu_431_clr.Rda")
load("Data/group_crcad.Rda")
source("Codes/functions.R")


# Setup the selected feature data
gene0 <- get(load("Data/uBiome_selected feat_MCCVType0.Rda"))
gene1 <- get(load("Data/uBiome_selected feat_MCCVType1.Rda"))
gene2 <- get(load("Data/uBiome_selected feat_MCCVType2.Rda"))
gene3 <- get(load("Data/uBiome_selected feat_MCCVType3.Rda"))


# Find the most common ubiome genes
gene0$rank <- seq(1:nrow(gene0))
gene1$rank <- seq(1:nrow(gene1))
gene2$rank <- seq(1:nrow(gene2))
gene3$rank <- seq(1:nrow(gene3))

gene.all <- merge(gene0, gene1, by = c("Feat", "Genus")) %>%
  rename("Freq.gene0" = Freq.x,
         "Ind.gene0" = Ind.x,
         "rank.gene0" = rank.x,
         "Freq.gene1" = Freq.y,
         "Ind.gene1" = Ind.y,
         "rank.gene1" = rank.y,)
gene.all <- merge(gene.all, gene2, by = c("Feat", "Genus")) %>%
  rename("Freq.gene2" = Freq,
         "Ind.gene2" = Ind,
         "rank.gene2" = rank)
gene.all <- merge(gene.all, gene3, by = c("Feat", "Genus")) %>%
  rename("Freq.gene3" = Freq,
         "Ind.gene3" = Ind,
         "rank.gene3" = rank) %>%
  mutate(sum.rank = rank.gene0 + rank.gene1 + rank.gene2 + rank.gene3) %>%
  arrange(sum.rank)


# Create a cluster of cores
cl <- makeCluster(getOption("cl.cores", 4))
clusterExport(cl=cl, 
              varlist=c("otu.clr", "group", 
                        "gene.all", "perf_eval"))


# Main algorithm
system.time({
  out <- parLapply(cl=cl,
                   X=1:500,
                   fun=function(i) {
                     require(glmnet)
                     require(limma)
                     require(magrittr)
                     require(mixOmics)
                     require(caret)
                     require(dplyr)
                     require(magrittr)
                     
                     
                     # Performance measure
                     set.seed(i)
                     model.perf <- perf_eval(k = 8, X = otu.clr, y = group, 
                                             df.feat = gene.all, rank.by = "sum.rank",
                                             MCCV.type = 3, fold = 3)
                     
                     return(model.perf)
                   }
  )
})


# Stop the cluster
stopCluster(cl)
gc() # Running time: 2744.876s

save(out, file="Data/uBiome_model performance_MCCV3.Rda")


## Calculate the data summary
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(magrittr)


# Original MCCV
load("Data/uBiome_model performance_MCCV0.Rda")

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
         k = factor(k, levels = c(2:8)))


# MCCV type 1
load("Data/uBiome_model performance_MCCV1.Rda")

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
         k = factor(k, levels = c(2:8)))

eval <- rbind(eval, tmp)


# MCCV type 2
load("Data/uBiome_model performance_MCCV2.Rda")

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
         k = factor(k, levels = c(2:8)))

eval <- rbind(eval, tmp)


# MCCV type 3
load("Data/uBiome_model performance_MCCV3.Rda")

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
         k = factor(k, levels = c(2:8)))

eval <- rbind(eval, tmp)

save(eval, file = "Data/lassoBiome_eval_differentMCCV.Rda")


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
