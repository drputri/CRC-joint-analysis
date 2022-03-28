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






#----------------------------------------------------------------------------------------
# FEATURES SELECTION
#----------------------------------------------------------------------------------------
# Create a cluster of cores
cl <- makeCluster(getOption("cl.cores", 4))
clusterExport(cl=cl, 
              varlist=c("rna", "group", "select_feature"))


# Main algorithm
system.time({
  out <- parLapply(cl=cl,
                   X=1:1000,
                   fun=function(i) {
                     require(glmnet)
                     
                     # sample 2/3 of the subjects at random
                     id_keep <- sample(x=rownames(rna),
                                       size=floor(2*nrow(rna)/3),
                                       replace=FALSE)
                     
                     # subset data
                     x <- rna[rownames(rna) %in% id_keep, ]
                     
                     # x <- rna #--old version
                     
                     # LASSO glmnet
                     feat_select <- select_feature(X=x, y=group[names(group) %in% id_keep], fold=5, seed=i)
                     # feat_select <- select_feature(X=x, y=group, fold=5, seed=i) #--old version
                     
                     return(feat_select)
                   }
  )
})


# Stop the cluster
stopCluster(cl)
gc() # Running time: 125.270s


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
  mutate(Ind=ifelse(Freq>=500, "1", "0"))

save(feat_select, file="Data/uArray_selected feat.Rda")


# Plot
ggplot(data=feat_select[1:25, ], aes(x=reorder(SYMBOL, -Freq), y=Freq)) +
  geom_col(aes(fill=Ind)) +
  scale_fill_manual(name="> 500",
                    labels=c("No", "Yes"),
                    values=c("1" = "tomato3", "0" = "grey54")) +
  labs(x="Features", y="Frequency", title="Number of repeats: 1,000") +
  lims(y=c(0, 1000)) +
  theme(panel.background=element_blank(),
        panel.border=element_blank(),
        legend.position="none", 
        axis.text.x=element_text(angle=90, hjust=1, size = 8))






#----------------------------------------------------------------------------------------
# FIXED GENES: 5-fold cross validation + performance assessment
#----------------------------------------------------------------------------------------
load("Data/uArray_selected feat.Rda")


# Create a cluster of cores
cl <- makeCluster(getOption("cl.cores", 4))
clusterExport(cl=cl, 
              varlist=c("rna", "group", 
                        "feat_select", "perf_eval"))


# Main algorithm
system.time({
  out <- parLapply(cl=cl,
                   X=1:1000,
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
                     model.perf <- perf_eval(k=8, X=rna, y=group, df.feat=feat_select, 
                                             fold=5, alpha=0)
                     
                     return(model.perf)
                   }
  )
})


# Stop the cluster
stopCluster(cl)
gc() # Running time: 125.270s

save(out, file="Data/uArray_model performance.Rda")


## Calculate the data summary
library(tibble)
library(tidyr)
load("Data/uArray_model performance.Rda")

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
  summarise(avg=mean(value),
            std=sd(value))
  
# Plot
pd <- position_dodge(0.1)
ggplot(eval[eval$stat=="BER", ], aes(x=k, y=avg)) + 
  geom_errorbar(aes(ymin=avg-std, ymax=avg+std), width=.1, position=pd) +
  #geom_jitter(data = eval[, c("BER", "k")], aes(x = k, y = BER), color = "lightblue", alpha = 0.5, width = 0.1)
  geom_line(position=pd) +
  geom_point(position=pd) +
  # ylim(0, .2) +
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
