#####################################################
## STATISTICAL ANALYSIS                            ##
## LASSO joint: uArray all MCCV types              ##
## Prepared by: Dea Putri (dputri@its.jnj.com)     ##
#####################################################



## Load data
out0 <- get(load("Data/jointLasso_uArray_MCCVType0.Rda"))
out1 <- get(load("Data/jointLasso_uArray_MCCVType1.Rda"))
out2 <- get(load("Data/jointLasso_uArray_MCCVType2.Rda"))
out3 <- get(load("Data/jointLasso_uArray_MCCVType3.Rda"))

out <- list(out0, out1, out2, out3)
names(out) <- c("0", "1", "2", "3")



## Calculate the data summary
library(tibble)
library(tidyr)
library(magrittr)
library(dplyr)

type <- c("0", "1", "2", "3")
perf <- lapply(type, function(t){
  
  o <- out[[t]]
  
  tmp <- lapply(o, function(x){
    
    as.data.frame(x$evaluation)
    
  })
  tmp <- do.call("rbind", tmp) %>%
    rownames_to_column(var = "i") %>%
    mutate(MCCVType = t)
  
})

perf <- do.call("rbind", perf)

save(perf, file="Data/jointLasso_uArray_perf.Rda")


# Density plot
perf <- perf %>%
  pivot_longer(cols = 2:8,
               names_to = "var",
               values_to = "value")

perf %>%
  filter(var != "kappa") %>%
  ggplot(aes(x = value, color = MCCVType)) +
  geom_density(adjust = 2) +
  facet_wrap(~ var) +
  scale_color_hue(labels = c("0" = "Test: 14-7",
                             "1" = "Test: 11-11",
                             "2" = "Test: 20-10",
                             "3" = "Test: 30-10")) +
  theme_bw()



## Lasso coefficients
coef <- lapply(type, function(t){
  
  o <- out[[t]]
  
  tmp <- lapply(o, function(x){
    
    as.data.frame(x$coefficients)
    
  })
  tmp <- do.call("rbind", tmp) %>%
    mutate(MCCVType = t)
  
})

coef <- do.call("rbind", coef) %>%
  filter(feat != "(Intercept)")
rownames(coef) <- NULL

save(coef, file="Data/jointLasso_uArray_coef.Rda")


# Get the information of the microbiome
load("Data/uBiome_taxtable.Rda")
coef <- merge(coef, taxTable, by.x = "feat", by.y = "row.names", all.x = TRUE) %>%
  filter(feat %in% grep("SVs", coef$feat, ignore.case = TRUE, value = TRUE)) 
sum.coef <- coef %>%
  group_by(MCCVType, Genus) %>%
  summarise(n = n()) %>%
  arrange(MCCVType, -n) %>%
  filter(n >= 160)

sum.coef %>%
  ggplot(aes(x = Genus, y = n, fill = MCCVType)) +
  geom_bar(stat = "identity",position = position_dodge(width = 0.5),
           alpha = .6) +
  scale_fill_hue(labels = c("0" = "Test: 14-7",
                             "1" = "Test: 11-11",
                             "2" = "Test: 20-10",
                             "3" = "Test: 30-10")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) 
  
