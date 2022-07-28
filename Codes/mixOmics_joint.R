## STATISTICAL ANALYSIS
## JOINT ANALYSIS
## Prepared by: Dea Putri (dputri@its.jnj.com)
#####################################################
setwd("/mnt/exports/shared/home/dputri/CRC microbiome/drp_crcmicrobiome/Joint analysis")






#----------------------------------------------------------------------------------------
# PACKAGES
#----------------------------------------------------------------------------------------
library(mixOmics)







#----------------------------------------------------------------------------------------
# DATASET PREPARATION
#----------------------------------------------------------------------------------------
load("Data/rna_crcad.Rda")         # expression matrix
load("Data/otu_431_clr.Rda")       # otu relative abundance matrix
load("Data/group_crcad.Rda")


# Data for mixOmics
data <- list(uArray = rna.crcad,
             uBiome = otu.clr)
save(data, file = "Data/mixomics_data.Rda")







#----------------------------------------------------------------------------------------
# MIXOMICS
#----------------------------------------------------------------------------------------
# Design matrix -- full weighted design matrix [USE THIS]
design = matrix(0.1, ncol = length(data)+1, nrow = length(data)+1, 
                dimnames = list(c(names(data), "group"), c(names(data), "group")))
design[, colnames(design) == "group"] <- 1
design[rownames(design) == "group",  ] <- 1
diag(design) = 0


# Design matrix -- Null design
design_null = matrix(0, ncol = length(data)+1, nrow = length(data)+1, 
                     dimnames = list(c(names(data), "group"), c(names(data), "group")))
design_null[, colnames(design_null) == "group"] <- 1
design_null[rownames(design_null) == "group",  ] <- 1



## Full weighted design
sgccda.res <- block.splsda(X=data, Y=group, ncomp=5, 
                           design=design)
set.seed(1)
t1 = proc.time()
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds=5, nrepeat=500)
t2 = proc.time()
running_time = t2 - t1; running_time   # 19.6 minutes

save(perf.diablo, file="Data/performance diablo.Rda")

source("Codes/funct_plotperf.R")
plot.perf.sgccda.mthd(perf.diablo, sd = TRUE)

ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "mahalanobis.dist"]


# Tuning keepX: tuning the number of variables to select in each data
set.seed(1)
test.keepX=list(uArray=seq(10, 100, 10),
                uBiome=seq(10, 20, 1))
BPPARAM <- BiocParallel::MulticoreParam(workers = parallel::detectCores()-1)
tune.data <- tune.block.splsda(X=data, Y=group, ncomp=2, 
                              test.keepX=test.keepX, design=design,
                              folds=5, nrepeat=500,
                              BPPARAM=BPPARAM)
list.keepX = tune.data$choice.keepX


# Final model
sgccda.res = block.splsda(X=data, Y=group, ncomp=2, 
                          keepX=list.keepX, design=design)
save(sgccda.res, file = "Output/mixomics_finalmodel.Rda")

selectVar(sgccda.res, block='uArray', comp=1)$uArray$name           # 15 genes
selectVar(sgccda.res, block='uBiome', comp=1)$uBiome$name           # 10 genus
plotDiablo(sgccda.res, ncomp=1)
plotIndiv(sgccda.res, ind.names = TRUE, legend = TRUE, 
          title='DIABLO')

png("mixomics_full_arrow.png")
plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')
dev.off()

plotVar(sgccda.res, var.names = FALSE, style = 'graphics', legend = TRUE, 
        pch = c(16, 17), cex = c(2,2), col = c('brown1', 'lightgreen'))
circosPlot(sgccda.res, cutoff = 0.7, line = TRUE, 
           color.blocks= c('brown1', 'lightgreen'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)
network(sgccda.res, blocks = c(1,2),
        color.node = c('brown1', 'lightgreen'), cutoff = 0.4)
plotLoadings(sgccda.res, comp = 2, contrib = 'max', method = 'median')







#----------------------------------------------------------------------------------------
# GET THE GENUS AND GENES NAME
#----------------------------------------------------------------------------------------
load("Output/mixomics_finalmodel.Rda")

library(Biobase)
library(magrittr)
library(tidyr)
library(dplyr)

#----- FUNCTIONS TO GET THE LOADING VECTORS -----#
get.loadings.ndisplay <- function(object,
                                  comp,
                                  block,
                                  name.var,
                                  name.var.complete,
                                  ndisplay)
{
  ##selectvar
  selected.var = selectVar(object, comp = comp, block = block) # gives name and values of the blocks in 'block'
  name.selected.var = selected.var[[1]]$name
  value.selected.var = selected.var[[1]]$value
  
  # ndisplay
  # ------
  # if null set by default to all variables from selectVar
  if (is.null(ndisplay))
  {
    ndisplay.temp = length(name.selected.var)
  } else if (ndisplay > length(name.selected.var)) {
    message("'ndisplay' value is larger than the number of selected variables! It has been reseted to ", length(name.selected.var), " for block ", block)
    ndisplay.temp = length(name.selected.var)
  } else {
    ndisplay.temp = ndisplay
  }
  
  name.selected.var = name.selected.var[1:ndisplay.temp]
  value.selected.var = value.selected.var[1:ndisplay.temp,]
  
  #comp
  # ----
  if (is(object, c("mixo_pls","mixo_spls", "rcc")))# cause pls methods just have 1 ncomp, block approaches have different ncomp per block
  {
    ncomp = object$ncomp
    object$X = list(X = object$X, Y = object$Y) # so that the data is in object$X, either it's a pls or block approach
  } else {
    ncomp = object$ncomp[block]
  }
  
  if (any(max(comp) > ncomp))
    stop(paste("Argument 'comp' should be less or equal to ", ncomp))
  
  names.block = as.character(names(selected.var)[1]) #it should be one block and ncomp, so we take the first one
  
  X = object$X[names.block][[1]]
  
  #name.var
  ind.match = match(name.selected.var, colnames(X)) # look at the position of the selected variables in the original data X
  if(!is.null(name.var))
  {
    if(length(name.var)!= ncol(X))
      stop("For block '", names.block,"', 'name.var' should be a vector of length ", ncol(X))
    
    colnames.X = as.character(name.var[ind.match]) # get the
  }else{
    colnames.X = as.character(colnames(X))[ind.match]
  }
  X = X[, name.selected.var, drop = FALSE] #reduce the problem to ndisplay
  
  #completing colnames.X by the original names of the variables when missing
  if (name.var.complete == TRUE)
  {
    ind = which(colnames.X == "")
    if (length(ind) > 0)
      colnames.X[ind] = colnames(X)[ind]
  }
  
  
  return(list(X = X, names.block = names.block, colnames.X = colnames.X, name.selected.var = name.selected.var, value.selected.var = value.selected.var))
}

#-----------------------------------------------------#





## Microarray features
res = get.loadings.ndisplay(object=sgccda.res, comp=1, block='uArray', 
                            name.var=NULL, name.var.complete=FALSE, ndisplay=NULL)
load1.uArray <- data.frame(var=res$name.selected.var, 
                           load.value=res$value.selected.var) %>%
  mutate(abs.load.val=abs(load.value))


# Get microarray features from LASSO
feat <- get(load("Data/uArray_selected feat.Rda"))
get.feat <- feat[1:6, ]

load1.uArray <- merge(load1.uArray, feat, by.x="var", by.y="Feat", all.x=TRUE)
load1.uArray$prop <- load1.uArray$Freq/1000
load1.uArray$Ind <- ifelse(load1.uArray$prop > 0.4, 1, 0)


# Get the genes name
load("Data/esetRna_new.Rda") 

fdata <- as(featureData(esetRna), "data.frame")
temp <- subset(fdata, select = c(SYMBOL))

load1.uArray <- merge(load1.uArray, temp, by.x="var", by.y="row.names", all.x=TRUE) %>%
  rename(name="SYMBOL.y") %>%
  select(-c(SYMBOL.x)) %>%
  arrange(desc(abs.load.val))


# Plot
ggplot(load1.uArray[1:40, ], aes(x=abs.load.val, y=prop, colour=droplevels(as.factor(Ind)), label=name))+
  geom_text(hjust=0, vjust=0) +
  scale_color_manual(na.translate=FALSE,
                     values=c("black", "red")) +
  labs(color="Proportion > 0.4",
       x="Abs loading weight",
       y="Proportion of selection") +
  lims(x=c(0.05, 0.3)) +
  theme_bw()



## Genera
res = get.loadings.ndisplay(object=sgccda.res, comp=1, block='uBiome', 
                            name.var=NULL, name.var.complete=FALSE, ndisplay=NULL)
load1.uBiome <- data.frame(var=res$name.selected.var, 
                           load.value=res$value.selected.var) %>%
  mutate(abs.load.val=abs(load.value))


res = get.loadings.ndisplay(object=sgccda.res, comp=2, block='uBiome', 
                            name.var=NULL, name.var.complete=FALSE, ndisplay=NULL)
load2.uBiome <- data.frame(var=res$name.selected.var, 
                           load.value=res$value.selected.var) %>%
  mutate(abs.load.val=abs(load.value))

load.uBiome <- rbind(load1.uBiome, load2.uBiome)

# Get microarray features from LASSO
feat <- get(load("Data/uBiome_selected feat.Rda"))

load.uBiome <- merge(load.uBiome, feat, by.x="var", by.y="Feat", all.x=TRUE)
load.uBiome$prop <- load.uBiome$Freq/1000
load.uBiome$Ind <- ifelse(load.uBiome$prop > 0.4, 1, 0) 
load.uBiome <- load.uBiome[order(-load.uBiome$abs.load.val), ]



# Plot
ggplot(load.uBiome, aes(x=abs.load.val, y=reorder(Genus, abs.load.val))) +
  geom_point(stat="identity", aes(col=as.factor(Ind), size=prop)) +
  scale_color_manual(na.translate=FALSE,
                     values=c("black", "red")) +
  labs(color="Lasso",
       size="Prop. of selection",
       x="Abs loading weight",
       y="Genera") +
  theme_bw()
