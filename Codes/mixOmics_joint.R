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
load("uBiome_tf20.Rda")       # top 20 genera of uBiome from univariate analysis
load("uArray_tf15.Rda")       # top 15 features of uArray from univariate analysis
load("rna_crcad.Rda")         # expression matrix
load("oturl_crcad.Rda")       # otu relative abundance matrix



## Microbiome
load("Data/Obj_Genus.rda")
metadata <- data.frame(sample_data(Obj_Genus))
otu <- data.frame(otu_table(Obj_Genus))



## Microarray 
load("Data/esetRna_new.Rda") 
rna <- exprs(esetRna)
pdata <- pData(esetRna)


# Subject selection
subj.rna <- unique(pdata$ID[which(pdata$biopt %in% c("6 : Polyp", "7 : Tumor"))])
subj.mb <- as.numeric(as.character(unique(metadata$Patient_ID[which(metadata$Location_type_colon %in% c("Adenoma", "Carcinoma"))])))
subject <- as.numeric(intersect(subj.rna, subj.mb))


# Subset uArray data
ftr <- as.character(tf_k15$Var1[which(tf_k15$select == "yes")])
rna.crcad <- as.data.frame.matrix(rna.crcad)
rna.crcad <- rna.crcad[, colnames(rna.crcad) %in% ftr]

subject <- as.character(subject)
rna.crcad <- rna.crcad[rownames(rna.crcad) %in% subject, ]
rna.crcad <- rna.crcad[order(row.names(rna.crcad)), ]  
rna.crcad <- as.matrix.data.frame(rna.crcad)


# Subset uBiome data
genus <- as.character(tf_k20$Var1[which(tf_k20$select == "yes")])
metadata$sampleID <- rownames(metadata)
var <- c("sampleID", "Patient_ID")
temp <- metadata[var]
temp <- temp[which(temp$Patient_ID %in% subject), ]

oturl.crcad <- as.data.frame.matrix(oturl.crcad)
oturl.crcad$sampleID <- rownames(oturl.crcad)
oturl.crcad <- merge(oturl.crcad, temp, by = "sampleID")
rownames(oturl.crcad) <- oturl.crcad$Patient_ID
oturl.crcad <- subset(oturl.crcad, select = -c(sampleID, Patient_ID))
oturl.crcad <- oturl.crcad[, colnames(oturl.crcad) %in% genus]
oturl.crcad <- oturl.crcad[rownames(oturl.crcad) %in% subject, ]
oturl.crcad <- oturl.crcad[order(row.names(oturl.crcad)), ] 
oturl.crcad <- as.matrix.data.frame(oturl.crcad)


# Outcome -- Y
var <- c("Patient_ID", "GroupAssignment")
temp <- metadata[var]
temp <- unique(temp[which(temp$Patient_ID %in% rownames(oturl.crcad)), ])
temp <- temp[order(temp$Patient_ID), ]
group <- factor(temp$GroupAssignment)


# Data for mixOmics
data <- list(uArray = rna.crcad,
             uBiome = oturl.crcad)
save(data, file = "Data/mixomics_data.Rda")







#----------------------------------------------------------------------------------------
# MIXOMICS
#----------------------------------------------------------------------------------------
# Design matrix -- full weighted design matrix
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



## FULL DESIGN
sgccda.res <- block.splsda(X = data, Y = group, ncomp = 10, 
                           design = design)
set.seed(1)
t1 = proc.time()
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 3, nrepeat = 1000)
t2 = proc.time()
running_time = t2 - t1; running_time   # 38 minutes

source("funct_plotperf.R")
plot.perf.sgccda.mthd(perf.diablo, sd = TRUE)

ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "max.dist"]


# Tuning keepX: tuning the number of variables to select in each data
set.seed(1)
test.keepX = list(rna = c(5:15),
                  mbiome = c(5:15))
tune.TCGA = tune.block.splsda(X = data, Y = group, ncomp = 2, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 3, nrepeat = 1000,
                              cpus = 4, dist = "max.dist")
list.keepX = tune.TCGA$choice.keepX


# Final model
sgccda.res = block.splsda(X = data, Y = group, ncomp = 2, 
                          keepX = list.keepX, design = design)
save(sgccda.res, file = "mixomics_finalmodel_full.Rda")

selectVar(sgccda.res, block = 'uArray', comp = 2)$uArray$name           # 15 genes
selectVar(sgccda.res, block = 'uBiome', comp = 3)$uBiome$name           # 5 genus
plotDiablo(sgccda.res, ncomp = 1)
plotIndiv(sgccda.res, ind.names = TRUE, legend = TRUE, 
          title = 'DIABLO')

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



## NULL DESIGN
sgccda.res <- block.splsda(X = data, Y = group, ncomp = 10, 
                           design = design_null)
set.seed(1)
t1 = proc.time()
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 3, nrepeat = 1000)
t2 = proc.time()
running_time = t2 - t1; running_time   # 36 minutes

source("funct_plotperf.R")
plot.perf.sgccda.mthd(perf.diablo, sd = TRUE)

ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "max.dist"]


# Tuning keepX: tuning the number of variables to select in each data
set.seed(1)
test.keepX = list(rna = c(5:15),
                  mbiome = c(5:15))
tune.TCGA = tune.block.splsda(X = data, Y = group, ncomp = 2, 
                              test.keepX = test.keepX, design = design_null,
                              validation = 'Mfold', folds = 3, nrepeat = 1000,
                              cpus = 4, dist = "max.dist")
list.keepX = tune.TCGA$choice.keepX


# Final model
sgccda.res = block.splsda(X = data, Y = group, ncomp = 2, 
                          keepX = list.keepX, design = design_null)
save(sgccda.res, file = "mixomics_finalmodel_null.Rda")

selectVar(sgccda.res, block = 'uArray', comp = 1)$uArray$name           # 15 genes
selectVar(sgccda.res, block = 'uBiome', comp = 2)$uBiome$name           # 5 genus
plotDiablo(sgccda.res, ncomp = 3)
png("mixomics_full_wo cma_3 comp.png")
plotIndiv(sgccda.res, ind.names = FALSE, legend = TRUE,
          style = 'ggplot2', title = 'DIABLO')
dev.off()

png("mixomics_null_arrow.png")
plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')
dev.off()

plotVar(sgccda.res, var.names = TRUE, style = 'graphics', legend = TRUE, 
        pch = c(16, 17), cex = c(1,1), col = c('brown1', 'lightgreen'))

png("mixomics_null_circos.png")
circosPlot(sgccda.res, cutoff = 0.6, line = TRUE, 
           color.blocks= c('brown1', 'lightgreen'),
           color.cor = c("chocolate3","grey20"), 
           size.variables = .6, size.labels = 1)
dev.off()

png("mixomics_null_network.png")
network(sgccda.res, blocks = c(1,2),
        color.node = c('brown1', 'lightgreen'), cutoff = 0.4)
dev.off()
plotLoadings(sgccda.res, comp = 2, contrib = 'max', method = 'median')







#----------------------------------------------------------------------------------------
# GET THE GENUS AND GENES NAME
#----------------------------------------------------------------------------------------
load("mixomics_finalmodel_full.Rda")


# Features
f.comp1 <- selectVar(sgccda.res, block = 'uArray', comp = 1)$uArray$name  
f.comp2 <- selectVar(sgccda.res, block = 'uArray', comp = 2)$uArray$name 

f.comp1 <- substring(f.comp1, 2)
f.comp2 <- substring(f.comp2, 2)

fdata <- as(featureData(esetRna), "data.frame")
fdata$SYMBOL[which(rownames(fdata) %in% f.comp1)]
fdata$SYMBOL[which(rownames(fdata) %in% f.comp2)]


# Genera
g.comp1 <- selectVar(sgccda.res, block = 'uBiome', comp = 1)$uBiome$name  
g.comp2 <- selectVar(sgccda.res, block = 'uBiome', comp = 2)$uBiome$name  

taxtable <- as.data.frame(tax_table(Obj_Genus))
taxtable$Genus[which(rownames(taxtable) %in% g.comp1)]
taxtable$Family[which(rownames(taxtable) == "SVs775")]
