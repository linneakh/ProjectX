# Linnea Honeker
# 2/15/22
# linneah@arizona.edu
# Objective: To perform WGCNA analysis on soil pyruvate metaT data paired with metabolomics/VOCs


library(WGCNA)
library(flashClust)
library(dplyr)
library(DESeq2)
library(pasilla)
library(qvalue)
library(dun.test)
library(DescTools)
library(ggplot2)
library(tidyverse)
library(nlme)
library(lmerTest)

options(stringsAsFactors = FALSE)

#load metaT gene copy data and column metadata
datExpr0 <-read.csv("./Data/Deseq2/Drought_vs_predrought_metaT_gene_copy.csv", header = TRUE)
rownames(datExpr0) <- datExpr0$Feature
datExpr0$Feature <- NULL

coldata <- read.csv("./Data/Deseq2/metaT-metadata.csv", header = TRUE)
rownames(coldata) <- coldata$SampleID
coldata$SampleID <- NULL
row.names(coldata) [7] <- "X3300045460_Site1_Drought_48hr"

#filter features that are not in at least 9samples
datExpr.f <- datExpr0 %>%
  mutate(count = rowSums(. > 0)) %>%
  filter(count > 9) %>%
  subset(select = -count )

#normalize filtered gene copy data using deseq2, variance stabilization transformation
dds <- DESeqDataSetFromMatrix(countData = datExpr.f,
                              colData = coldata,
                              design= ~ Condition_Time)

#vst transformation#
vstcounts <- vst(dds, blind=TRUE)
vst.export<- assay(vstcounts)
vst.export<- as.data.frame(vst.export)
write.csv(vst.export, "./Output/WGCNA/revision2/vst-metaT.csv")


datExpr <- t(vst.export)
names(datExpr) <- rownames(datExpr.f)
row.names(datExpr) [7] <- "X3300045460_Site1_Drought_48hr"

#not sure what this does
gsg = goodSamplesGenes(datExpr, verbose=3)
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# make a cluster tree of samples to detect outliers
sampleTree = hclust(dist(datExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 160, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 160, minSize = 5)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr[keepSamples, ]
names(datExpr) <- rownames(datExpr.f)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#load trait data - this is the output from merging VOC flux with expression data
traitData = read.csv("./Output/WGCNA/trait-data-VOC-flux-physicochemical.csv", header = TRUE)
traitData$X <- NULL
traitData$flux.acetate <- as.numeric(traitData$flux.acetate)
traitData$flux.acetone <- as.numeric(traitData$flux.acetone)
traitData$flux.diacetyl <- as.numeric(traitData$flux.diacetyl)
traitData$SM <- as.numeric(traitData$SM)
traitData$ST <- as.numeric(traitData$ST)
traitData$pH <- as.numeric(traitData$pH)

#merge trait and exp data in order to create a trait table with nas for missing values
datExpr <- as.data.frame(datExpr)
traitData <- as.data.frame(traitData)

merged_traitData <- datExpr %>%
  rownames_to_column(var = "SampleID") %>%
  merge(., traitData, by= "SampleID", all = TRUE) %>%
  dplyr::select(SampleID, flux.acetate, flux.acetone, flux.diacetyl, SM, ST, pH) %>%
  mutate(Condition = case_when(
    grepl("PreDrought", SampleID) ~ "0",
    grepl("Drought", SampleID) ~ "1",
    grepl("Pre_drought", SampleID) ~ "0"
  )) 

merged_traitData$Condition <- as.numeric(merged_traitData$Condition)

# Form a data frame analogous to expression data that will hold the clinical traits.
Samples = rownames(datExpr)
traitRows = match(Samples, merged_traitData$SampleID)
datTraits = merged_traitData[traitRows, -1]
rownames(datTraits) = merged_traitData[traitRows, 1]
collectGarbage()

# Re-cluster samples
sampleTree2 = flashClust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = TRUE, naColor = "grey");
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

# save input for wgcna analysis
save(datExpr, datTraits, file = "./Output/WGCNA/revision2/soil-pyruvate-input.RData")
load(file = "./Output/WGCNA/revision2/soil-pyruvate-input.RData")

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up # certain calculations
# and is optional.
# Any error here may be ignored but you may want to update WGCNA if you see one.
allowWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "./Output/WGCNA/revision2/soil-pyruvate-input.RData");
#The variable lnames contains the names of loaded variables.
lnames
## Automatic network construction and module detection

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#create WGCNA TOM network, power = 6 chosen becuase this is where R2 passed 0.8
net = blockwiseModules(datExpr, checkMissingData = TRUE, power = 6, minModuleSize = 80 ,
                       numericLabels = TRUE,
                       pamRespectsDendro = FALSE,
                       corTypr = "bicor",
                       maxPOutliers = 0.1,
                       networkType = "signed",
                       #loadTOM = TRUE,
                       TOMType = "signed",
                       saveTOMs = TRUE,
                       saveTOMFileBase = "./Output/WGCNA/revision2/SoilPyruvateTOM", 
                       minKMEtoStay = 0.35,
                       verbose = 3)
table(net$colors)

#load TOM created previously
#create WGCNA TOM network, power = 6 chosen becuase this is where R2 passed 0.8
net = blockwiseModules(datExpr, checkMissingData = TRUE, power = 6, minModuleSize = 80 ,
                       numericLabels = TRUE,
                       pamRespectsDendro = FALSE,
                       corTypr = "bicor",
                       maxPOutliers = 0.1,
                       networkType = "signed",
                       loadTOM = TRUE,
                       TOMType = "signed",
                       #saveTOMs = TRUE,
                       saveTOMFileBase = "./Output/WGCNA/revision2/SoilPyruvateTOM", 
                       minKMEtoStay = 0.35,
                       verbose = 3)


# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath, save image (Fig S7a)
filename=paste("./Figures/Fig4-S8-WGCNA-ClusterProfiler/WGCNA/revision2/dendogram.pdf",sep="")
pdf(filename ,width=7, height=6)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#save network construction data
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "./Output/WGCNA/revision2/networkConstruction-auto.RData")
load(file="./Output/WGCNA/revision2/networkConstruction-auto.RData")

##Relating modules to external information and identifying important genes
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

#Erase grey column from MEs 
MEs$MEgrey <- NULL

######All modules#####
moduleTraitCor = cor(MEs, datTraits, use = "p");

#calculate correlations with subset of traitdata
datTraits.subset <- datTraits[,c(7, 1, 2)]
moduleTraitCor.subset <- cor(MEs, datTraits.subset, use = "p")

#calculate correlations and pvalues, create each set of p values as a vector to caclulate fdr-corrected q values
moduleTraitPvalue = corPvalueStudent(moduleTraitCor.subset, nSamples);
moduleTraitPvector= c(moduleTraitPvalue)
moduleTraitQobj = qvalue(moduleTraitPvector, lambda = 0,  fdr.level = 0.05)

moduleTraitQvector= moduleTraitQobj$qvalues

moduleTraitQvalues = matrix(moduleTraitQvector,
                            nrow=nrow(moduleTraitCor.subset), ncol=ncol(moduleTraitCor.subset), 
                            byrow=FALSE, 
                            dimnames=list(names(MEs),
                                          names(datTraits.subset)))

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor.subset, 2), "\n(",
                   signif(moduleTraitQvalues, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor.subset)
par(mar = c(6, 8.5, 3, 3));

#Display the correlation values within a heatmap plot
filename=paste("./Figures/Fig4-S8-WGCNA-ClusterProfiler/WGCNA/revision2/modules-trait-relationships_all.pdf",sep="")
pdf(filename ,width=6, height=6)
labeledHeatmap(Matrix = moduleTraitCor.subset,
               #xLabels = c("acetate-C2", "acetone-C2", "diacetyl-C2", "SM", "ST", "pH"),
               xLabels = colnames(moduleTraitCor.subset),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               cex.text = 1.0,
               cex.lab.y = 1.2,
               cex.lab.x = 1.2,
               xLabelsAngle = 0,
               xLabelsAdj = 0.5,
               xLabelsPosition = "bottom",
               zlim = c(-1,1))
dev.off()

######Subset modules to those discussed in paper (Fig S7)#####
MEs_order = MEs[,c(9,8,7,1)]

#correlate traits and MEs
moduleTraitCor.subset.fig = cor(MEs_order, datTraits.subset, use = "p");

#calculate correlations and pvalues, create each set of p values as a vector to caclulate fdr-corrected q values
moduleTraitPvalue = corPvalueStudent(moduleTraitCor.subset.fig, nSamples);
moduleTraitPvector= c(moduleTraitPvalue)
moduleTraitQobj = qvalue(moduleTraitPvector, lambda = 0,  fdr.level = 0.05)

moduleTraitQvector= moduleTraitQobj$qvalues

moduleTraitQvalues = matrix(moduleTraitQvector,
                            nrow=nrow(moduleTraitCor.subset.fig), ncol=ncol(moduleTraitCor.subset.fig), 
                            dimnames=list(names(MEs_order),
                                          names(datTraits.subset)))


# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor.subset.fig, 2), "\n(",
                   signif(moduleTraitQvalues, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor.subset.fig)
par(mar = c(6, 8.5, 3, 3));

#Display the correlation values within a heatmap plot
filename=paste("./Figures/Fig4-S8-WGCNA-ClusterProfiler/WGCNA/revision2/modules-trait-relationships_subset_FigS8b.pdf",sep="")
pdf(filename ,width=6, height=6)
labeledHeatmap(Matrix = moduleTraitCor.subset.fig,
               xLabels = c("Condition",  
                           "Acetate-C2", "Acetone-C2"),
               yLabels = c("Pink", "Green", "Magenta", "Brown"),
               #yLabels = names(MEs_order),
               ySymbols = names(MEs_order),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               #setStdMargins = FALSE,
               cex.text = 1.0,
               cex.lab.y = 1.0,
               cex.lab.x = 1.0,
               xLabelsAngle = 0,
               xLabelsAdj = 0.5,
               xLabelsPosition = "bottom",
               zlim = c(-1,1))
dev.off()




# Plot the relationships among the eigengenes as png
filename=paste("./Figures/Fig5-S7-WGCNA-ClusterProfiler/WGCNA/revision2/modules_eigengene3-subset.png",sep="")
png(filename, width=12, height=9, unit = "in", res = 1000)
plotEigengeneNetworks(MEs_order,"",marDendro=c(0,4,1,6), 
                      marHeatmap=c(3,4,1,2), cex.lab=0.8,xLabelsAngle=90)
dev.off()


#####Condition#####
#Define variable Condition
Condition = as.data.frame(datTraits$Condition);
names(Condition) = ("Condition")
# names (colors) of the modules
modNames = substring(names(MEs0), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs0, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, Condition, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
#GSPvalue = corPvalueStudent(as.matrix(geneTraitSignificance), nSamples)
GSPvalue.0 = as.vector(GSPvalue[,'Condition'])
GSQvalue.0 = qvalue(p=GSPvalue.0, lambda = 0, fdr.level = 0.05)
GSQvalue <- as.data.frame(GSQvalue.0$qvalues, row.names=names(datExpr))
summary(GSQvalue.0)

pi0 <- GSQvalue.0$pi0
lfdr <- GSQvalue.0$lfdr 

moduleColor=names(MEs0)
probes=names(datExpr)

# Create the starting data frame
geneInfo0 = data.frame(datExpr = probes, ModuleColor = moduleColors,
                       geneTraitSignificance,
                       GSQvalue)

# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs0, Condition, use = "p")));

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
write.csv(geneInfo0, file = "./Output/WGCNA/geneInfo-Condition-not-ordered.csv")

#gene screening method bsed on detailied definition of module membership
NS1=networkScreening(y=datTraits$Condition, datME=MEs0, datExpr=datExpr,
                     oddPower=3, blockSize=1000, minimumSampleSize=4,
                     addMEy=TRUE, removeDiag=FALSE, weightESy=0.5, corOptions = "use = 'p', method = 'spearman'")
write.csv(NS1, "./Output/WGCNA/gene-screening-condition-spearman.csv")

#####Acetate-c2#####
#Define variable AcetateC2 containing the Acetate.C2 column of datTrait
AcetateC2 = as.data.frame(datTraits$Acetate.C2);
names(AcetateC2) = ("AcetateC2")
# names (colors) of the modules
modNames = substring(names(MEs0), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs0, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, AcetateC2, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
#GSPvalue = corPvalueStudent(as.matrix(geneTraitSignificance), nSamples)
GSPvalue.0 = as.vector(GSPvalue[,'AcetateC2'])
GSQvalue.0 = qvalue(p=GSPvalue.0, lambda = 0, fdr.level = 0.05)
GSQvalue <- as.data.frame(GSQvalue.0$qvalues, row.names=names(datExpr))
summary(GSQvalue.0)

pi0 <- GSQvalue.0$pi0
lfdr <- GSQvalue.0$lfdr 

# Create the starting data frame
geneInfo0 = data.frame(datExpr = probes, ModuleColor = moduleColors,
                       geneTraitSignificance,
                       GSQvalue)

# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, AcetateC2, use = "p")));

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
write.csv(geneInfo0, file = "geneInfo-not-Acetatec2-ordered-q.csv")

#gene screening method bsed on detailied definition of module membership
NS1=networkScreening(y=datTraits$Acetate.C2, datME=MEs0, datExpr=datExpr,
                     oddPower=3, blockSize=1000, minimumSampleSize=4,
                     addMEy=TRUE, removeDiag=FALSE, weightESy=0.5, corOptions = "use = 'p', method = 'spearman'")
write.csv(NS1, "./Output/WGCNA/gene-screening-AcetateC2-spearman.csv")


#####acetone-c2#####
#Define variable AcetateC2 containing the Acetate.C2 column of datTrait
AcetoneC2 = as.data.frame(datTraits$Acetone.C2);
names(AcetoneC2) = ("AcetoneC2")
# names (colors) of the modules
modNames = substring(names(MEs0), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs0, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, AcetoneC2, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
#GSPvalue = corPvalueStudent(as.matrix(geneTraitSignificance), nSamples)
GSPvalue.0 = as.vector(GSPvalue[,'AcetoneC2'])
GSQvalue.0 = qvalue(p=GSPvalue.0, lambda = 0, fdr.level = 0.05)
GSQvalue <- as.data.frame(GSQvalue.0$qvalues, row.names=names(datExpr))
summary(GSQvalue.0)


pi0 <- GSQvalue.0$pi0
lfdr <- GSQvalue.0$lfdr 

# Create the starting data frame
geneInfo0 = data.frame(datExpr = probes, ModuleColor = moduleColors,
                       geneTraitSignificance,
                       GSQvalue)

# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs0, AcetoneC2, use = "p")));

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod
                                                                    ]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance

write.csv(geneInfo0, file = "./Output/WGCNA/geneInfo-not-AcetoneC2-ordered-q.csv")

#gene screening method bsed on detailied definition of module membership
NS1=networkScreening(y=datTraits$Acetone.C2, datME=MEs0, datExpr=datExpr,
                     oddPower=3, blockSize=1000, minimumSampleSize=4,
                     addMEy=TRUE, removeDiag=FALSE, weightESy=0.5, corOptions = "use = 'p', method = 'spearman'")
write.csv(NS1, "gene-screening-AcetoneC2-spearman.csv")


##plot eigengene expression per condition_time (Fig 5)
coldata$Condition_Time <- factor(coldata$Condition_Time, c("PreDrought_0hr", "PreDrought_6hr", "PreDrought_48hr",
                                                           "Drought_0hr", "Drought_6hr", "Drought_48hr"))
coldata$Condition <- factor(coldata$Condition, c("PreDrought", "Drought"))
coldata$Time <- factor(coldata$Time, c("0hr", "6hr", "48hr"))

coldata$SampleID <- rownames(coldata)

###plot eigengene expression per sample,
filename=paste("./Figures/Fig5-S7-WGCNA-ClusterProfiler/WGCNA/black-EGexp.png", sep = "")
png(filename ,width=3, height=2.5, unit='in', res = 300)
which.module="black"
signif(cor(MEs, use="p"), 2)
ME.black=MEs[, paste("ME",which.module, sep="")]
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME.black, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
dev.off()


filename=paste("./Figures/Fig5-S7-WGCNA-ClusterProfiler/WGCNA/magenta-EGexp.png", sep = "")
png(filename ,width=3, height=2.5, unit='in', res = 300)
which.module="magenta"
signif(cor(MEs, use="p"), 2)
ME.magenta=MEs[, paste("ME",which.module, sep="")]
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME.magenta, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
dev.off()

filename=paste("./Figures/Fig5-S7-WGCNA-ClusterProfiler/WGCNA/blue-EGexp.png", sep = "")
png(filename ,width=3, height=2.5, unit='in', res = 300)
which.module="blue"
signif(cor(MEs, use="p"), 2)
ME.blue=MEs[, paste("ME",which.module, sep="")]
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME.blue, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
dev.off()

filename=paste("./Figures/Fig5-S7-WGCNA-ClusterProfiler/WGCNA/yellow-EGexp.png", sep = "")
png(filename ,width=3, height=2.5, unit='in', res = 300)
which.module="yellow"
signif(cor(MEs, use="p"), 2)
ME.yellow=MEs[, paste("ME",which.module, sep="")]
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME.yellow, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
dev.off()

filename=paste("./Figures/Fig5-S7-WGCNA-ClusterProfiler/WGCNA/brown-EGexp.png", sep = "")
png(filename ,width=3, height=2.5, unit='in', res = 300)
which.module="brown"
signif(cor(MEs, use="p"), 2)
ME.brown=MEs[, paste("ME",which.module, sep="")]
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME.brown, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
dev.off()

filename=paste("./Figures/Fig5-S7-WGCNA-ClusterProfiler/WGCNA/turquoise-EGexp.png", sep = "")
png(filename ,width=3, height=2.5, unit='in', res = 300)
which.module="turquoise"
signif(cor(MEs, use="p"), 2)
ME.turquoise=MEs[, paste("ME",which.module, sep="")]
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME.turquoise, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
dev.off()

filename=paste("./Figures/Fig5-S7-WGCNA-ClusterProfiler/WGCNA/pink-EGexp.png", sep = "")
png(filename ,width=3, height=2.5, unit='in', res = 300)
which.module="pink"
signif(cor(MEs, use="p"), 2)
ME.pink=MEs[, paste("ME",which.module, sep="")]
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME.pink, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
dev.off()

filename=paste("./Figures/Fig5-S7-WGCNA-ClusterProfiler/WGCNA/green-EGexp.png", sep = "")
png(filename ,width=3, height=2.5, unit='in', res = 300)
which.module="green"
signif(cor(MEs, use="p"), 2)
ME.green=MEs[, paste("ME",which.module, sep="")]
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME.green, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
dev.off()

#black
ME.black <- as.data.frame(ME.black)
rownames(ME.black) <- rownames(datExpr)
ME.black$SampleID <- rownames(ME.black)
ME.black.m <- ME.black %>%
  merge(.,coldata, by = "SampleID")


ggplot(ME.black.m, aes(x=Time, y= ME.black)) +
  geom_point() + geom_boxplot(fill = "darkgrey") +
  facet_grid(~Condition) +
  #ylim(-0.35, 0.4) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  labs(title = "black module", y = "Eigengene Expression", x= "Time post pyruvate addition (h)") +
  theme_bw() +
  theme(
        axis.title.x = element_text(size = 10, color = "black", face = "bold"),
        axis.title.y = element_text(size = 10, color = "black", face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
filename=paste("./Figures/Fig5-S7-WGCNA-ClusterProfiler/WGCNA/black-EGexp-facet.pdf", sep = "")
ggsave(filename,width=4,height=2,units="in",dpi=300)



#magenta
ME.magenta <- as.data.frame(ME.magenta)
rownames(ME.magenta) <- rownames(datExpr)
ME.magenta$SampleID <- rownames(ME.magenta)
ME.magenta.m <- ME.magenta %>%
  merge(.,coldata, by = "SampleID")

write.csv(ME.magenta.m, "./Output/WGCNA/ME-magenta.csv")

ggplot(ME.magenta.m, aes(x=Time, y= ME.magenta)) +
  geom_boxplot(fill = "magenta") + geom_point() +
  facet_grid(~Condition) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  ylim(-0.3, 0.7) +
  labs(title = "magenta module", y = "Eigengene Expression", x= "Time post pyruvate addition (h)") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 10, color = "black", face = "bold"),
    axis.title.y = element_text(size = 10, color = "black", face = "bold"),
    axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10))
filename=paste("./Figures/Fig5-S7-WGCNA-ClusterProfiler/WGCNA/Fig5-magenta-EGexp-facet.pdf", sep = "")
ggsave(filename,width=4,height=2,units="in",dpi=300)

#blue
ME.blue <- as.data.frame(ME.blue)
rownames(ME.blue) <- rownames(datExpr)
ME.blue$SampleID <- rownames(ME.blue)
ME.blue.m <- ME.blue %>%
  merge(.,coldata, by = "SampleID")

ggplot(ME.blue.m, aes(x=Time, y= ME.blue)) +
  geom_point() + geom_boxplot(fill = "blue") +
  facet_grid(~Condition) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  labs(title = "blue module", y = "Eigengene Expression", x= "Time post pyruvate addition (h)") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 10, color = "black", face = "bold"),
    axis.title.y = element_text(size = 10, color = "black", face = "bold"),
    axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10))
filename=paste("./Figures/Fig5-S7-WGCNA-ClusterProfiler/WGCNA/blue-EGexp-facet.pdf", sep = "")
ggsave(filename,width=4,height=2,units="in",dpi=300)

#yellow
ME.yellow <- as.data.frame(ME.yellow)
rownames(ME.yellow) <- rownames(datExpr)
ME.yellow$SampleID <- rownames(ME.yellow)
ME.yellow.m <- ME.yellow %>%
  merge(.,coldata, by = "SampleID")

ggplot(ME.yellow.m, aes(x=Time, y= ME.yellow)) +
  geom_point() + geom_boxplot(fill = "yellow") +
  facet_grid(~Condition) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  labs(title = "yellow module", y = "Eigengene Expression", x= "Time post pyruvate addition (h)") +
  ylim(-0.3, 0.4) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 10, color = "black", face = "bold"),
    axis.title.y = element_text(size = 10, color = "black", face = "bold"),
    axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10))
filename=paste("./Figures/Fig5-S7-WGCNA-ClusterProfiler/WGCNA/yellow-EGexp-facet.pdf", sep = "")
ggsave(filename,width=4,height=2,units="in",dpi=300)

#brown
ME.brown <- as.data.frame(ME.brown)
rownames(ME.brown) <- rownames(datExpr)
ME.brown$SampleID <- rownames(ME.brown)
ME.brown.m <- ME.brown %>%
  merge(.,coldata, by = "SampleID")

write.csv(ME.brown.m, "./Output/WGCNA/ME-brown.csv")


ggplot(ME.brown.m, aes(x=Time, y= ME.brown)) +
   geom_boxplot(fill = "brown") +geom_point() +
  facet_grid(~Condition) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  ylim(-.5, .4) +
  labs(title = "brown module", y = "Eigengene Expression", x= "Time post pyruvate addition (h)") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 10, color = "black", face = "bold"),
    axis.title.y = element_text(size = 10, color = "black", face = "bold"),
    axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10))
filename=paste("./Figures/Fig5-S7-WGCNA-ClusterProfiler/WGCNA/Fig5-brown-EGexp-facet.pdf", sep = "")
ggsave(filename,width=4,height=2,units="in",dpi=300)

#turquoise
ME.turquoise <- as.data.frame(ME.turquoise)
rownames(ME.turquoise) <- rownames(datExpr)
ME.turquoise$SampleID <- rownames(ME.turquoise)
ME.turquoise.m <- ME.turquoise %>%
  merge(.,coldata, by = "SampleID")

ggplot(ME.turquoise.m, aes(x=Time, y= ME.turquoise)) +
  geom_point() + geom_boxplot(fill = "turquoise") +
  facet_grid(~Condition) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  labs(title = "turquoise module", y = "Eigengene Expression", x= "Time post pyruvate addition (h)") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 10, color = "black", face = "bold"),
    axis.title.y = element_text(size = 10, color = "black", face = "bold"),
    axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10))
filename=paste("./Figures/Fig5-S7-WGCNA-ClusterProfiler/WGCNA/turquoise-EGexp-facet.pdf", sep = "")
ggsave(filename,width=4,height=2,units="in",dpi=300)

#pink
ME.pink <- as.data.frame(ME.pink)
rownames(ME.pink) <- rownames(datExpr)
ME.pink$SampleID <- rownames(ME.pink)
ME.pink.m <- ME.pink %>%
  merge(.,coldata, by = "SampleID")

write.csv(ME.pink.m, "./Output/WGCNA/ME-pink.csv")


ggplot(ME.pink.m, aes(x=Time, y= ME.pink)) +
  geom_boxplot(fill = "pink") + geom_point() +
  facet_grid(~Condition) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  ylim(-0.4, 0.4) +
  labs(title = "pink module", y = "Eigengene Expression", x= "Time post pyruvate addition (h)") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 10, color = "black", face = "bold"),
    axis.title.y = element_text(size = 10, color = "black", face = "bold"),
    axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10))
filename=paste("./Figures/Fig5-S7-WGCNA-ClusterProfiler/WGCNA/Fig5-pink-EGexp-facet.pdf", sep = "")
ggsave(filename,width=4,height=2,units="in",dpi=300)

#green
ME.green <- as.data.frame(ME.green)
rownames(ME.green) <- rownames(datExpr)
ME.green$SampleID <- rownames(ME.green)
ME.green.m <- ME.green %>%
  merge(.,coldata, by = "SampleID")

write.csv(ME.green.m, "./Output/WGCNA/ME-green.csv")


ggplot(ME.green.m, aes(x=Time, y= ME.green)) +
  geom_boxplot(fill = "green") + geom_point() + 
  facet_grid(~Condition) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  ylim(-0.4, 0.5) +
  labs(title = "green module", y = "Eigengene Expression", x= "Time post pyruvate addition (h)") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 10, color = "black", face = "bold"),
    axis.title.y = element_text(size = 10, color = "black", face = "bold"),
    axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10))
filename=paste("./Figures/Fig5-S7-WGCNA-ClusterProfiler/WGCNA/Fig5-green-EGexp-facet.pdf", sep = "")
ggsave(filename,width=4,height=2,units="in",dpi=300)

##statistics on eigenggene expression
##### linear mixed effect models############
###Pink#####
# Condition
lme.pink <- lme(ME.pink ~ Condition,
                random = list(Site = ~1),
                data = ME.pink.m,
                weights =  varIdent(form = ~1|Condition)
)
sum <- summary(lme.pink)
tabl = sum$tTable
tabl
#Value  Std.Error DF   t-value      p-value
#(Intercept)      -0.137531 0.02126949 30 -6.466116 3.818768e-07
#ConditionDrought  0.275062 0.03622325 30  7.593521 1.813810e-08

# Time within pre-drought,
ME.pink.m.p.6 <- ME.pink.m %>%
  filter(Condition == "PreDrought" )

lme.pink.p.6 <- lme(ME.pink ~ Time,
                  random = list(Site = ~1),
                  data = ME.pink.m.p.6,
                  weights =  varIdent(form = ~1|Time)
)
sum <- summary(lme.pink.p.6)
tabl = sum$tTable
tabl



# Time within drought,
ME.pink.m.d.6 <- ME.pink.m %>%
  filter(Condition == "Drought" )

lme.pink.d.6 <- lme(ME.pink ~ Time,
                    random = list(Site = ~1),
                    data = ME.pink.m.d.6,
                    weights =  varIdent(form = ~1|Time)
)
sum <- summary(lme.pink.d.6)
tabl = sum$tTable
tabl


###magenta#####
# Condition
lme.magenta <- lme(ME.magenta ~ Condition,
                   random = list(Site = ~1),
                   data = ME.magenta.m,
                   weights =  varIdent(form = ~1|Condition)
)
sum <- summary(lme.magenta)
tabl = sum$tTable
tabl
#                     Value  Std.Error DF   t-value    p-value
#(Intercept)      -0.05963747 0.03799614 29 -1.569567 0.12736418
#ConditionDrought  0.13369675 0.05304906 29  2.520247 0.01748625

# Time within pre-drought,
ME.magenta.m.p.6 <- ME.magenta.m %>%
  filter(Condition == "PreDrought")

lme.magenta.p.6 <- lme(ME.magenta ~ Time,
                    random = list(Site = ~1),
                    data = ME.magenta.m.p.6,
                    weights =  varIdent(form = ~1|Time)
)
sum <- summary(lme.magenta.p.6)
tabl = sum$tTable
tabl
#                Value  Std.Error DF    t-value    p-value
#(Intercept) -0.00220941 0.05278143 12 -0.0418596 0.96729905
#Time6hr     -0.03285108 0.02266492 12 -1.4494238 0.17284584
#Time48hr    -0.12927503 0.04994979 12 -2.5880995 0.02374089



# Time within drought,
ME.magenta.m.d.48 <- ME.magenta.m %>%
  filter(Condition == "Drought")

lme.magenta.d.48 <- lme(ME.magenta ~ Time,
                     random = list(Site = ~1),
                     data = ME.magenta.m.d.48,
                     weights =  varIdent(form = ~1|Time)
)
sum <- summary(lme.magenta.d.48)
tabl = sum$tTable
tabl

#                Value Std.Error DF   t-value   p-value
#(Intercept)  0.2185654 0.1255910 11  1.740295 0.1096692
#Time6hr     -0.2075054 0.1319343 11 -1.572793 0.1440684
#Time48hr    -0.2114322 0.1266594 11 -1.669298 0.1232385

###green###
# Condition
lme.green <- lme(ME.green ~ Condition,
                 random = list(Site = ~1),
                 data = ME.green.m,
                 weights =  varIdent(form = ~1|Condition)
)
sum <- summary(lme.green)
tabl = sum$tTable
tabl
#(Intercept)      -0.09530071 0.03027621 29 -3.147709 0.003791218
#ConditionDrought  0.17796105 0.04963319 29  3.585525 0.001216596

# Time within pre-drought, 
ME.green.m.p.6 <- ME.green.m %>%
  filter(Condition == "PreDrought")

lme.green.p.6 <- lme(ME.green ~ Time,
                       random = list(Site = ~1),
                       data = ME.green.m.p.6,
                       weights =  varIdent(form = ~1|Time)
)
sum <- summary(lme.green.p.6)
tabl = sum$tTable
tabl

#                 Value  Std.Error DF    t-value     p-value
#(Intercept) -0.12780851 0.03449529 12 -3.7051001 0.003007612
#Time6hr      0.07526111 0.03155674 12  2.3849453 0.034451097
#Time48hr     0.03105081 0.03196725 12  0.9713317 0.350555474

# Time within drought,
ME.green.m.d.6 <- ME.green.m %>%
  filter(Condition == "Drought")

lme.green.d.6 <- lme(ME.green ~ Time,
                       random = list(Site = ~1),
                       data = ME.green.m.d.6,
                       weights =  varIdent(form = ~1|Time)
)
sum <- summary(lme.green.d.6)
tabl = sum$tTable
tabl

#               Value  Std.Error DF    t-value    p-value
#(Intercept)  0.01616096 0.07690994 11  0.2101283 0.83740895
#Time6hr     -0.02253652 0.07990887 11 -0.2820278 0.78315768
#Time48hr     0.23177612 0.12585236 11  1.8416510 0.09262781

###brown###
# Condition
lme.brown <- lme(ME.brown ~ Condition,
                 random = list(Site = ~1),
                 data = ME.brown.m,
                 weights =  varIdent(form = ~1|Condition)
)
sum <- summary(lme.brown)
tabl = sum$tTable
tabl
#Value  Std.Error DF   t-value      p-value
#(Intercept)       0.1028595 0.03595102 29  2.861100 0.0077527328
#ConditionDrought -0.1960638 0.04888601 29 -4.010632 0.0003886758

# Time within pre-drought
ME.brown.m.p.6 <- ME.brown.m %>%
  filter(Condition == "PreDrought")

lme.brown.p.6 <- lme(ME.brown ~ Time,
                     random = list(Site = ~1),
                     data = ME.brown.m.p.6,
                     weights =  varIdent(form = ~1|Time)
)
sum <- summary(lme.brown.p.6)
tabl = sum$tTable
tabl


# Time within drought,
ME.brown.m.d.6 <- ME.brown.m %>%
  filter(Condition == "Drought")

lme.brown.d.6 <- lme(ME.brown ~ Time,
                     random = list(Site = ~1),
                     data = ME.brown.m.d.6,
                     weights =  varIdent(form = ~1|Time)
)
sum <- summary(lme.brown.d.6)
tabl = sum$tTable
tabl




### red ###
pairwise.wilcox.test(ME.red.m$ME.red, ME.red.m$Condition_Time, p.adjust.method = "BH")
#data:  ME.red.m$ME.red and ME.red.m$Condition_Time 

#PreDrought_0hr PreDrought_6hr PreDrought_48hr Drought_0hr Drought_6hr
#PreDrought_6hr  0.0519         -              -               -           -          
#PreDrought_48hr 0.1320         0.7922         -               -           -          
#  Drought_0hr     0.9307         0.1508         0.4286          -           -          
#  Drought_6hr     0.1797         0.4286         0.6991          0.3290      -          
#  Drought_48hr    0.0087         0.1255         0.2403          0.1255      0.1797     

#P value adjustment method: none 
#data:  ME.red.m$ME.red and ME.red.m$Condition_Time 

#             PreDrought_0hr PreDrought_6hr PreDrought_48hr Drought_0hr Drought_6hr
#PreDrought_6hr  0.34           -              -               -           -          
#  PreDrought_48hr 0.34           0.85           -               -           -          
#  Drought_0hr     0.93           0.34           0.54            -           -          
#  Drought_6hr     0.34           0.54           0.81            0.49        -          
#  Drought_48hr    0.13           0.34           0.40            0.34        0.34       

#P value adjustment method: BH 

pairwise.wilcox.test(ME.red.m$ME.red, ME.red.m$Condition, p.adjust.method = "none")
#data:  ME.red.m$ME.red and ME.red.m$Condition 

#        PreDrought
#Drought 0.45     

#P value adjustment method: none 

## brown ##
pairwise.wilcox.test(ME.brown.m$ME.brown, ME.brown.m$Condition_Time, p.adjust.method = "BH")
#data:  ME.brown.m$ME.brown and ME.brown.m$Condition_Time 

#               PreDrought_0hr PreDrought_6hr PreDrought_48hr Drought_0hr Drought_6hr
#PreDrought_6hr  0.0043         -              -               -           -          
#  PreDrought_48hr 0.0087         0.7922         -               -           -          
#  Drought_0hr     0.0043         0.0079         0.1775          -           -          
#  Drought_6hr     0.0022         0.0173         0.3095          0.5368      -          
#  Drought_48hr    0.0022         0.0173         0.1797          0.6623      0.5887     

#P value adjustment method: none

#data:  ME.brown.m$ME.brown and ME.brown.m$Condition_Time 

#             PreDrought_0hr PreDrought_6hr PreDrought_48hr Drought_0hr Drought_6hr
#PreDrought_6hr  0.016          -              -               -           -          
#  PreDrought_48hr 0.022          0.792          -               -           -          
#  Drought_0hr     0.016          0.022          0.269           -           -          
#  Drought_6hr     0.016          0.032          0.422           0.671       -          
#  Drought_48hr    0.016          0.032          0.269           0.710       0.679      

#P value adjustment method: BH 

pairwise.wilcox.test(ME.brown.m$ME.brown, ME.brown.m$Condition, p.adjust.method = "none")
#data:  ME.brown.m$ME.brown and ME.brown.m$Condition 

#      PreDrought
#Drought 2.8e-05   

#P value adjustment method: none 



## green ##
#pairwise.wilcox.test(ME.green.m$ME.green, ME.green.m$Condition_Time, p.adjust.method = "BH")

#data:  ME.green.m$ME.green and ME.green.m$Condition_Time 

#               PreDrought_0hr PreDrought_6hr PreDrought_48hr Drought_0hr Drought_6hr
#PreDrought_6hr  0.052          -              -               -           -          
#  PreDrought_48hr 0.310          0.177          -               -           -          
##  Drought_0hr     0.126          0.548          0.177           -           -          
#  Drought_6hr     0.026          0.429          0.065           0.792       -          
#  Drought_48hr    0.015          0.082          0.041           0.082       0.093      

#P value adjustment method: none 

#data:  ME.green.m$ME.green and ME.green.m$Condition_Time 

#            PreDrought_0hr PreDrought_6hr PreDrought_48hr Drought_0hr Drought_6hr
#PreDrought_6hr  0.17           -              -               -           -          
#  PreDrought_48hr 0.39           0.24           -               -           -          
#  Drought_0hr     0.21           0.59           0.24            -           -          
#  Drought_6hr     0.17           0.49           0.17            0.79        -          
#  Drought_48hr    0.17           0.17           0.17            0.17        0.17       

#P value adjustment method: BH

pairwise.wilcox.test(ME.green.m$ME.green, ME.green.m$Condition, p.adjust.method = "none")
#data: ME.green.m$ME.green and ME.green.m$Condition

#      PreDrought
#Drought 0.0015   

#P value adjustment method: none 



### blue ###
pairwise.wilcox.test(ME.blue.m$ME.blue, ME.blue.m$Condition_Time, p.adjust.method = "BH")
#data:  ME.blue.m$ME.blue and ME.blue.m$Condition_Time 

#                PreDrought_0hr PreDrought_6hr PreDrought_48hr Drought_0hr Drought_6hr
#PreDrought_6hr    0.0519         -              -               -           -          
# PreDrought_48hr 0.5887         0.3290         -               -           -          
#  Drought_0hr     0.5368         0.0079         0.3290          -           -          
#  Drought_6hr     0.6991         0.0303         0.6991          0.3290      -          
#  Drought_48hr    0.9372         0.0087         0.5887          0.5368      0.5887     

#P value adjustment method: none 

#data:  ME.blue.m$ME.blue and ME.blue.m$Condition_Time 

#          PreDrought_0hr PreDrought_6hr PreDrought_48hr Drought_0hr Drought_6hr
#PreDrought_6hr  0.195          -              -               -           -          
#  PreDrought_48hr 0.736          0.705          -               -           -          
#  Drought_0hr     0.736          0.065          0.705           -           -          
#  Drought_6hr     0.749          0.152          0.749           0.705       -          
#  Drought_48hr    0.937          0.065          0.736           0.736       0.736      

#P value adjustment method: BH 

pairwise.wilcox.test(ME.blue.m$ME.blue, ME.blue.m$Condition, p.adjust.method = "none")
#data:  ME.blue.m$ME.blue and ME.blue.m$Condition 

#        PreDrought
#Drought 0.079   

#P value adjustment method: none 

### magenta ####
pairwise.wilcox.test(ME.magenta.m$ME.magenta, ME.magenta.m$Condition_Time, p.adjust.method = "BH")
#data:  ME.magenta.m$ME.magenta and ME.magenta.m$Condition_Time 

#               PreDrought_0hr PreDrought_6hr PreDrought_48hr Drought_0hr Drought_6hr
#PreDrought_6hr  0.537          -              -               -           -          
#  PreDrought_48hr 0.041          0.177          -               -           -          
#  Drought_0hr     0.177          0.095          0.052           -           -          
#  Drought_6hr     0.699          0.329          0.065           0.247       -          
#  Drought_48hr    0.699          0.931          0.093           0.126       0.589      

#P value adjustment method: none 

#data:  ME.magenta.m$ME.magenta and ME.magenta.m$Condition_Time 

#PreDrought_0hr PreDrought_6hr PreDrought_48hr Drought_0hr Drought_6hr
#PreDrought_6hr  0.73           -              -               -           -          
#  PreDrought_48hr 0.29           0.33           -               -           -          
#  Drought_0hr     0.33           0.29           0.29            -           -          
#  Drought_6hr     0.75           0.49           0.29            0.41        -          
#  Drought_48hr    0.75           0.93           0.29            0.31        0.74       

#P value adjustment method: BH 

pairwise.wilcox.test(ME.magenta.m$ME.magenta, ME.magenta.m$Condition, p.adjust.method = "none")
#data:  ME.magenta.m$ME.magenta and ME.magenta.m$Condition 

#PreDrought
#Drought 0.049     

#P value adjustment method: none 

### pink ###
#pairwise wilcox text
pairwise.wilcox.test(ME.pink.m$ME.pink, ME.pink.m$Condition_Time, p.adjust.method = "BH")
#data:  ME.pink.m$ME.pink and ME.pink.m$Condition_Time 

#              PreDrought_0hr PreDrought_6hr PreDrought_48hr Drought_0hr Drought_6hr
#PreDrought_6hr  0.0303         -              -               -           -          
#  PreDrought_48hr 0.0931         0.9307         -               -           -          
#  Drought_0hr     0.0043         0.0079         0.0087          -           -          
#  Drought_6hr     0.0043         0.0303         0.1320          0.0173      -          
#  Drought_48hr    0.0022         0.0043         0.0043          0.1255      0.1797     

#P value adjustment method: none 

#data:  ME.pink.m$ME.pink and ME.pink.m$Condition_Time 

#               PreDrought_0hr PreDrought_6hr PreDrought_48hr Drought_0hr Drought_6hr
#PreDrought_6hr  0.045          -              -               -           -          
#  PreDrought_48hr 0.127          0.931          -               -           -          
#  Drought_0hr     0.013          0.019          0.019           -           -          
#  Drought_6hr     0.013          0.045          0.152           0.032       -          
#  Drought_48hr    0.013          0.013          0.013           0.152       0.192      

# P value adjustment method: BH 

pairwise.wilcox.test(ME.pink.m$ME.pink, ME.pink.m$Condition, p.adjust.method = "none")
#data:  ME.pink.m$ME.pink and ME.pink.m$Condition 

#        PreDrought
#Drought 4.4e-07   

#P value adjustment method: none 



### black ###

pairwise.wilcox.test(ME.black.m$ME.black, ME.black.m$Condition_Time, p.adjust.method = "BH")
#data:  ME.black.m$ME.black and ME.black.m$Condition_Time 

#                  PreDrought_0hr PreDrought_6hr PreDrought_48hr Drought_0hr Drought_6hr
#PreDrought_6hr  0.329          -              -               -           -          
#  PreDrought_48hr 1.000          0.537          -               -           -          
#  Drought_0hr     0.247          0.690          0.247           -           -          
#  Drought_6hr     0.394          0.931          0.485           0.537       -          
#  Drought_48hr    0.065          0.429          0.180           0.792       0.394      

#P value adjustment method: none 

#Data:  ME.black.m$ME.black and ME.black.m$Condition_Time 

#            PreDrought_0hr PreDrought_6hr PreDrought_48hr Drought_0hr Drought_6hr
#PreDrought_6hr  0.73           -              -               -           -          
#  PreDrought_48hr 1.00           0.73           -               -           -          
#  Drought_0hr     0.73           0.86           0.73            -           -          
#  Drought_6hr     0.73           1.00           0.73            0.73        -          
#  Drought_48hr    0.73           0.73           0.73            0.91        0.73       

#P value adjustment method: BH 

pairwise.wilcox.test(ME.black.m$ME.black, ME.black.m$Condition, p.adjust.method = "none")
#data:  ME.black.m$ME.black and ME.black.m$Condition 

#        PreDrought
#Drought 0.073   

#P value adjustment method: none 


####### Kruskal Wallace Dunn test################
#separate into pre-drought and drought, then run tests on each
#magenta
ME.magenta.m.PD <- ME.magenta.m %>%
  filter(Condition == "PreDrought")

ME.magenta.m.D <- ME.magenta.m %>%
  filter(Condition == "Drought")

DunnTest(ME.magenta ~ Time, data = ME.magenta.m.PD, method = "BH") #not sig
DunnTest(ME.magenta ~ Time, data = ME.magenta.m.D, method = "BH") #not sig

#pink
ME.pink.m.PD <- ME.pink.m %>%
  filter(Condition == "PreDrought")

ME.pink.m.D <- ME.pink.m %>%
  filter(Condition == "Drought")

DunnTest(ME.pink ~ Time, data = ME.pink.m.PD, method = "BH") #not sig
DunnTest(ME.pink ~ Time, data = ME.pink.m.D, method = "BH") #0 vs 6hr, p=0.0250

#brown
ME.brown.m.PD <- ME.brown.m %>%
  filter(Condition == "PreDrought")

ME.brown.m.D <- ME.brown.m %>%
  filter(Condition == "Drought")

DunnTest(ME.brown ~ Time, data = ME.brown.m.PD, method = "BH") #0 vs 6hr, p=0.0222; 0 vs 48hr, p=0.0153
DunnTest(ME.brown ~ Time, data = ME.brown.m.D, method = "BH") #not sig

#green
ME.green.m.PD <- ME.green.m %>%
  filter(Condition == "PreDrought")

ME.green.m.D <- ME.green.m %>%
  filter(Condition == "Drought")

DunnTest(ME.green ~ Time, data = ME.green.m.PD, method = "BH") #not sig
DunnTest(ME.green ~ Time, data = ME.green.m.D, method = "BH") #not sig
