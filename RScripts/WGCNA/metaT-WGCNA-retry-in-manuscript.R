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

setwd("~/Documents/R/soil-pyruvate-metagenomics/deseq2/more_samples/only-missing-one-sample/WGCNA/power8-min30/only-black-red-brown")

options(stringsAsFactors = FALSE)

#load metaT gene copy data and column metadata
datExpr0 <-read.csv("../Drought_vs_predrought_metaT_missing_1_no_feature.csv", header = TRUE)
rownames(datExpr0) <- datExpr0$Feature
datExpr0$Feature <- NULL

coldata <- read.csv("../metaT-metadata.csv", header = TRUE)
rownames(coldata) <- coldata$SampleID
coldata$SampleID <- NULL
row.names(coldata) [7] <- "X3300045460_Site1_Drought_48hr"

#filter features that are not in at least x number of samples
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
write.csv(vst.export, "vst-metaT.csv")


datExpr <- t(vst.export)
names(datExpr) <- rownames(datExpr.f)
row.names(datExpr) [7] <- "X3300045460_Site1_Drought_48hr"


head(datExpr, n=5)

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

sampleTree = hclust(dist(datExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
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
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


traitData = read.csv("../metaT-metadata-for-WGCNA.csv", header = TRUE)
traitData$Acetone.C2 <- as.numeric(traitData$Acetone.C2)
traitData$Acetate.C1 <- NULL
traitData$Acetate.C2 <- as.numeric(traitData$Acetate.C2)
traitData$Diacetyl.C2 <- as.numeric(traitData$Diacetyl.C2)
traitData$Time <- NULL
traitData$NOSC <- NULL
traitData$Amino.sugar <- NULL
traitData$Condensed.hydrocarbon <- NULL
traitData$Lignin <- NULL
traitData$Protein <- NULL
traitData$Tannin <- NULL
traitData$Unsaturated.hydrocarbon <- NULL
traitData$Carbohydrate <- NULL
traitData$Lipid <- NULL


dim(traitData)
names(traitData)
#traitData[24:46] <- NULL

# Form a data frame analogous to expression data that will hold the clinical traits.
Samples = rownames(datExpr)
traitRows = match(Samples, traitData$SampleID)
datTraits = traitData[traitRows, -1]
rownames(datTraits) = traitData[traitRows, 1]
collectGarbage()



#datTraits = as.numeric(datTraits$Sex, datTraits$Esonophls_per, datTraits$Neutrophils_per,
#                       datTraits$Macrophages_per, datTraits$Lymphocytes_per,
#                        datTraits$Eosinophils_total, datTraits$Neutrophils_tot, 
#                        datTraits$Lymphocytes_total, datTraits$Treatment, datTraits$Treatment,
#                        datTraits$Hutt_Amish)

#charData = apply(as.matrix(datTraits), 2, as.character)
#numData = apply(charData, 2, as.numeric)

# Re-cluster samples
sampleTree2 = flashClust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = TRUE, naColor = "grey");
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")


save(datExpr, datTraits, file = "germ-free-input.RData")

#s2

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up # certain calculations
# and is optional.
# Any error here may be ignored but you may want to update WGCNA if you see one.
allowWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "germ-free-input.RData");
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

#set integer values to numeric
#datExpr <- lapply(datExpr, as.numeric)
class(datExpr)

#standard "gene" screening based on marginal correlation
#GS1 = as.numeric(cor(y, datExpr, use="p"))

#datExpr = as.numeric(datExpr)
net = blockwiseModules(datExpr, checkMissingData = TRUE, power = 6, minModuleSize = 80 ,
                       numericLabels = TRUE,
                       pamRespectsDendro = FALSE,
                       corTypr = "bicor",
                       maxPOutliers = 0.1,
                       networkType = "signed",
                       #loadTOM = TRUE,
                       TOMType = "signed",
                       saveTOMs = TRUE,
                       saveTOMFileBase = "SoilPyruvateTOM", 
                       minKMEtoStay = 0.35,
                       verbose = 3)
table(net$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "germfree-02-networkConstruction-auto.RData")

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
#MEs_order = MEs[,c(9,8,2,4,5,1)]
MEs_order <- MEs

moduleTraitCor = cor(MEs_order, datTraits, use = "p");

#calculate correlations and pvalues, create each set of p values as a vector to caclulate fdr-corrected q values
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
moduleTraitPvector= c(moduleTraitPvalue)
moduleTraitQobj = qvalue(moduleTraitPvector, lambda = 0,  fdr.level = 0.05)

moduleTraitQvector= moduleTraitQobj$qvalues


ncol(MEs_order)
nrow(MEs_order)

moduleTraitQvalues = matrix(moduleTraitQvector,
                            nrow=9, ncol=4, byrow=FALSE, 
                            dimnames=list(names(MEs_order),
                                          names(datTraits)))




#colnames(MEs) <- sub("")

#sizeGrWindow(20,20)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitQvalues, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));


#change order of rows


#Display the correlation values within a heatmap plot
filename=paste("modules-trait-relationships_all.png",sep="")
png(filename ,width=12, height=20,unit = "in", res = 1000)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = c("Condition", "Acetate-c2", "Acetone-c2", "Diacetyl-c2"),
               #yLabels = c("BLL", "GRN", "RED", "BLU", "BRN", "YLW"),
               yLabels = names(MEs_order),
               ySymbols = names(MEs_order),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               #setStdMargins = FALSE,
               cex.text = 1.0,
               cex.lab.y = 1.2,
               cex.lab.x = 1.2,
               xLabelsAngle = 0,
               xLabelsAdj = 0.5,
               xLabelsPosition = "bottom",
               zlim = c(-1,1))
dev.off()

######Subset modules#####
MEs_order = MEs[,c(9,8,7,1)]

#Remove diacetyl for figure
datTraits$Diacetyl.C2 <- NULL

#correlate traits and MEs
moduleTraitCor = cor(MEs_order, datTraits, use = "p");

#calculate correlations and pvalues, create each set of p values as a vector to caclulate fdr-corrected q values
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
moduleTraitPvector= c(moduleTraitPvalue)
moduleTraitQobj = qvalue(moduleTraitPvector, lambda = 0,  fdr.level = 0.05)

moduleTraitQvector= moduleTraitQobj$qvalues


ncol(MEs_order)
nrow(MEs_order)

moduleTraitQvalues = matrix(moduleTraitQvector,
                            nrow=4, ncol=3, byrow=FALSE, 
                            dimnames=list(names(MEs_order),
                                          names(datTraits)))




#colnames(MEs) <- sub("")

#sizeGrWindow(20,20)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitQvalues, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));


#change order of rows


#Display the correlation values within a heatmap plot
filename=paste("modules-trait-relationships_subset.png",sep="")
png(filename ,width=6, height=6,unit = "in", res = 1000)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = c("Condition", "Acetate-C2", "Acetone-C2"),
               yLabels = c("Pink", "Green", "Magenta", "Brown"),
               #yLabels = names(MEs_order),
               ySymbols = names(MEs_order),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               #setStdMargins = FALSE,
               cex.text = 1,
               cex.lab.y = 1,
               cex.lab.x = 1,
               xLabelsAngle = 0,
               xLabelsAdj = 0.6,
               xLabelsPosition = "bottom",
               zlim = c(-1,1))
dev.off()

# Plot the relationships among the eigengenes 
filename=paste("modules_eigengene3-no-pink-turq.pdf",sep="")
pdf(filename, width=8, height=6)
plotEigengeneNetworks(MEs_order,"",marDendro=c(0,4,1,6), 
                      marHeatmap=c(3,4,1,2), cex.lab=0.8,xLabelsAngle=90)
dev.off()

# Plot the relationships among the eigengenes as png
filename=paste("modules_eigengene3-no-pink-turq.png",sep="")
png(filename, width=12, height=9, unit = "in", res = 1000)
plotEigengeneNetworks(MEs_order,"",marDendro=c(0,4,1,6), 
                      marHeatmap=c(3,4,1,2), cex.lab=0.8,xLabelsAngle=90)
dev.off()


#####Condition#####
#Define variable MT.ND1 containing the MT.ND1 column of datTrait
Condition = as.data.frame(datTraits$Condition);
names(Condition) = ("Condition")
# names (colors) of the modules
modNames = substring(names(MEs_order), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs_order, use = "p"));
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

names(geneTraitSignificance) = paste("GS.", names(Condition), sep="");
names(GSQvalue) = paste("q.GS.", names(Condition), sep="");
module = "black"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Condition",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#write.csv(geneModuleMembership, file= "geneModuleMembeship_updated_trait_data-q.csv")
write.csv(geneTraitSignificance, file= "geneTraitsSignif-Condition-black-updated-trait-data-q-no-pink-turq.csv")

#Summary output of network analysis results
names(datExpr) <- rownames(datExpr.f)
x <- names(datExpr)[moduleColors=="black"]
write.csv(x, file = "datExpr_black-no-pink-tuq.csv")

moduleColor=names(MEs_order)
probes=names(datExpr)
# Create the starting data frame
geneInfo0 = data.frame(datExpr = probes, ModuleColor = moduleColors,
                       geneTraitSignificance,
                       GSQvalue)

# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs_order, Condition, use = "p")));

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
geneOrder = order(geneInfo0$ModuleColor, -abs(geneInfo0$GS.Condition));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo_Condition_no-pink-tur.csv")
write.csv(geneInfo0, file = "geneInfo-not-Condition-ordered-q-no-pink-turq.csv")

#gene screening method bsed on detailied definition of module membership
NS1=networkScreening(y=datTraits$Condition, datME=MEs0, datExpr=datExpr,
                     oddPower=3, blockSize=1000, minimumSampleSize=4,
                     addMEy=TRUE, removeDiag=FALSE, weightESy=0.5)
write.csv(NS1, "condition-gene-screening-no-pink-turq.csv")


#####Acetate-c2#####
#Define variable AcetateC2 containing the Acetate.C2 column of datTrait
AcetateC2 = as.data.frame(datTraits$Acetate.C2);
names(AcetateC2) = ("AcetateC2")
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
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

names(geneTraitSignificance) = paste("GS.", names(AcetateC2), sep="");
names(GSQvalue) = paste("q.GS.", names(AcetateC2), sep="");
module = "green"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for AcetateC2",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#write.csv(geneModuleMembership, file= "geneModuleMembeship_updated_trait_data-q.csv")
write.csv(geneTraitSignificance, file= "geneTraitsSignif-AcetateC2-green-updated-trait-data-q.csv")

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
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod
                                                                    ]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$ModuleColor, -abs(geneInfo0$GS.AcetateC2));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo_Acetatec2.csv")
write.csv(geneInfo0, file = "geneInfo-not-Acetatec2-ordered-q.csv")

#gene screening method bsed on detailied definition of module membership
NS1=networkScreening(y=datTraits$Acetate.C2, datME=MEs0, datExpr=datExpr,
                     oddPower=3, blockSize=1000, minimumSampleSize=4,
                     addMEy=TRUE, removeDiag=FALSE, weightESy=0.5)
write.csv(NS1, "gene-screening-acetatec2.csv")


#####acetone-c2#####
#Define variable AcetateC2 containing the Acetate.C2 column of datTrait
AcetoneC2 = as.data.frame(datTraits$Acetone.C2);
names(AcetoneC2) = ("AcetoneC2")
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
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

names(geneTraitSignificance) = paste("GS.", names(AcetoneC2), sep="");
names(GSQvalue) = paste("q.GS.", names(AcetoneC2), sep="");
module = "green"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for AcetoneC2",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#write.csv(geneModuleMembership, file= "geneModuleMembeship_updated_trait_data-q.csv")
write.csv(geneTraitSignificance, file= "geneTraitsSignif-AcetoneC2-green-updated-trait-data-q.csv")

# Create the starting data frame
geneInfo0 = data.frame(datExpr = probes, ModuleColor = moduleColors,
                       geneTraitSignificance,
                       GSQvalue)

# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, AcetoneC2, use = "p")));

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
geneOrder = order(geneInfo0$ModuleColor, -abs(geneInfo0$GS.Diacetyl.C2));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo_AcetoneC2.csv")
write.csv(geneInfo0, file = "geneInfo-not-AcetoneC2-ordered-q.csv")

#gene screening method bsed on detailied definition of module membership
NS1=networkScreening(y=datTraits$Acetone.C2, datME=MEs0, datExpr=datExpr,
                     oddPower=3, blockSize=1000, minimumSampleSize=4,
                     addMEy=TRUE, removeDiag=FALSE, weightESy=0.5)
write.csv(NS1, "gene-screening-AcetoneC2.csv")

#####Diacetyl-c2#####
#Define variable AcetateC2 containing the Acetate.C2 column of datTrait
DiacetylC2 = as.data.frame(datTraits$Diacetyl.C2);
names(DiacetylC2) = ("DiacetylC2")
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, DiacetylC2, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
#GSPvalue = corPvalueStudent(as.matrix(geneTraitSignificance), nSamples)
GSPvalue.0 = as.vector(GSPvalue[,'DiacetylC2'])
GSQvalue.0 = qvalue(p=GSPvalue.0, lambda = 0, fdr.level = 0.05)
GSQvalue <- as.data.frame(GSQvalue.0$qvalues, row.names=names(datExpr))
summary(GSQvalue.0)


pi0 <- GSQvalue.0$pi0
lfdr <- GSQvalue.0$lfdr 

names(geneTraitSignificance) = paste("GS.", names(DiacetylC2), sep="");
names(GSQvalue) = paste("q.GS.", names(DiacetylC2), sep="");
module = "green"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for DiacetylC2",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#write.csv(geneModuleMembership, file= "geneModuleMembeship_updated_trait_data-q.csv")
write.csv(geneTraitSignificance, file= "geneTraitsSignif-DiacetylC2-green-updated-trait-data-q.csv")

# Create the starting data frame
geneInfo0 = data.frame(datExpr = probes, ModuleColor = moduleColors,
                       geneTraitSignificance,
                       GSQvalue)

# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, DiacetylC2, use = "p")));

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
geneOrder = order(geneInfo0$ModuleColor, -abs(geneInfo0$GS.DiacetylC2));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo_DiacetylC2.csv")
write.csv(geneInfo0, file = "geneInfo-not-DiacetylC2-ordered-q.csv")

#gene screening method bsed on detailied definition of module membership
NS1=networkScreening(y=datTraits$Acetone.C2, datME=MEs0, datExpr=datExpr,
                     oddPower=3, blockSize=1000, minimumSampleSize=4,
                     addMEy=TRUE, removeDiag=FALSE, weightESy=0.5)
write.csv(NS1, "gene-screening-DiacetylC2.csv")

####calculate intramodular connectivity
ADJ1=abs(cor(datExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ1,moduleColors)
head(Alldegrees1)
write.csv(Alldegrees1, file = "module-interconnectivity.csv")

#calculate gene significance Il13_slope
GS_cond=as.numeric(cor(datTraits$Condition,datExpr,use="p"))
geneSignificance=abs(GS_cond)

GS_acetatec2 = as.numeric(cor(datTraits$Acetate.C2, datExpr, use="p"))
geneSignificance = abs(GS_acetatec2)

GS_Prg3 = as.numeric(cor(datTraits$Prg3, datExpr, use="p"))
geneSignificance = abs(GS_Prg3)

GS_reos = as.numeric(cor(datTraits$Resident_Eos, datExpr, use="p"))
geneSignificance = abs(GS_reos)

GS_reos = as.numeric(cor(datTraits$Resident_Eos, datExpr, use="p"))
geneSignificance = abs(GS_reos)

moduleSignificance=tapply(geneSignificance, moduleColors, mean, na.rm=T)

##relationship between gene significance and intramodular connectivity
colorlevels=unique(moduleColors)
sizeGrWindow(9,6)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))

for (i in c(1:length(colorlevels)))
{
  whichmodule=colorlevels[[i]];
  restrict1 = (moduleColors==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     geneSignificance[restrict1], col=moduleColors[restrict1],
                     main=whichmodule,
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}

##relationship between gene significance and module membership
ModMem <- read.csv("geneInfo-not-Condition-ordered-q.csv", header = TRUE)
colorlevels=unique(moduleColors)
sizeGrWindow(9,6)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))

for (i in c(1:length(colorlevels)))
{
  whichmodule=colorlevels[[i]];
  restrict1 = (moduleColors==whichmodule);
  verboseScatterplot(ModMem$MM[restrict1],
                     geneSignificance[restrict1], col=moduleColors[restrict1],
                     main=whichmodule,
                     xlab = "Module membership", ylab = "Gene Significance", abline = TRUE)
}


datKME=signedKME(datExpr, MEs0, outputColumnName="MM.")
head(datKME)


#Find genes with high gene significance and high intramodular connectivity in ineresting modules
FilterGenes = abs(GS_cond)> .5 & datKME$MM.green> .8
table(FilterGenes)



x <- dimnames(data.frame(datExpr))[[2]][FilterGenes]
write.csv(x,file="filtered_cond_green.csv")


#relationship between module membership and intramodular connectivity
sizeGrWindow(8,6)
par(mfrow=c(2,2))
# We choose 4 modules to plot: turquoise, blue, brown, green.
# For simplicity we write the code out explicitly for each module.
which.color="green";
restrictGenes=moduleColors==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                   col=which.color,
                   xlab="Intramodular Connectivity",
                   ylab="(Module Membership)^6")
which.color="black";
restrictGenes=moduleColors==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                   col=which.color,
                   xlab="Intramodular Connectivity",
                   ylab="(Module Membership)^6")
which.color="yellow";
restrictGenes=moduleColors==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                   col=which.color,
                   xlab="Intramodular Connectivity",
                   ylab="(Module Membership)^6")

which.color="brown";
restrictGenes=moduleColors==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                   col=which.color,
                   xlab="Intramodular Connectivity",
                   ylab="(Module Membership)^6")



##plot eigengene expression per condition_time
coldata$Condition_Time <- factor(coldata$Condition_Time, c("PreDrought_0hr", "PreDrought_6hr", "PreDrought_48hr",
                                                           "Drought_0hr", "Drought_6hr", "Drought_48hr"))
coldata$Condition <- factor(coldata$Condition, c("PreDrought", "Drought"))
coldata$Time <- factor(coldata$Time, c("0hr", "6hr", "48hr"))

coldata$SampleID <- rownames(coldata)

###plot eigengene expression per sample, 
filename=paste("black-EGexp.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
which.module="black"
signif(cor(MEs, use="p"), 2)
ME.black=MEs[, paste("ME",which.module, sep="")]
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME.black, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
dev.off()

filename=paste("red-EGexp.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
which.module="red"
signif(cor(MEs, use="p"), 2)
ME.red=MEs[, paste("ME",which.module, sep="")]
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME.red, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
dev.off()

filename=paste("magenta-EGexp.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
which.module="magenta"
signif(cor(MEs, use="p"), 2)
ME.magenta=MEs[, paste("ME",which.module, sep="")]
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME.magenta, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
dev.off()

filename=paste("blue-EGexp.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
which.module="blue"
signif(cor(MEs, use="p"), 2)
ME.blue=MEs[, paste("ME",which.module, sep="")]
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME.blue, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
dev.off()

filename=paste("yellow-EGexp.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
which.module="yellow"
signif(cor(MEs, use="p"), 2)
ME.yellow=MEs[, paste("ME",which.module, sep="")]
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME.yelow, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
dev.off()

filename=paste("brown-EGexp.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
which.module="brown"
signif(cor(MEs, use="p"), 2)
ME.brown=MEs[, paste("ME",which.module, sep="")]
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME.brown, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
dev.off()

filename=paste("turquoise-EGexp.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
which.module="turquoise"
signif(cor(MEs, use="p"), 2)
ME.turquoise=MEs[, paste("ME",which.module, sep="")]
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME.turquoise, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
dev.off()

filename=paste("pink-EGexp.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
which.module="pink"
signif(cor(MEs, use="p"), 2)
ME.pink=MEs[, paste("ME",which.module, sep="")]
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME.pink, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")
dev.off()

filename=paste("green-EGexp.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
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

filename=paste("black-EGexp-facet.png", sep = "")
png(filename ,width=4, height=2, unit='in', res = 1000)
ggplot(ME.black.m, aes(x=Time, y= ME.black)) +
  geom_point() + geom_boxplot(fill = "darkgrey") +
  facet_grid(~Condition) +
  ylim(-0.35, 0.4) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  labs(title = "black module", y = "Eigengene Expression", x= "Time post pyruvate addition (hr)") +
  theme_bw() +
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"))
dev.off()

#red
ME.red <- as.data.frame(ME.red)
rownames(ME.red) <- rownames(datExpr)
ME.red$SampleID <- rownames(ME.red)
ME.red.m <- ME.red %>%
  merge(.,coldata, by = "SampleID")

filename=paste("red-EGexp-facet.png", sep = "")
png(filename ,width=4, height=2, unit='in', res = 1000)
ggplot(ME.red.m, aes(x=Time, y= ME.red)) +
  geom_point() + geom_boxplot(fill = "red") +
  facet_grid(~Condition) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  labs(title = "red module", y = "Eigengene Expression", x= "Time post pyruvate addition (hr)") +
  theme_bw() +
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"))
dev.off()

#magenta
ME.magenta <- as.data.frame(ME.magenta)
rownames(ME.magenta) <- rownames(datExpr)
ME.magenta$SampleID <- rownames(ME.magenta)
ME.magenta.m <- ME.magenta %>%
  merge(.,coldata, by = "SampleID")

filename=paste("magenta-EGexp-facet.png", sep = "")
png(filename ,width=4, height=2, unit='in', res = 1000)
ggplot(ME.magenta.m, aes(x=Time, y= ME.magenta)) +
  geom_point() + geom_boxplot(fill = "magenta") +
  facet_grid(~Condition) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  labs(title = "magenta module", y = "Eigengene Expression", x= "Time post pyruvate addition (hr)") +
  theme_bw() +
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"))
dev.off()

#blue
ME.blue <- as.data.frame(ME.blue)
rownames(ME.blue) <- rownames(datExpr)
ME.blue$SampleID <- rownames(ME.blue)
ME.blue.m <- ME.blue %>%
  merge(.,coldata, by = "SampleID")

filename=paste("blue-EGexp-facet.png", sep = "")
png(filename ,width=4, height=2, unit='in', res = 1000)
ggplot(ME.blue.m, aes(x=Time, y= ME.blue)) +
  geom_point() + geom_boxplot(fill = "blue") +
  facet_grid(~Condition) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  labs(title = "blue module", y = "Eigengene Expression", x= "Time post pyruvate addition (hr)") +
  theme_bw() +
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"))
dev.off()

#yellow
ME.yellow <- as.data.frame(ME.yellow)
rownames(ME.yellow) <- rownames(datExpr)
ME.yellow$SampleID <- rownames(ME.yellow)
ME.yellow.m <- ME.yellow %>%
  merge(.,coldata, by = "SampleID")

filename=paste("yellow-EGexp-facet.png", sep = "")
png(filename ,width=4, height=2, unit='in', res = 1000)
ggplot(ME.yellow.m, aes(x=Time, y= ME.yellow)) +
  geom_point() + geom_boxplot(fill = "yellow") +
  facet_grid(~Condition) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  labs(title = "yellow module", y = "Eigengene Expression", x= "Time post pyruvate addition (hr)") +
  ylim(-0.3, 0.4) +
  theme_bw() +
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"))
dev.off()

#brown
ME.brown <- as.data.frame(ME.brown)
rownames(ME.brown) <- rownames(datExpr)
ME.brown$SampleID <- rownames(ME.brown)
ME.brown.m <- ME.brown %>%
  merge(.,coldata, by = "SampleID")

filename=paste("brown-EGexp-facet.png", sep = "")
png(filename ,width=4, height=2, unit='in', res = 1000)
ggplot(ME.brown.m, aes(x=Time, y= ME.brown)) +
  geom_point() + geom_boxplot(fill = "brown") +
  facet_grid(~Condition) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  ylim(-.45, .3) +
  labs(title = "brown module", y = "Eigengene Expression", x= "Time post pyruvate addition (hr)") +
  theme_bw() +
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"))
dev.off()

#turquoise
ME.turquoise <- as.data.frame(ME.turquoise)
rownames(ME.turquoise) <- rownames(datExpr)
ME.turquoise$SampleID <- rownames(ME.turquoise)
ME.turquoise.m <- ME.turquoise %>%
  merge(.,coldata, by = "SampleID")

filename=paste("turquoise-EGexp-facet.png", sep = "")
png(filename ,width=4, height=2, unit='in', res = 1000)
ggplot(ME.turquoise.m, aes(x=Time, y= ME.turquoise)) +
  geom_point() + geom_boxplot(fill = "turquoise") +
  facet_grid(~Condition) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  labs(title = "turquoise module", y = "Eigengene Expression", x= "Time post pyruvate addition (hr)") +
  theme_bw() +
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"))
dev.off()

#pink
ME.pink <- as.data.frame(ME.pink)
rownames(ME.pink) <- rownames(datExpr)
ME.pink$SampleID <- rownames(ME.pink)
ME.pink.m <- ME.pink %>%
  merge(.,coldata, by = "SampleID")

filename=paste("pink-EGexp-facet.png", sep = "")
png(filename ,width=4, height=2, unit='in', res = 1000)
ggplot(ME.pink.m, aes(x=Time, y= ME.pink)) +
  geom_point() + geom_boxplot(fill = "pink") +
  facet_grid(~Condition) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  labs(title = "pink module", y = "Eigengene Expression", x= "Time post pyruvate addition (hr)") +
  theme_bw() +
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"))
dev.off()

#green
ME.green <- as.data.frame(ME.green)
rownames(ME.green) <- rownames(datExpr)
ME.green$SampleID <- rownames(ME.green)
ME.green.m <- ME.green %>%
  merge(.,coldata, by = "SampleID")

filename=paste("green-EGexp-facet.png", sep = "")
png(filename ,width=4, height=2, unit='in', res = 1000)
ggplot(ME.green.m, aes(x=Time, y= ME.green)) +
  geom_point() + geom_boxplot(fill = "green") +
  facet_grid(~Condition) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  labs(title = "green module", y = "Eigengene Expression", x= "Time post pyruvate addition (hr)") +
  theme_bw() +
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"))
dev.off()

##statistics on eigenggene expression
pairwise.wilcox.test(ME.pink.m$ME.pink, ME.pink.m$Condition_Time, p.adjust.method = "none")
#data:  ME.pink.m$ME.pink and ME.pink.m$Condition_Time 

#              PreDrought_0hr PreDrought_6hr PreDrought_48hr Drought_0hr Drought_6hr
#PreDrought_6hr  0.0303         -              -               -           -          
#  PreDrought_48hr 0.0931         0.9307         -               -           -          
#  Drought_0hr     0.0043         0.0079         0.0087          -           -          
#  Drought_6hr     0.0043         0.0303         0.1320          0.0173      -          
#  Drought_48hr    0.0022         0.0043         0.0043          0.1255      0.1797     

#P value adjustment method: none 

pairwise.wilcox.test(ME.pink.m$ME.pink, ME.pink.m$Condition, p.adjust.method = "none")
#data:  ME.pink.m$ME.pink and ME.pink.m$Condition 

#        PreDrought
#Drought 4.4e-07   

#P value adjustment method: none 

pairwise.wilcox.test(ME.black.m$ME.black, ME.black.m$Condition_Time, p.adjust.method = "none")
#data:  ME.black.m$ME.black and ME.black.m$Condition_Time 

#                  PreDrought_0hr PreDrought_6hr PreDrought_48hr Drought_0hr Drought_6hr
#PreDrought_6hr  0.329          -              -               -           -          
#  PreDrought_48hr 1.000          0.537          -               -           -          
#  Drought_0hr     0.247          0.690          0.247           -           -          
#  Drought_6hr     0.394          0.931          0.485           0.537       -          
#  Drought_48hr    0.065          0.429          0.180           0.792       0.394      

#P value adjustment method: none 

pairwise.wilcox.test(ME.black.m$ME.black, ME.black.m$Condition, p.adjust.method = "none")
#data:  ME.black.m$ME.black and ME.black.m$Condition 

#        PreDrought
#Drought 0.073   

#P value adjustment method: none 

pairwise.wilcox.test(ME.red.m$ME.red, ME.red.m$Condition_Time, p.adjust.method = "none")
#data:  ME.red.m$ME.red and ME.red.m$Condition_Time 

#PreDrought_0hr PreDrought_6hr PreDrought_48hr Drought_0hr Drought_6hr
#PreDrought_6hr  0.0519         -              -               -           -          
#PreDrought_48hr 0.1320         0.7922         -               -           -          
#  Drought_0hr     0.9307         0.1508         0.4286          -           -          
#  Drought_6hr     0.1797         0.4286         0.6991          0.3290      -          
#  Drought_48hr    0.0087         0.1255         0.2403          0.1255      0.1797     

#P value adjustment method: none 

pairwise.wilcox.test(ME.red.m$ME.red, ME.red.m$Condition, p.adjust.method = "none")
#data:  ME.red.m$ME.red and ME.red.m$Condition 

#        PreDrought
#Drought 0.45     

#P value adjustment method: none 


pairwise.wilcox.test(ME.brown.m$ME.brown, ME.brown.m$Condition_Time, p.adjust.method = "none")
#data:  ME.brown.m$ME.brown and ME.brown.m$Condition_Time 

#               PreDrought_0hr PreDrought_6hr PreDrought_48hr Drought_0hr Drought_6hr
#PreDrought_6hr  0.0043         -              -               -           -          
#  PreDrought_48hr 0.0087         0.7922         -               -           -          
#  Drought_0hr     0.0043         0.0079         0.1775          -           -          
#  Drought_6hr     0.0022         0.0173         0.3095          0.5368      -          
#  Drought_48hr    0.0022         0.0173         0.1797          0.6623      0.5887     

#P value adjustment method: none

pairwise.wilcox.test(ME.brown.m$ME.brown, ME.brown.m$Condition, p.adjust.method = "none")
#data:  ME.brown.m$ME.brown and ME.brown.m$Condition 

#      PreDrought
#Drought 2.8e-05   

#P value adjustment method: none 

pairwise.wilcox.test(ME.green.m$ME.green, ME.green.m$Condition_Time, p.adjust.method = "none")
#data:  ME.green.m$ME.green and ME.green.m$Condition_Time 

#               PreDrought_0hr PreDrought_6hr PreDrought_48hr Drought_0hr Drought_6hr
#PreDrought_6hr  0.052          -              -               -           -          
#  PreDrought_48hr 0.310          0.177          -               -           -          
##  Drought_0hr     0.126          0.548          0.177           -           -          
#  Drought_6hr     0.026          0.429          0.065           0.792       -          
#  Drought_48hr    0.015          0.082          0.041           0.082       0.093      

#P value adjustment method: none 

pairwise.wilcox.test(ME.green.m$ME.green, ME.green.m$Condition, p.adjust.method = "none")
#data: ME.green.m$ME.green and ME.green.m$Condition

#      PreDrought
#Drought 0.0015   

#P value adjustment method: none 

pairwise.wilcox.test(ME.blue.m$ME.blue, ME.blue.m$Condition_Time, p.adjust.method = "none")
#data:  ME.blue.m$ME.blue and ME.blue.m$Condition_Time 

#                PreDrought_0hr PreDrought_6hr PreDrought_48hr Drought_0hr Drought_6hr
#PreDrought_6hr    0.0519         -              -               -           -          
 # PreDrought_48hr 0.5887         0.3290         -               -           -          
#  Drought_0hr     0.5368         0.0079         0.3290          -           -          
#  Drought_6hr     0.6991         0.0303         0.6991          0.3290      -          
#  Drought_48hr    0.9372         0.0087         0.5887          0.5368      0.5887     

#P value adjustment method: none 

pairwise.wilcox.test(ME.blue.m$ME.blue, ME.blue.m$Condition, p.adjust.method = "none")
#data:  ME.blue.m$ME.blue and ME.blue.m$Condition 

#        PreDrought
#Drought 0.079   

#P value adjustment method: none 

pairwise.wilcox.test(ME.magenta.m$ME.magenta, ME.magenta.m$Condition_Time, p.adjust.method = "none")
#data:  ME.magenta.m$ME.magenta and ME.magenta.m$Condition_Time 

#               PreDrought_0hr PreDrought_6hr PreDrought_48hr Drought_0hr Drought_6hr
#PreDrought_6hr  0.537          -              -               -           -          
#  PreDrought_48hr 0.041          0.177          -               -           -          
#  Drought_0hr     0.177          0.095          0.052           -           -          
#  Drought_6hr     0.699          0.329          0.065           0.247       -          
#  Drought_48hr    0.699          0.931          0.093           0.126       0.589      

#P value adjustment method: none 

pairwise.wilcox.test(ME.magenta.m$ME.magenta, ME.magenta.m$Condition, p.adjust.method = "none")
#data:  ME.magenta.m$ME.magenta and ME.magenta.m$Condition 

#PreDrought
#Drought 0.049     

#P value adjustment method: none 

## Kruskal Wallace Dunn test
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


