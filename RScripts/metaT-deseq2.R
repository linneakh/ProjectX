###Linnea Honeker
###9/14/21
###Analyze metaT data for differential gene abundance. Raw input is gene copy numbers obtained from IMG.


library("DESeq2")
packageVersion("DESeq2")
library(dplyr)
library(knitr)
library(apeglm)
library("ggplot2")
library("pheatmap")
library("stringr")

setwd("~/Documents/R/soil-pyruvate-metagenomics/deseq2/more_samples/only-missing-one-sample")

#import count data, transform to matrix, import column metadata, check that its a data.frame
cts.0 <- read.csv("Drought_vs_predrought_metaT_missing_1_no_feature.csv", header = TRUE, sep =",")
coldata <- read.csv("metaT-metadata.csv", header = TRUE, sep = ",")
cts.02 <- as.matrix(cts.0)
row.names(cts.02) <- cts.02[,1]
cts.02 <- cts.02[,-1]

cts <- matrix(as.numeric(cts.02),    # Convert to numeric matrix
                  ncol = ncol(cts.02))
cts.0$Feature = NULL
colnames(cts) <- colnames(cts.0)

row.names(cts) <- row.names(cts.02)

row.names(coldata) <- coldata$SampleID
coldata$SampleID <- NULL


###create a deseq2 object with Condition, Time, and Site as design###
coldata$Condition <- factor(coldata$Condition, c("PreDrought", "Drought"))

coldata$Condition <- relevel(coldata$Condition, ref = "PreDrought")


dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ Condition + Site)

###data exploration###
vstcounts <- vst(dds, blind=TRUE)
vst.export<- assay(vstcounts)
vst.export<- as.data.frame(vst.export)
write.csv(vst.export, "vst-metaT.csv")

plotPCA(vstcounts, intgroup=c("Condition"))

select <- order(rowMeans(counts(dds)),
                decreasing=TRUE)[1:500]
df <- as.data.frame(colData(dds)[,c("Condition","Site")])

filename=paste("figs/heatmap-500-vst.png", sep = "")
png(filename ,width=10, height=10, unit='in', res = 1000)
pheatmap(assay(vstcounts)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
dev.off()

###deifferential expression analysis###

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

##drought vs predrought time 0
res <- results(dds, contrast = c("Condition","Drought","PreDrought"))
plotMA(res, ylim=c(-40,40))

summary(res)

#make volcano plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-40,40), ylim=c(0,40)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#shrink effect size
resLFC <- lfcShrink(dds, coef="Condition_Drought_vs_PreDrought", type="apeglm")
resLFC
plotMA(resLFC, ylim=c(-2,2))

summary(resLFC)

#make volcano plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(resLFC, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-15,10), ylim=c(0,25)))
with(subset(resLFC, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(resLFC, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#export results
write.csv(as.data.frame(res), 
          file="drought-vs-predrought.csv")
write.csv(as.data.frame(resLFC),
          file="drought-vs-predrought-lfc.csv")

#cmopare results with WGCNA 
resLFC.e <- as.data.frame(resLFC)
resLFC.e$datExpr <- rownames(resLFC.e)

T0.WGCNA.lfc.dp <- geneInfo0 %>%
  filter(datExpr != "<NA>") %>%
  merge(.,resLFC.e, by = "datExpr", all = TRUE) %>%
  arrange(pvalue)

write.csv(T0.WGCNA.lfc.dp, "T0.WGCNA.dp.lfc.csv")

res.e <- as.data.frame(res)
res.e$datExpr <- rownames(res.e)

T0.WGCNA.dp <- geneInfo0 %>%
  filter(datExpr != "<NA>") %>%
  merge(.,res.e, by = "datExpr", all = TRUE) %>%
  arrange(pvalue)

write.csv(T0.WGCNA.dp, "T0.WGCNA.dp.csv")

###create a deseq2 object with Condition, Time, and Site as design###
coldata$Condition_Time <- factor(coldata$Condition_Time, c("PreDrought_0hr", "PreDrought_6hr", "PreDrought_48hr",
                                                           "Drought_0hr", "Drought_6hr", "Drought_48hr"))

coldata$Condition_Time <- relevel(coldata$Condition_Time, ref = "PreDrought_0hr")


dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ Condition_Time + Site)

###data exploration###
vstcounts <- vst(dds, blind=TRUE)
vst.export<- assay(vstcounts)
vst.export<- as.data.frame(vst.export)
write.csv(vst.export, "vst-metaT.csv")

plotPCA(vstcounts, intgroup=c("Condition"))

select <- order(rowMeans(counts(dds)),
                decreasing=TRUE)[1:500]
df <- as.data.frame(colData(dds)[,c("Condition_Time","Site")])

filename=paste("figs/heatmap-500-vst.png", sep = "")
png(filename ,width=10, height=10, unit='in', res = 1000)
pheatmap(assay(vstcounts)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
dev.off()

###deifferential expression analysis###

dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

##drought vs predrought time 0
res.T0 <- results(dds, contrast = c("Condition_Time","Drought_0hr","PreDrought_0hr"))
plotMA(res, ylim=c(-40,40))

summary(res.T0)

#make volcano plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res.T0, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-40,40), ylim=c(0,40)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res.T0, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res.T0, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#shrink effect size
resLFC.T0 <- lfcShrink(dds, coef="Condition_Time_Drought_0hr_vs_PreDrought_0hr", type="apeglm")
resLFC.T0
plotMA(resLFC.T0, ylim=c(-2,2))

summary(resLFC.T0)

#make volcano plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(resLFC.T0, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-15,10), ylim=c(0,25)))
with(subset(resLFC.T0, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(resLFC.T0, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#export results
write.csv(as.data.frame(res.T0), 
          file="drought-0hr-vs-predrought-0hr.csv")
write.csv(as.data.frame(resLFC.T0),
          file="drought-0hr-vs-predrought-0hr-lfc.csv")

#cmopare results with WGCNA 
resLFC.T0e <- as.data.frame(resLFC.T0)
resLFC.T0e$datExpr <- rownames(resLFC.T0e)

T0.LFC.WGCNA <- geneInfo0 %>%
  filter(datExpr != "<NA>") %>%
  merge(.,resLFC.T0e, by = "datExpr", all = TRUE) %>%
  arrange(pvalue)

write.csv(T0.LFC.WGCNA, "T0.LFC.WGCNA.csv")

res.T0e <- as.data.frame(res.T0)
res.T0e$datExpr <- rownames(res.T0e)

T0.WGCNA <- geneInfo0 %>%
  filter(datExpr != "<NA>") %>%
  merge(.,res.T0e, by = "datExpr", all = TRUE) %>%
  arrange(pvalue)

write.csv(T0.WGCNA, "T0.WGCNA.csv")

##Predrought time 6 vs predrought time 0
res.P6.0 <- results(dds, contrast = c("Condition_Time","PreDrought_6hr","PreDrought_0hr"))
plotMA(res, ylim=c(-40,40))

summary(res.P6.0)

#make volcano plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res.P6.0, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-40,40), ylim=c(0,40)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res.P6.0, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res.P6.0, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#shrink effect size
resLFC.P6.0 <- lfcShrink(dds, coef="Condition_Time_PreDrought_6hr_vs_PreDrought_0hr", type="apeglm")
resLFC.P6.0
plotMA(resLFC.P6.0, ylim=c(-2,2))

summary(resLFC)

#make volcano plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(resLFC.P6.0, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-15,10), ylim=c(0,25)))
with(subset(resLFC.P6.0, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(resLFC.P6.0, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#export results
write.csv(as.data.frame(res.P6.0), 
          file="pre-drought-6hr-vs-predrought-0hr.csv")
write.csv(as.data.frame(resLFC.P6.0),
          file="pre-drought-6hr-vs-predrought-0hr-lfc.csv")

#compare DEG with xx module


##Predrought time 48 vs predrought time 0
res <- results(dds, contrast = c("Condition_Time","PreDrought_48hr","PreDrought_0hr"))
plotMA(res, ylim=c(-40,40))

summary(res)

#make volcano plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-40,40), ylim=c(0,40)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#shrink effect size
resLFC <- lfcShrink(dds, coef="Condition_Time_PreDrought_48hr_vs_PreDrought_0hr", type="apeglm")
resLFC
plotMA(resLFC, ylim=c(-2,2))

summary(resLFC)

#make volcano plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(resLFC, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-15,10), ylim=c(0,25)))
with(subset(resLFC, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(resLFC, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#export results
write.csv(as.data.frame(res), 
          file="pre-drought-48hr-vs-predrought-0hr.csv")
write.csv(as.data.frame(resLFC),
          file="pre-drought-48hr-vs-predrought-0hr-lfc.csv")

###drought time 6 vs drought time 0#####
###create a deseq2 object with Condition, Time, and Site as design###
coldata$Condition_Time <- relevel(coldata$Condition_Time, ref = "Drought_0hr")


dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ Condition_Time + Site)

###deifferential expression analysis###
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

res.d6 <- results(dds, contrast = c("Condition_Time","Drought_6hr","Drought_0hr"))
plotMA(res, ylim=c(-40,40))

summary(res.d6)

#make volcano plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res.d6, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-40,40), ylim=c(0,40)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res.d6, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res.d6, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#shrink effect size
resLFC.d6 <- lfcShrink(dds, coef="Condition_Time_Drought_6hr_vs_Drought_0hr", type="apeglm")
resLFC.d6
plotMA(resLFC.d6, ylim=c(-2,2))

summary(resLFC.d6)

#make volcano plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(resLFC.d6, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-15,10), ylim=c(0,25)))
with(subset(resLFC.d6, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(resLFC.d6, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#export results
write.csv(as.data.frame(res.d6), 
          file="drought-6hr-vs-drought-0hr.csv")
write.csv(as.data.frame(resLFC.d6),
          file="drought-6hr-vs-drought-0hr-lfc.csv")

#cmopare results with WGCNA 
resLFC.d6e <- as.data.frame(resLFC.d6)
resLFC.d6e$datExpr <- rownames(resLFC.d6e)

T0.WGCNA.d6 <- geneInfo0 %>%
  filter(datExpr != "<NA>") %>%
  merge(.,resLFC.d6e, by = "datExpr", all = TRUE) %>%
  arrange(pvalue)

write.csv(T0.WGCNA.d6, "T0.WGCNA.d6.csv")

##drought time 48 vs drought time 0##
res <- results(dds, contrast = c("Condition_Time","Drought_48hr","Drought_0hr"))
plotMA(res, ylim=c(-40,40))

summary(res)

#make volcano plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-40,40), ylim=c(0,40)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#shrink effect size
resLFC <- lfcShrink(dds, coef="Condition_Time_Drought_48hr_vs_Drought_0hr", type="apeglm")
resLFC
plotMA(resLFC, ylim=c(-2,2))

summary(resLFC)

#make volcano plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(resLFC, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-15,10), ylim=c(0,25)))
with(subset(resLFC, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(resLFC, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#export results
write.csv(as.data.frame(res), 
          file="drought-48hr-vs-drought-0hr.csv")
write.csv(as.data.frame(resLFC),
          file="drought-48hr-vs-drought-0hr-lfc.csv")

###drought time 6 vs predrought time 6#####
###create a deseq2 object with Condition, Time, and Site as design###
coldata$Condition_Time <- relevel(coldata$Condition_Time, ref = "PreDrought_6hr")


dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ Condition_Time + Site)

###deifferential expression analysis###
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

res.T6 <- results(dds, contrast = c("Condition_Time","Drought_6hr","PreDrought_6hr"))
plotMA(res.T6, ylim=c(-40,40))

summary(res.T6)

#make volcano plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res.T6, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-40,40), ylim=c(0,40)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res.T6, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res.T6, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#shrink effect size
resLFC.T6 <- lfcShrink(dds, coef="Condition_Time_Drought_6hr_vs_PreDrought_6hr", type="apeglm")
resLFC.T6
plotMA(resLFC.T6, ylim=c(-2,2))

summary(resLFC.T6)

#make volcano plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(resLFC.T6, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-15,10), ylim=c(0,25)))
with(subset(resLFC.T6, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(resLFC.T6, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#export results
write.csv(as.data.frame(res.T6), 
          file="drought-6hr-vs-predrought-6hr.csv")
write.csv(as.data.frame(resLFC.T6),
          file="drought-6hr-vs-predrought-6hr-lfc.csv")

#cmopare results with WGCNA 
resLFC.T6e <- as.data.frame(resLFC.T6)
resLFC.T6e$datExpr <- rownames(resLFC.T6e)

T0.WGCNA.LFC.T6 <- geneInfo0 %>%
  filter(datExpr != "<NA>") %>%
  merge(.,resLFC.T6e, by = "datExpr", all = TRUE) %>%
  arrange(pvalue)

write.csv(T0.WGCNA.LFC.T6, "T0.WGCNA.T6.LFC.csv")

res.T6e <- as.data.frame(res.T6)
res.T6e$datExpr <- rownames(res.T6e)

T0.WGCNA.T6 <- geneInfo0 %>%
  filter(datExpr != "<NA>") %>%
  merge(.,res.T6e, by = "datExpr", all = TRUE) %>%
  arrange(pvalue)

write.csv(T0.WGCNA.T6, "T0.WGCNA.T6.csv")

###drought time 48 vs predrought time 48#####
###create a deseq2 object with Condition, Time, and Site as design###
coldata$Condition_Time <- relevel(coldata$Condition_Time, ref = "PreDrought_48hr")


dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ Condition_Time + Site)

###deifferential expression analysis###
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients

res.t48 <- results(dds, contrast = c("Condition_Time","Drought_48hr","PreDrought_48hr"))
plotMA(res, ylim=c(-40,40))

summary(res)

#make volcano plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res.t48, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-40,40), ylim=c(0,40)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res.t48, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res.t48, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#shrink effect size
resLFC.t48 <- lfcShrink(dds, coef="Condition_Time_Drought_48hr_vs_PreDrought_48hr", type="apeglm")
resLFC.t48
plotMA(resLFC, ylim=c(-2,2))

summary(resLFC.t48)

#make volcano plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(resLFC.t48, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-15,10), ylim=c(0,25)))
with(subset(resLFC.t48, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(resLFC.t48, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#export results
write.csv(as.data.frame(res.t48), 
          file="drought-48hr-vs-predrought-48hr.csv")
write.csv(as.data.frame(resLFC.t48),
          file="drought-48hr-vs-predrought-48hr-lfc.csv")

#cmopare results with WGCNA 
resLFC.t48 <- as.data.frame(resLFC.t48)
resLFC.t48$datExpr <- rownames(resLFC.t48)

T0.WGCNA.LFC.T48 <- geneInfo0 %>%
  filter(datExpr != "<NA>") %>%
  merge(.,resLFC.t48, by = "datExpr", all = TRUE) %>%
  arrange(pvalue)

write.csv(T0.WGCNA.LFC.T48, "T0.WGCNA.T48.LFC.csv")

res.t48 <- as.data.frame(res.t48)
res.t48$datExpr <- rownames(res.t48)

T0.WGCNA.T48 <- geneInfo0 %>%
  filter(datExpr != "<NA>") %>%
  merge(.,res.t48, by = "datExpr", all = TRUE) %>%
  arrange(pvalue)

write.csv(T0.WGCNA.T48, "T0.WGCNA.T48.csv")
####heatmaps and clustering#############
library('pheatmap')
library("RColorBrewer")
library("ggplot2")

col_list = c("lightblue", "blue", "darkblue", "#FFEA33", "#DFCF46", "#BFB559")
#col_list2 = c("#0165F8", "#FFEA33")
col_list2 = c("#17C137", "#C117A1")


#heatmap
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:500]
df <- as.data.frame(colData(dds)[,c("Condition_Time","Site")])
pheatmap(assay(vstcounts)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)

#distance matrix
sampleDists <- dist(t(assay(vstcounts)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vstcounts$Condition_Time, vstcounts$Site, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotDispEsts(dds)

#PCA plot
pcaData <- plotPCA(vstcounts, intgroup=c("Condition_Time", "Site"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

filename=paste("figs/pca-vst.png", sep = "")
png(filename ,width=5, height=4, unit='in', res = 1000)
ggplot(pcaData, aes(PC1, PC2, color=Condition_Time, shape=Site)) +
  geom_point(size=3) +
  scale_color_manual(values = col_list) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw() +
  coord_fixed()
dev.off()

#######plot genes##########
#create dds object based on condition and time for plotting purposes
coldata$Condition <- factor(coldata$Condition, c("PreDrought", "Drought"))
coldata$Time <- factor(coldata$Time, c("0hr", "6hr", "48hr"))

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ Condition + Time)

d <- plotCounts(dds, gene="KO:K00163", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K00163.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K00163") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K01647", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K01647.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K01647") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K01574", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K01574.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K01574") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K01035", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K01035.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K01035") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K01653", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K01653.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K01653") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K01067", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K01067.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
d1 <- ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K01067") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K00001", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K00001.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K00001") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K14028", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K14028.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K14028") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K00138", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K00138.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K00138") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K00132", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K00132.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K00132") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()



d <- plotCounts(dds, gene="KO:K01034", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K01034.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K01034") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()


d <- plotCounts(dds, gene="KO:K00626", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K00626.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K00626") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K01512", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K01512.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K01512") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K00925", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K00925.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K00925") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()



d <- plotCounts(dds, gene="KO:K00625", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K00625.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K00625") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K10854", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K10854.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K10854") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K00171", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K00171.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K00171") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K00169", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K00169.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K00169") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K00170", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K00170.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K00170") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K00048", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K00048.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K00048") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K18371", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K18371.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K18371") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K18382", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K18382.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K18372") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K00128", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K00128.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K00128") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K01067", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K01067.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K01067") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1, hjust = 1, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()


d <- plotCounts(dds, gene="KO:K00925", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K00925.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K00925") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()



d <- plotCounts(dds, gene="KO:K00132", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K00132.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K00132") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K01895", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K01895.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K01895") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K00138", intgroup="Condition_Time", 
                returnData=TRUE)

filename=paste("figs/K00138.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Condition_Time, y=count, fill = Condition_Time)) + 
  geom_point() + geom_boxplot() + 
  scale_y_log10() +
  scale_fill_manual (values = col_list) +
  labs(title = "K00138") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

#for figures
d <- plotCounts(dds, gene="KO:K00156", intgroup=c("Condition", "Time"),
                returnData=TRUE)
filename=paste("figs/K00156.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Time, y=count, fill = Condition)) + 
  geom_point() + geom_boxplot() + 
  facet_grid(~Condition) +
  scale_y_log10() +
  scale_fill_manual (values = col_list2) +
  labs(title = "K00156:PDH-quinone") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K00158", intgroup=c("Condition", "Time"),
                returnData=TRUE)

filename=paste("figs/K00158.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Time, y=count, fill = Condition)) + 
  geom_point() + geom_boxplot() + 
  facet_grid(~Condition) +
  scale_y_log10() +
  scale_fill_manual (values = col_list2) +
  labs(title = "K00158:pyruvate oxidase") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K01067", intgroup=c("Condition", "Time"),
                returnData=TRUE)

filename=paste("figs/K01067.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Time, y=count, fill = Condition)) + 
  geom_point() + geom_boxplot() + 
  facet_grid(~Condition) +
  scale_y_log10() +
  scale_fill_manual (values = col_list2) +
  labs(title = "K01067:acetyl-CoA hydrolase") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K01026", intgroup=c("Condition", "Time"),
                returnData=TRUE)

filename=paste("figs/K01026.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Time, y=count, fill = Condition)) + 
  geom_point() + geom_boxplot() + 
  facet_grid(~Condition) +
  scale_y_log10() +
  scale_fill_manual (values = col_list2) +
  labs(title = "K01026:propionate CoA-transferase") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K18118", intgroup=c("Condition", "Time"), 
                returnData=TRUE)

filename=paste("figs/K18118.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Time, y=count, fill = Condition)) + 
  geom_point() + geom_boxplot() + 
  facet_grid(~Condition) +
  scale_y_log10() +
  scale_fill_manual (values = col_list2) +
  labs(title = "K18118:acetate CoA-transferase") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K01895", intgroup=c("Condition", "Time"), 
                returnData=TRUE)

filename=paste("figs/K01895.png", sep = "")
png(filename ,width=6, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Time, y=count, fill = Condition)) + 
  geom_point() + geom_boxplot() + 
  facet_grid(~Condition) +
  scale_y_log10() +
  scale_fill_manual (values = col_list2) +
  labs(title = "K01895:acetyl-CoA synthetase") +
  theme_bw() +
  theme(text = element_text(size = 16,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black",face = "bold"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K01574", intgroup=c("Condition", "Time"), 
                returnData=TRUE)

filename=paste("figs/K01574.png", sep = "")
png(filename ,width=7, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Time, y=count, fill = Condition)) + 
  geom_point() + geom_boxplot() + 
  facet_grid(~Condition) +
  scale_y_log10() +
  scale_fill_manual (values = col_list2) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  labs(title = "K01574:acetoacetate decarboxylase", y = "Count", x = "Time post pyruvate addition (hr)") +
  theme_bw() +
  theme(text = element_text(size = 14,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K18382", intgroup=c("Condition", "Time"), 
                returnData=TRUE)

filename=paste("figs/K18382.png", sep = "")
png(filename ,width=7, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Time, y=count, fill = Condition)) + 
  geom_point() + geom_boxplot() + 
  facet_grid(~Condition) +
  scale_y_log10() +
  scale_fill_manual (values = col_list2) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  labs(title = "K18382:secondary alcohol dehydrogenase", y = "Count", x = "Time post pyruvate addition (hr)") +
  theme_bw() +
  theme(text = element_text(size = 14,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K18371", intgroup=c("Condition", "Time"), 
                returnData=TRUE)

filename=paste("figs/K18371.png", sep = "")
png(filename ,width=7, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Time, y=count, fill = Condition)) + 
  geom_point() + geom_boxplot() + 
  facet_grid(~Condition) +
  scale_y_log10() +
  scale_fill_manual (values = col_list2) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  labs(title = "K18371:acetone monooxygenase", y = "Count", x = "Time post pyruvate addition (hr)") +
  theme_bw() +
  theme(text = element_text(size = 14,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K01907", intgroup=c("Condition", "Time"), 
                returnData=TRUE)

filename=paste("figs/K01907.png", sep = "")
png(filename ,width=7, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Time, y=count, fill = Condition)) + 
  geom_point() + geom_boxplot() + 
  facet_grid(~Condition) +
  scale_y_log10() +
  scale_fill_manual (values = col_list2) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  labs(title = "K01907:acetoacetyl-CoA synthetase", y = "Count", x = "Time post pyruvate addition (hr)") +
  theme_bw() +
  theme(text = element_text(size = 14,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K01027", intgroup=c("Condition", "Time"), 
                returnData=TRUE)

filename=paste("figs/K01027.png", sep = "")
png(filename ,width=7, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Time, y=count, fill = Condition)) + 
  geom_point() + geom_boxplot() + 
  facet_grid(~Condition) +
  scale_y_log10() +
  scale_fill_manual (values = col_list2) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  labs(title = "K01027:3-oxoacid ClA-transferase", y = "Count", x = "Time post pyruvate addition (hr)") +
  theme_bw() +
  theme(text = element_text(size = 14,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K01034", intgroup=c("Condition", "Time"), 
                returnData=TRUE)

filename=paste("figs/K01034.png", sep = "")
png(filename ,width=7, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Time, y=count, fill = Condition)) + 
  geom_point() + geom_boxplot() + 
  facet_grid(~Condition) +
  scale_y_log10() +
  scale_fill_manual (values = col_list2) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  labs(title = "K01034:acetoacetate CoA transferase (α)", y = "Count", x = "Time post pyruvate addition (hr)") +
  theme_bw() +
  theme(text = element_text(size = 14,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K00656", intgroup=c("Condition", "Time"), 
                returnData=TRUE)

filename=paste("figs/K00656.png", sep = "")
png(filename ,width=7, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Time, y=count, fill = Condition)) + 
  geom_point() + geom_boxplot() + 
  facet_grid(~Condition) +
  scale_y_log10() +
  scale_x_discrete(labels = c("0", "6", "48")) +
  scale_fill_manual (values = col_list2) +
  labs(title = "K00656:formate C-acetyltransferase", y = "Count", x = "Time post pyruvate addition (hr)") +
  theme_bw() +
  theme(text = element_text(size = 14,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K00169", intgroup=c("Condition", "Time"), 
                returnData=TRUE)

filename=paste("figs/K00169.png", sep = "")
png(filename ,width=7, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Time, y=count, fill = Condition)) + 
  geom_point() + geom_boxplot() + 
  facet_grid(~Condition) +
  scale_y_log10() +
  scale_fill_manual (values = col_list2) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  labs(title = "K00169:pyruvate ferredoxin oxidoreductase (α)", y = "Count", x = "Time post pyruvate addition (hr)") +
  theme_bw() +
  theme(text = element_text(size = 14,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"))
dev.off()

d <- plotCounts(dds, gene="KO:K00174", intgroup=c("Condition", "Time"), 
                returnData=TRUE)

filename=paste("figs/K00174.png", sep = "")
png(filename ,width=7, height=5, unit='in', res = 1000)
ggplot(d, aes(x=Time, y=count, fill = Condition)) + 
  geom_point() + geom_boxplot() + 
  facet_grid(~Condition) +
  scale_y_log10() +
  scale_fill_manual (values = col_list2) +
  scale_x_discrete(labels = c("0", "6", "48")) +
  labs(title = "K00174:2-oxoacid ferredoxin oxidoreductase (α)", y = "Count", x = "Time post pyruvate addition (hr)") +
  theme_bw() +
  theme(text = element_text(size = 14,
                            family = "Arial",
                            color = "black"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black"),
        axis.text.y = element_text(color = "black"))
dev.off()


#####Extract acetate, acetone, and diacetyl genes of interest from combined WGCNA files created above

interest <- read.csv("KO_interest.csv", header = TRUE)

##all times

T0.WGCNA.dp.f <- T0.WGCNA.dp %>%
  merge(., interest, by = "datExpr") %>%
  select(-contains("MM")) %>%
  filter(Priority == "yes") %>%
  mutate(Significance =
           case_when(padj < 0.05 ~ "Sig", 
                     padj > 0.05 ~ "NotSig")) %>%
  mutate(datExpr = str_replace(datExpr, 'KO:', '')) %>%
  arrange(.,log2FoldChange)


gene.order <- T0.WGCNA.dp.f$gene_name
T0.WGCNA.dp.f$gene_name <- factor(T0.WGCNA.dp.f$gene_name, gene.order)



filename=paste("figs/log2fc-dp-modules.png", sep = "")
png(filename ,width=4, height=5, unit='in', res = 1000)
ggplot(T0.WGCNA.dp.f,aes(x=datExpr, y=log2FoldChange, fill = ModuleColor, color = Activity)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  facet_wrap(~VOC, scales = "free", ncol = 1) +
  scale_fill_manual (breaks = c("red", "blue", "black" , "brown", "turquoise", "yellow","pink", "grey"),
                     values = c("red", "blue", "black" , "brown", "turquoise", "yellow", "pink", "grey"))
dev.off()

filename=paste("figs/log2fc-dp-v2.png", sep = "")
png(filename ,width=3, height=5, unit='in', res = 1000)
ggplot(T0.WGCNA.dp.f,aes(x=gene_name, y=log2FoldChange, fill = Sig_Cat)) +
  geom_col() +
  ylim(-1.5, 6.5) +
  coord_flip() +
  theme_bw() +
  labs(y="Log2FC", x="") +
  facet_wrap(~VOC, scales = "free", ncol = 1) +
  scale_fill_manual (values =  c("#2ECC71" ,"#F39C12"),
                     breaks = c("Pre-Drought","produce")) +
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.title = element_blank() ,
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
dev.off()

#time 0
T0.WGCNA.f <- T0.WGCNA %>%
  merge(., interest, by = "datExpr") %>%
  select(-contains("MM")) %>%
  filter(Priority == "yes") %>%
  mutate(Significance =
           case_when(padj < 0.05 ~ "Sig", 
                     padj > 0.05 ~ "NotSig")) %>%
  mutate(datExpr = str_replace(datExpr, 'KO:', '')) %>%
  arrange(.,log2FoldChange)


gene.order <- T0.WGCNA.f$gene_name
T0.WGCNA.f$gene_name <- factor(T0.WGCNA.f$gene_name, gene.order)

filename=paste("figs/log2fc-modules.png", sep = "")
png(filename ,width=4, height=5, unit='in', res = 1000)
ggplot(T0.WGCNA.f,aes(x=datExpr, y=log2FoldChange, fill = ModuleColor, color = Activity)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  facet_wrap(~VOC,  scales = "free", ncol = 1) +
  scale_fill_manual (breaks = c("red", "blue", "black" , "brown", "turquoise", "yellow", "pink", "grey"),
                     values = c("red", "blue", "black" , "brown", "turquoise", "yellow", "pink", "grey"))
dev.off()

filename=paste("figs/log2fc-t0-v2.png", sep = "")
png(filename ,width=3, height=5, unit='in', res = 1000)
ggplot(T0.WGCNA.f,aes(x=gene_name, y=log2FoldChange, fill = Activity)) +
  geom_col() +
  ylim(-1.5, 6.5) +
  coord_flip() +
  theme_bw() +
  labs(y="Log-2FC") +
  facet_wrap(~VOC, scales = "free", ncol = 1) +
  scale_fill_manual (values = c("lightgray","black"),
                     breaks = c("consume", "produce")) +
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.title = element_blank() ,
        legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.title.x = element_text(face= "bold"), 
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
dev.off()

##time 6

T0.WGCNA.T6.f <- T0.WGCNA.T6 %>%
  merge(., interest, by = "datExpr") %>%
  select(-contains("MM")) %>%
  filter(Priority == "yes") %>%
  mutate(Significance =
           case_when(padj < 0.05 ~ "Sig", 
                     padj > 0.05 ~ "NotSig")) %>%
  mutate(datExpr = str_replace(datExpr, 'KO:', '')) %>%
  arrange(.,log2FoldChange)

#gene.order <- T0.WGCNA.T6.f$datExpr
T0.WGCNA.T6.f$gene_name <- factor(T0.WGCNA.T6.f$gene_name, gene.order)

filename=paste("figs/log2fc-t6-modules.png", sep = "")
png(filename ,width=4, height=5, unit='in', res = 1000)
ggplot(T0.WGCNA.T6.f,aes(x=datExpr, y=log2FoldChange, fill = ModuleColor, color = activity)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  facet_wrap(~VOC, scales = "free", ncol = 1) +
  scale_fill_manual (breaks = c("red", "blue", "black" , "brown", "turquoise", "yellow", "grey"),
                     values = c("red", "blue", "black" , "brown", "turquoise", "yellow", "grey"))
dev.off()

filename=paste("figs/log2fc-t6-v2.png", sep = "")
png(filename ,width=3, height=5, unit='in', res = 1000)
ggplot(T0.WGCNA.T6.f,aes(x=gene_name, y=log2FoldChange, fill = Activity)) +
  geom_col() +
  ylim(-1.5, 6.5) +
  coord_flip() +
  theme_bw() +
  labs(y="Log-2FC") +
  facet_wrap(~VOC, scales = "free", ncol = 1) +
  scale_fill_manual (values = c("lightgray","black"),
                     breaks = c("consume", "produce")) +
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.title = element_blank() ,
        legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.title.x = element_text(face= "bold"), 
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))

##time 48
T0.WGCNA.T48.f <- T0.WGCNA.T48 %>%
  merge(., interest, by = "datExpr") %>%
  select(-contains("MM")) %>%
  filter(Priority == "yes") %>%
  mutate(Significance =
           case_when(padj < 0.05 ~ "Sig", 
                     padj > 0.05 ~ "NotSig")) %>%
  mutate(datExpr = str_replace(datExpr, 'KO:', '')) %>%
  arrange(.,log2FoldChange)


#gene.order <- T0.WGCNA.T48.f$datExpr
T0.WGCNA.T48.f$gene_name <- factor(T0.WGCNA.T48.f$gene_name, gene.order)

filename=paste("figs/log2fc-t48-modules.png", sep = "")
png(filename ,width=4, height=5, unit='in', res = 1000)
ggplot(T0.WGCNA.T48.f,aes(x=datExpr, y=log2FoldChange, fill = ModuleColor, color = activity)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  facet_wrap(~VOC, scales = "free", ncol = 1) +
  scale_fill_manual (breaks = c("red", "blue", "black" , "brown", "turquoise", "yellow", "grey"),
                     values = c("red", "blue", "black" , "brown", "turquoise", "yellow", "grey"))
dev.off()

filename=paste("figs/log2fc-t48-v2.png", sep = "")
png(filename ,width=3, height=5, unit='in', res = 1000)
geom_col() +
  ylim(-1.5, 6.5) +
  coord_flip() +
  theme_bw() +
  labs(y="Log-2FC") +
  facet_wrap(~VOC, scales = "free", ncol = 1) +
  scale_fill_manual (values = c("lightgray","black"),
                     breaks = c("consume", "produce")) +
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.title = element_blank() ,
        legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.title.x = element_text(face= "bold"), 
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
dev.off()

#####Extract acetate, acetone, and diacetyl genes of interest from combined WGCNA files created above - repeated with out using lfc

interest <- read.csv("KO_interest.csv", header = TRUE)

##all times

T0.WGCNA.dp.f <- T0.WGCNA.dp %>%
  merge(., interest, by = "datExpr") %>%
  select(-contains("MM")) %>%
  filter(Priority == "yes") %>%
  mutate(Significance =
           case_when(padj < 0.05 ~ "Sig", 
                     padj > 0.05 ~ "NotSig")) %>%
  mutate(datExpr = str_replace(datExpr, 'KO:', '')) %>%
  arrange(.,log2FoldChange)

gene.order <- T0.WGCNA.dp.f$gene_name
T0.WGCNA.dp.f$gene_name <- factor(T0.WGCNA.dp.f$gene_name, gene.order)



filename=paste("figs/log2fc-dp-modules.png", sep = "")
png(filename ,width=4, height=5, unit='in', res = 1000)
ggplot(T0.WGCNA.dp.f,aes(x=datExpr, y=log2FoldChange, fill = ModuleColor, color = Activity)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  facet_wrap(~VOC, scales = "free", ncol = 1) +
  scale_fill_manual (breaks = c("red", "blue", "black" , "brown", "turquoise", "yellow","pink", "grey"),
                     values = c("red", "blue", "black" , "brown", "turquoise", "yellow", "pink", "grey"))
dev.off()

filename=paste("figs/log2fc-dp-v2.png", sep = "")
png(filename ,width=3, height=5, unit='in', res = 1000)
ggplot(T0.WGCNA.dp.f,aes(x=gene_name, y=log2FoldChange, fill = Activity)) +
  geom_col() +
  ylim(-6, 6) +
  coord_flip() +
  theme_bw() +
  labs(y="Log-2FC") +
  facet_wrap(~VOC, scales = "free", ncol = 1) +
  scale_fill_manual (values = c("lightgray","black"),
                     breaks = c("consume", "produce")) +
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.title = element_blank() ,
        legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.title.x = element_text(face= "bold"), 
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
dev.off()

#time 0
T0.WGCNA.f <- T0.WGCNA %>%
  merge(., interest, by = "datExpr") %>%
  select(-contains("MM")) %>%
  filter(Priority == "yes") %>%
  mutate(Significance =
           case_when(padj < 0.05 ~ "Sig", 
                     padj > 0.05 ~ "NotSig")) %>%
  mutate(datExpr = str_replace(datExpr, 'KO:', '')) %>%
  arrange(.,log2FoldChange)


gene.order <- T0.WGCNA.f$gene_name
T0.WGCNA.f$gene_name <- factor(T0.WGCNA.f$gene_name, gene.order)

filename=paste("figs/log2fc-modules.png", sep = "")
png(filename ,width=4, height=5, unit='in', res = 1000)
ggplot(T0.WGCNA.f,aes(x=datExpr, y=log2FoldChange, fill = ModuleColor, color = Activity)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  facet_wrap(~VOC,  scales = "free", ncol = 1) +
  scale_fill_manual (breaks = c("red", "blue", "black" , "brown", "turquoise", "yellow", "pink", "grey"),
                     values = c("red", "blue", "black" , "brown", "turquoise", "yellow", "pink", "grey"))
dev.off()

filename=paste("figs/log2fc-t0-v2.png", sep = "")
png(filename ,width=3, height=5, unit='in', res = 1000)
ggplot(T0.WGCNA.f,aes(x=gene_name, y=log2FoldChange, fill = Activity)) +
  geom_col() +
  ylim(-20, 20) +
  coord_flip() +
  theme_bw() +
  labs(y="Log2-FC") +
  facet_wrap(~VOC, scales = "free", ncol = 1) +
  scale_fill_manual (values = c("lightgray","black"),
                     breaks = c("consume", "produce")) +
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.title = element_blank() ,
        legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.title.x = element_text(face= "bold"), 
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
dev.off()

##time 6

T0.WGCNA.T6.f <- T0.WGCNA.T6 %>%
  merge(., interest, by = "datExpr") %>%
  select(-contains("MM")) %>%
  filter(Priority == "yes") %>%
  mutate(Significance =
           case_when(padj < 0.05 ~ "Sig", 
                     padj > 0.05 ~ "NotSig")) %>%
  mutate(datExpr = str_replace(datExpr, 'KO:', '')) %>%
  arrange(.,log2FoldChange)

#gene.order <- T0.WGCNA.T6.f$datExpr
T0.WGCNA.T6.f$gene_name <- factor(T0.WGCNA.T6.f$gene_name, gene.order)

filename=paste("figs/log2fc-t6-modules.png", sep = "")
png(filename ,width=4, height=5, unit='in', res = 1000)
ggplot(T0.WGCNA.T6.f,aes(x=datExpr, y=log2FoldChange, fill = ModuleColor, color = activity)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  facet_wrap(~VOC, scales = "free", ncol = 1) +
  scale_fill_manual (breaks = c("red", "blue", "black" , "brown", "turquoise", "yellow", "grey"),
                     values = c("red", "blue", "black" , "brown", "turquoise", "yellow", "grey"))
dev.off()

filename=paste("figs/log2fc-t6-v2.png", sep = "")
png(filename ,width=3, height=5, unit='in', res = 1000)
ggplot(T0.WGCNA.T6.f,aes(x=gene_name, y=log2FoldChange, fill = Activity)) +
  geom_col() +
  ylim(-6, 6) +
  coord_flip() +
  theme_bw() +
  labs(y="Log2-FC") +
  facet_wrap(~VOC, scales = "free", ncol = 1) +
  scale_fill_manual (values = c("lightgray","black"),
                     breaks = c("consume", "produce")) +
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.title = element_blank() ,
        legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.title.x = element_text(face= "bold"), 
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
dev.off()

##time 48
T0.WGCNA.T48.f <- T0.WGCNA.T48 %>%
  merge(., interest, by = "datExpr") %>%
  select(-contains("MM")) %>%
  filter(Priority == "yes") %>%
  mutate(Significance =
           case_when(padj < 0.05 ~ "Sig", 
                     padj > 0.05 ~ "NotSig")) %>%
  mutate(datExpr = str_replace(datExpr, 'KO:', '')) %>%
  arrange(.,log2FoldChange)


#gene.order <- T0.WGCNA.T48.f$datExpr
T0.WGCNA.T48.f$gene_name <- factor(T0.WGCNA.T48.f$gene_name, gene.order)

filename=paste("figs/log2fc-t48-modules.png", sep = "")
png(filename ,width=4, height=5, unit='in', res = 1000)
ggplot(T0.WGCNA.T48.f,aes(x=datExpr, y=log2FoldChange, fill = ModuleColor, color = activity)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  facet_wrap(~VOC, scales = "free", ncol = 1) +
  scale_fill_manual (breaks = c("red", "blue", "black" , "brown", "turquoise", "yellow", "grey"),
                     values = c("red", "blue", "black" , "brown", "turquoise", "yellow", "grey"))
dev.off()

filename=paste("figs/log2fc-t48-v2.png", sep = "")
png(filename ,width=3, height=5, unit='in', res = 1000)
ggplot(T0.WGCNA.T48.f,aes(x=gene_name, y=log2FoldChange, fill = Activity)) +
  geom_col() +
  ylim(-6,6) +
  coord_flip() +
  theme_bw() +
  labs(y="Log2-FC") +
  facet_wrap(~VOC, scales = "free", ncol = 1) +
  scale_fill_manual (values = c("lightgray","black"),
                     breaks = c("consume", "produce")) +
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.title = element_blank() ,
        legend.position = "bottom",
        axis.title.y = element_blank(),
        axis.title.x = element_text(face= "bold"), 
        axis.text.x = element_text(vjust = 1, hjust = 0.5, color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10))
dev.off()


