##Linnea Honeker
##PCA on MNR data from root pyruvate experiment
##9/25/20

library(ggfortify)
library(tidyverse)
library(factoextra)
library(ggnewscale)
library(ggrepel)
library(vegan)
library(dplyr)
library(pairwiseAdonis)
library(remotes)
library(EcolUtils)

# load source code
source('./RScripts/FigS4-S7-PCA/extra_functions.R')


###parameters for plots
h = 6.5
w = 9
res = 300
size = 13
my_colors_d = c("#2ECC71" ,"#F39C12")
my_colors_st = c("#7CDFA6","#2ECC71", "#1D8047", "#F7C16B", "#F39C12", "#98620B")

list_of_shapes <- c(15,16,17)
########################################## PCA on log-transformed metaT############

df.0 <- read.csv(file="./Data/PCA/Drought_vs_predrought_metaT_gene_copies.csv", header = TRUE, blank.lines.skip = TRUE)
metadata <- read.csv(file="./Data/PCA/metaT-metadata.csv", header = TRUE)
rownames(metadata) <- metadata[,1]
names(metadata) [1] <- "SampleID"

#reformat - transpose and take log
rownames(df.0) <- df.0$Feature
df.0 <- df.0[,-1]
df <- t(df.0)
df <- as.data.frame(df)
df[is.na(df)] <- 0 #convert any na's to zeros

#remove row that doesn't match between tables
#df <- df[-25,]

#delete any columns (KOs) that sum to zero
df <- df %>%
  dplyr::select(which(!colSums(., na.rm=TRUE) %in% 0))
colSums(df)

df <- df + 1 #add psuedocount
df <- log(df) #log transform


########################################## PCA 
# Calculate PCA with prcomp()
pca <- prcomp(df, scale = FALSE, center = TRUE)
# Extract eigenvalues and variances
eigen <- get_eigenvalue(pca) # this function is from factoextra
dimensions <- c(1:dim(eigen)[1]) # this is for the plot


scree <- make_screeplot(eigen, dimensions) + ggtitle('Scree plot, PCA Class')
filename <- paste0("./Figures/FS4-S7-PCA/genomic/PCA_screeplot-log.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,scree)
cumvar <-  make_cumvar(eigen, dimensions) + ggtitle('Cumulative variance plot, PCA Class')
filename <- paste0("./Figures/FS4-S7-PCA/genomic/PCA_cumulative_var-log.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,cumvar)

# extract coordinates for PC1 and PC2
pca_results <- get_pca_ind(pca) #pca[["x"]]
pca_coordinates <- as.data.frame(pca_results$coord[,c(1,2)])
colnames(pca_coordinates) <- c('PC1','PC2')

# merge metadata with pca_coordinates:
pca_coordinates$SampleID <- rownames(pca_coordinates)
pca_coordinates <- merge(pca_coordinates, metadata, by = "SampleID")
pca_coordinates$Site <- as.factor(pca_coordinates$Site)
pca_coordinates$Time <- as.factor(pca_coordinates$Time)
pca_coordinates$Condition <- factor(pca_coordinates$Condition, levels = c("PreDrought", "Drought"))
pca_coordinates$Time <- factor(pca_coordinates$Time, levels = c("0hr", "6hr", "48hr"))
pca_coordinates$Condition_Time <- factor(pca_coordinates$Condition_Time, 
                                         levels = c("PreDrought_0hr", "PreDrought_6hr", "PreDrought_48hr",
                                                    "Drought_0hr", "Drought_6hr", "Drought_48hr"))

write.csv(pca_coordinates,file=paste0("./Output/PCA_genomic/PCA_individual_coordinates-log.csv"),row.names=TRUE)

# prepare label for graph
pc1 <- paste0('PC1 (',round(eigen$variance.percent[1],digits=1),'%)')
pc2 <- paste0('PC2 (',round(eigen$variance.percent[2],digits=1),'%)')

# arrows
arrows <- get_arrows(pca, pca_coordinates)
write.csv(arrows,file=paste0("./Output/PCA_genomic/PCA_vector_coordinates-log.csv"), row.names=TRUE)


# pca biplot
#pca_biplot <- ggplot(mapping = aes(x, y)) +
#  #new_scale_color() +
#  geom_segment(data=arrows, aes(x=0, y=0, xend=xend, yend=yend),
#               arrow=arrow(length = unit(0.1,"cm")), size=0.7, color = "grey") +
#  geom_text_repel(data=arrows, aes(x=xend, y=yend), color = "black",
##                  label = arrows$name,
#                  size=size-6, show.legend = FALSE)

#pca_biplot 



#pca condition
pca_plot <-  #pca_biplot +
  ggplot() +
  geom_point(data=pca_coordinates, aes(x=PC1, y=PC2, col=Condition), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= pc1, y=pc2) +
  scale_color_manual(values = my_colors_d) +
  #scale_shape_manual(values= list_of_shapes) +
  theme( legend.text = element_text(size=size+3, face="bold"),
         legend.title = element_blank(),
         legend.key.size = unit(0.6, "cm"),
         legend.key.width = unit(0.6,"cm"),
         #legend.position = "bottom",
         panel.grid = element_blank(),
         axis.title.x = element_text(size=size+3,face="bold"),
         axis.title.y = element_text(size=size+3,face="bold"),
         plot.title = element_text(size=size+3,face="bold")) 
pca_plot
filename <- paste0("./Figures/PCA_genomic/metaT-PCA-time.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

#pca time_site
pca_plot <-  #pca_biplot +
  ggplot() +
  geom_point(data=pca_coordinates, aes(x=PC1, y=PC2, col=Condition_Time, shape = Site), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= pc1, y=pc2) +
  scale_color_manual(values = my_colors_st) +
  scale_shape_manual(values= list_of_shapes) +
  theme( legend.text = element_text(size=size, face="bold"),
         legend.title = element_blank(),
         #legend.key.size = unit(0.6, "cm"),
         #legend.key.width = unit(0.6,"cm"),
         legend.position = "bottom",
         panel.grid = element_blank(),
         axis.title.x = element_text(size=size,face="bold"),
         axis.title.y = element_text(size=size,face="bold"),
         plot.title = element_text(size=size,face="bold")) 
pca_plot
filename <- paste0("./Figures/FigS4-S7-PCA/genomic/FigS4-metaT-PCA-condition-time.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

###determine if pre-drought vs drought are significantly different using permanova
df.dist <- vegdist(df, method="bray")

adonis(df.dist ~ Condition * Time * Site, metadata, permutations = 999, method = "bray",
       strata = NULL)#, contr.unordered = "contr.sum",
       #contr.ordered = "contr.poly")

#Results:
#Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
#Condition            1   0.16555 0.165555 2.41385 0.06508  0.038 *
#  Time                 2   0.29554 0.147771 2.15455 0.11617  0.024 *
#  Site                 2   0.19870 0.099349 1.44855 0.07811  0.122  
#Condition:Time       2   0.14027 0.070135 1.02259 0.05514  0.397  
#Condition:Site       2   0.15045 0.075227 1.09684 0.05914  0.303  
#Time:Site            4   0.26527 0.066318 0.96694 0.10428  0.462  
#Condition:Time:Site  4   0.16220 0.040550 0.59124 0.06376  0.951  
#Residuals           17   1.16595 0.068585         0.45832         
#Total               34   2.54394                  1.00000         
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(df.dist ~ Condition_Time, metadata)

#Results:
#Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
#Condition_Time  5   0.60705 0.12141  1.8178 0.23862  0.017 *
#  Residuals      29   1.93690 0.06679         0.76138         
#Total          34   2.54394                 1.00000         


########################################## PCA on log-transformed metaG############

df <- read.csv("./Data/PCA/metaG-gene-copies-no-ctrl.csv", header = TRUE)
metadata <- read.csv("./Data/PCA/metadata-metaG-no-ctrl.csv", header = TRUE)

#reformat - transpose and take log
rownames(df) <- df$Feature
df <- df[,-1]
df <- t(df)
df <- as.data.frame(df)
df[is.na(df)] <- 0 #convert any na's to zeros

#remove row that doesn't match between tables
#df <- df[-25,]

#delete any columns (KOs) that sum to zero
df <- df %>%
  select(which(!colSums(., na.rm=TRUE) %in% 0))
colSums(df)

df <- df + 1 #add psuedocount
df <- log(df) #log transform


########################################## PCA 
# Calculate PCA with prcomp()
pca <- prcomp(df, scale = FALSE, center = TRUE)
# Extract eigenvalues and variances
eigen <- get_eigenvalue(pca) # this function is from factoextra
dimensions <- c(1:dim(eigen)[1]) # this is for the plot

#


scree <- make_screeplot(eigen, dimensions) + ggtitle('Scree plot, PCA Class')
filename <- paste0("./Figures/FigS4-S7-PCA/genomic/PCA_screeplot-log-metaG.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,scree)
cumvar <-  make_cumvar(eigen, dimensions) + ggtitle('Cumulative variance plot, PCA Class')
filename <- paste0("./Figures/FigS4-S7-PCA/genomic/PCA_cumulative_var-log-metaG.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,cumvar)

# extract coordinates for PC1 and PC2
pca_results <- get_pca_ind(pca) #pca[["x"]]
pca_coordinates <- as.data.frame(pca_results$coord[,c(1,2)])
colnames(pca_coordinates) <- c('PC1','PC2')

# merge metadata with pca_coordinates:
pca_coordinates$Sample_ID <- rownames(pca_coordinates)
pca_coordinates <- merge(pca_coordinates, metadata, by = "Sample_ID")
pca_coordinates$Site <- as.factor(pca_coordinates$Site)
pca_coordinates$Time <- as.factor(pca_coordinates$Time)
pca_coordinates$Condition <- factor(pca_coordinates$Condition, levels = c("PreDrought", "Drought"))
pca_coordinates$Condition_time <- factor(pca_coordinates$Condition_time, 
                                         levels = c("PreDrought_0hr", "PreDrought_6hr", "PreDrought_48hr",
                                                    "Drought_0hr", "Drought_6hr", "Drought_48hr"))

write.csv(pca_coordinates,file=paste0("./Output/PCA_genomic/PCA_individual_coordinates-log-metaG.csv"),row.names=TRUE)

# prepare label for graph
pc1 <- paste0('PC1 (',round(eigen$variance.percent[1],digits=1),'%)')
pc2 <- paste0('PC2 (',round(eigen$variance.percent[2],digits=1),'%)')

# arrows
arrows <- get_arrows(pca, pca_coordinates)
write.csv(arrows,file=paste0("./Output/PCA_genomic/PCA_vector_coordinates-log-metaG.csv"), row.names=TRUE)


# pca biplot
#pca_biplot <- ggplot(mapping = aes(x, y)) +
#  #new_scale_color() +
#  geom_segment(data=arrows, aes(x=0, y=0, xend=xend, yend=yend),
#               arrow=arrow(length = unit(0.1,"cm")), size=0.7, color = "grey") +
#  geom_text_repel(data=arrows, aes(x=xend, y=yend), color = "black",
##                  label = arrows$name,
#                  size=size-6, show.legend = FALSE)

#pca_biplot 



#pca condition
pca_plot <-  #pca_biplot +
  ggplot() +
  geom_point(data=pca_coordinates, aes(x=PC1, y=PC2, col=Condition), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= pc1, y=pc2) +
  scale_color_manual(values = my_colors_d) +
  #scale_shape_manual(values= list_of_shapes) +
  theme( legend.text = element_text(size=size+3, face="bold"),
         legend.title = element_blank(),
         legend.key.size = unit(0.6, "cm"),
         legend.key.width = unit(0.6,"cm"),
         legend.position = "bottom",
         panel.grid = element_blank(),
         axis.title.x = element_text(size=size+3,face="bold"),
         axis.title.y = element_text(size=size+3,face="bold"),
         plot.title = element_text(size=size+3,face="bold")) 
pca_plot
filename <- paste0("./Figures/FigS4-S7-PCA/genomic/metaG-PCA-condition.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

#pca time_site
pca_plot <-  #pca_biplot +
  ggplot() +
  geom_point(data=pca_coordinates, aes(x=PC1, y=PC2, col=Condition_time, shape = Site), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= pc1, y=pc2) +
  scale_color_manual(values = my_colors_st) +
  scale_shape_manual(values= list_of_shapes) +
  theme( legend.text = element_text(size=size, face="bold"),
         legend.title = element_blank(),
         #legend.key.size = unit(0.6, "cm"),
         #legend.key.width = unit(0.6,"cm"),
         legend.position = "bottom",
         panel.grid = element_blank(),
         axis.title.x = element_text(size=size,face="bold"),
         axis.title.y = element_text(size=size,face="bold"),
         plot.title = element_text(size=size,face="bold")) 
pca_plot
filename <- paste0("./Figures/FigS4-S7-PCA/genomic/FigS4-metaG-PCA-site.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

###determine if pre-drought vs drought are significantly different using permanova
df.dist <- vegdist(df, method="bray")

adonis(df.dist ~ Condition * Time * Site, metadata, permutations = 999, method = "bray",
       strata = NULL, contr.unordered = "contr.sum",
       contr.ordered = "contr.poly")

#Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)    
#Df SumsOfSqs   MeanSqs F.Model      R2 Pr(>F)   
#Condition            1  0.001740 0.0017404  0.6524 0.01562  0.676   
#Time                 2  0.004765 0.0023824  0.8930 0.04277  0.499   
#Site                 2  0.018694 0.0093469  3.5035 0.16779  0.003 **
#  Condition:Time       2  0.007137 0.0035684  1.3375 0.06406  0.205   
#Condition:Site       2  0.008427 0.0042136  1.5794 0.07564  0.112   
#Time:Site            4  0.007473 0.0018682  0.7003 0.06707  0.834   
#Condition:Time:Site  4  0.017822 0.0044554  1.6700 0.15996  0.062 . 
#Residuals           17  0.045354 0.0026679         0.40709          
#Total               34  0.111411                   1.00000             
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
