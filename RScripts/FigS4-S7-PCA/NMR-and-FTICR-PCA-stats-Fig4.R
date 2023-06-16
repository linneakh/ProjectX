##Linnea Honeker
##PCA on MNR data from root pyruvate experiment
##9/25/20

library(ggfortify)
library(tidyverse)
library(factoextra)
library(ggnewscale)
library(ggrepel)
library(viridis)
library(vegan)


source('./RScripts/Fig4-S3-PCA/extra_functions.R')

colors = c("#2ECC71" ,"#F39C12")

list_of_shapes <- c(15,16,17,18,3,7,8,9,10,4,11,12,13,14)

# parameters for plots
h = 5
w = 5
res = 300
size = 10


########################################## PCA on all samples############

#NMR
df.0 = read.csv("./Data/nmr_raw/soil_pyruvate_nmr.csv", header = TRUE)
rownames(df.0) <- df.0$X
df.0$X <- NULL

df <- t(df.0)
df <- as.data.frame(df)

#FTICR masses
df.f.0 = read.csv("./Data/FTICR_norm/metabodirect_matrix_features_norm.csv")
rownames(df.f.0) <- df.f.0$Mass
df.f.0$Mass <- NULL

df.f.0 <- df.f.0 %>%
  mutate(Ctrl_Drought_P12_Site1_0 = Soil_Drought_P12_Site1_0_CTRL) %>%
  select(-Soil_Drought_P12_Site1_0_CTRL) %>%
  mutate(Ctrl_Drought_P12_Site1_6 = Soil_Drought_P12_Site1_6_CTRL) %>%
  select(-Soil_Drought_P12_Site1_6_CTRL) %>%
  mutate(Ctrl_Drought_P12_Site1_48 = Soil_Drought_P12_Site1_48_CTRL) %>%
  select(-Soil_Drought_P12_Site1_48_CTRL)
  
df.f <- t(df.f.0)
df.f <- as.data.frame(df.f)

#FTICR classes
df.c.0 = read.csv("./Data/FTICR_norm/metabodirect_class_composition.csv")
rownames(df.c.0) <- df.c.0$SampleID
df.c.0$SampleID <- NULL
df.c.0 <- t(df.c.0)
df.c.0 <- as.data.frame(df.c.0)

df.c.0 <- df.c.0 %>%
  mutate(Ctrl_Drought_P12_Site1_0 = Soil_Drought_P12_Site1_0_CTRL) %>%
  dplyr::select(-Soil_Drought_P12_Site1_0_CTRL) %>%
  mutate(Ctrl_Drought_P12_Site1_6 = Soil_Drought_P12_Site1_6_CTRL) %>%
  dplyr::select(-Soil_Drought_P12_Site1_6_CTRL) %>%
  mutate(Ctrl_Drought_P12_Site1_48 = Soil_Drought_P12_Site1_48_CTRL) %>%
  dplyr::select(-Soil_Drought_P12_Site1_48_CTRL) %>%
  replace(is.na(.),0)

df.c <- t(df.c.0)
df.c <- as.data.frame(df.c)

########################################## PCA on NMR#####################
# Calculate PCA with prcomp()
pca <- prcomp(df, scale = TRUE, center = TRUE)
# Extract eigenvalues and variances
eigen <- get_eigenvalue(pca) # this function is from factoextra
dimensions <- c(1:dim(eigen)[1]) # this is for the plot

scree <- make_screeplot(eigen, dimensions) + ggtitle('Scree plot, PCA Class')
filename <- paste0("./Figures/PCA_meta/PCA_screeplot.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,scree)
cumvar <-  make_cumvar(eigen, dimensions) + ggtitle('Cumulative variance plot, PCA Class')
filename <- paste0("./Figures/PCA_meta/PCA_cumulative_var.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,cumvar)

# extract coordinates for PC1 and PC2
pca_results <- get_pca_ind(pca) #pca[["x"]]
pca_coordinates <- as.data.frame(pca_results$coord[,c(1,2)])
colnames(pca_coordinates) <- c('PC1','PC2')

# merge metadata with pca_coordinates:
pca_coordinates$SampleID <- rownames(pca_coordinates)
pca_coordinates <- pca_coordinates %>% 
  separate(., SampleID, into = c('unique','Condition',	'Time',	'Site',	'Label'), sep = "_", remove = T)

pca_coordinates$Condition <- factor(pca_coordinates$Condition, c("pre", "drought"))

write.csv(pca_coordinates,file=paste0("./Output/PCA_meta/PCA_individual_coordinates.csv"),row.names=TRUE)

# prepare label for graph
pc1 <- paste0('PC1 (',round(eigen$variance.percent[1],digits=1),'%)')
pc2 <- paste0('PC2 (',round(eigen$variance.percent[2],digits=1),'%)')

# arrows
arrows <- get_arrows(pca, pca_coordinates)
write.csv(arrows,file=paste0("./Output/PCA_meta/PCA_vector_coordinates.csv"), row.names=TRUE)

### pca with no shape 
pca_plot <-  ggplot(mapping = aes(x, y)) +
  geom_point(data=pca_coordinates, aes(x=PC1, y=PC2, col=Condition), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= pc1, y=pc2) +
  scale_color_manual(values = colors,
                     breaks = c("pre", "drought"),
                     labels = c("Pre Drought", "Drought")) +
  scale_shape_manual(values= list_of_shapes) +
  
  theme( legend.text = element_text(size=size+3, face="bold"),
         legend.title = element_blank(),
         legend.key.size = unit(0.6, "cm"),
         legend.key.width = unit(0.6,"cm"),
         legend.position = "bottom",
         axis.title.x = element_text(size=size+3,face="bold"),
         axis.title.y = element_text(size=size+3,face="bold"),
         plot.title = element_text(size=size+3,face="bold")) 
pca_plot
filename <- paste0("./Figures/FigS4-S7-PCA/metabol/NMR_PCA_plot.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

# pca biplot
pca_biplot <- pca_plot +
  new_scale_color() +
  geom_segment(data=arrows, aes(x=0, y=0, xend=xend, yend=yend),
               arrow=arrow(length = unit(0.1,"cm")), size=0.7, color = "grey") +
  geom_text_repel(data=arrows, aes(x=xend, y=yend), color = "black",
                  label=arrows$name, size=size/3, show.legend = FALSE)

pca_biplot 

filename <- paste0("./Figures/FigS4-S7-PCA/metabol/FigS7-NMR_PCA_biplot.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_biplot)


########################################## PCA on FTICR#####################
# Calculate PCA with prcomp()
pca <- prcomp(df.f, scale = TRUE, center = TRUE)
# Extract eigenvalues and variances
eigen <- get_eigenvalue(pca) # this function is from factoextra
dimensions <- c(1:dim(eigen)[1]) # this is for the plot

scree <- make_screeplot(eigen, dimensions) + ggtitle('Scree plot, PCA Class')
filename <- paste0("./Figures/PCA_meta/PCA_screeplot_FTICR.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,scree)
cumvar <-  make_cumvar(eigen, dimensions) + ggtitle('Cumulative variance plot, PCA Class')
filename <- paste0("./Figures/PCA_meta/PCA_cumulative_var_FTICR.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,cumvar)

# extract coordinates for PC1 and PC2
pca_results <- get_pca_ind(pca) #pca[["x"]]
pca_coordinates <- as.data.frame(pca_results$coord[,c(1,2)])
colnames(pca_coordinates) <- c('PC1','PC2')

# merge metadata with pca_coordinates:
pca_coordinates$SampleID <- rownames(pca_coordinates)
pca_coordinates <- pca_coordinates %>% 
  separate(., SampleID, into = c('Soil','Condition',	'Collar',	'Site',	'Time'), sep = "_", remove = T)

pca_coordinates$Condition <- factor(pca_coordinates$Condition, c("PreDrought", "Drought"))

write.csv(pca_coordinates,file=paste0("./Output/PCA_meta/PCA_individual_coordinates_FTICR.csv"),row.names=TRUE)

# prepare label for graph
pc1 <- paste0('PC1 (',round(eigen$variance.percent[1],digits=1),'%)')
pc2 <- paste0('PC2 (',round(eigen$variance.percent[2],digits=1),'%)')

# arrows
arrows <- get_arrows(pca, pca_coordinates)
write.csv(arrows,file=paste0("./Output/PCA_meta/PCA_vector_coordinates_FTICR.csv"), row.names=TRUE)

#filter arrows to top 25% loadings. Can modify "0.75" to "0.85" for 15% top loadings, etc.
top <- quantile(arrows$coord, 0.99) 
arrows.f <- arrows %>%
  subset(coord > top)
length(arrows.f$coord)

#pca
pca_plot <-  ggplot(mapping = aes(x, y)) +
  geom_point(data=pca_coordinates, aes(x=PC1, y=PC2, col=Condition, shape = Site), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= pc1, y=pc2) +
  scale_color_manual(values = colors,
                     breaks = c("PreDrought", "Drought"),
                     labels = c("Pre Drought", "Drought")) +
  scale_shape_manual(values= list_of_shapes) +
  
  theme( legend.text = element_text(size=size+3, face="bold"),
         legend.title = element_blank(),
         legend.key.size = unit(0.6, "cm"),
         legend.key.width = unit(0.6,"cm"),
         legend.position = "bottom",
         axis.title.x = element_text(size=size+3,face="bold"),
         axis.title.y = element_text(size=size+3,face="bold"),
         plot.title = element_text(size=size+3,face="bold")) 
pca_plot
filename <- paste0("./Figures/FigS4-S7-PCA/metabol/PCA_plot_FTICR.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

# pca biplot
pca_biplot <- pca_plot +
  new_scale_color() +
  geom_segment(data=arrows.f, aes(x=0, y=0, xend=xend, yend=yend),
               arrow=arrow(length = unit(0.1,"cm")), size=0.7, color = "grey") +
  geom_text_repel(data=arrows.f, aes(x=xend, y=yend), color = "black",
                  label=arrows.f$name, size=size/3, show.legend = FALSE)

pca_biplot 

filename <- paste0("./Figures/FigS4-S7-PCA/metabol/PCA_biplot_FTICR-top-one-perc.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_biplot)

########################################## PCA on FTICR classes#####################
# Calculate PCA with prcomp()
pca <- prcomp(df.c, scale = TRUE, center = TRUE)
# Extract eigenvalues and variances
eigen <- get_eigenvalue(pca) # this function is from factoextra
dimensions <- c(1:dim(eigen)[1]) # this is for the plot

scree <- make_screeplot(eigen, dimensions) + ggtitle('Scree plot, PCA Class')
filename <- paste0("./Figures/PCA_meta/PCA_screeplot_FTICR_class.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,scree)
cumvar <-  make_cumvar(eigen, dimensions) + ggtitle('Cumulative variance plot, PCA Class')
filename <- paste0("./Figures/PCA_meta/PCA_cumulative_var_FTICR_class.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,cumvar)

# extract coordinates for PC1 and PC2
pca_results <- get_pca_ind(pca) #pca[["x"]]
pca_coordinates <- as.data.frame(pca_results$coord[,c(1,2)])
colnames(pca_coordinates) <- c('PC1','PC2')

# merge metadata with pca_coordinates:
pca_coordinates$SampleID <- rownames(pca_coordinates)
pca_coordinates <- pca_coordinates %>% 
  separate(., SampleID, into = c('Soil','Condition',	'Collar',	'Site',	'Time'), sep = "_", remove = T)

pca_coordinates$Condition <- factor(pca_coordinates$Condition, c("PreDrought", "Drought"))
pca_coordinates$Time <- factor(pca_coordinates$Time, c("0", "6", "48"))

write.csv(pca_coordinates,file=paste0("./Output/PCA_meta/PCA_individual_coordinates_FTICR_class.csv"),row.names=TRUE)

# prepare label for graph
pc1 <- paste0('PC1 (',round(eigen$variance.percent[1],digits=1),'%)')
pc2 <- paste0('PC2 (',round(eigen$variance.percent[2],digits=1),'%)')

# arrows
arrows <- get_arrows(pca, pca_coordinates)
write.csv(arrows,file=paste0("./Output/PCA_meta/PCA_vector_coordinates_FTICR_class.csv"), row.names=TRUE)

#pca with shapes as time since pyruvate addition
pca_plot <-  ggplot(mapping = aes(x, y)) +
  geom_point(data=pca_coordinates, aes(x=PC1, y=PC2, col=Condition, shape = Time), 
             size=size/2, show.legend = TRUE) +
  theme_linedraw(base_size = size) + labs(x= pc1, y=pc2) +
  scale_color_manual(values = colors,
                     breaks = c("PreDrought", "Drought"),
                     labels = c("Pre Drought", "Drought")) +
  scale_shape_manual(values= list_of_shapes) +
  
  theme( legend.text = element_text(size=size+3, face="bold"),
         legend.title = element_blank(),
         legend.key.size = unit(0.6, "cm"),
         legend.key.width = unit(0.6,"cm"),
         legend.position = "bottom",
         axis.title.x = element_text(size=size+3,face="bold"),
         axis.title.y = element_text(size=size+3,face="bold"),
         plot.title = element_text(size=size+3,face="bold")) 
pca_plot
filename <- paste0("./Figures/FigS4-S7-PCA/metabol/VOC_PCA_plot_FTICR_class.png")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_plot)

# pca biplot
pca_biplot <- pca_plot +
  new_scale_color() +
  geom_segment(data=arrows, aes(x=0, y=0, xend=xend, yend=yend),
               arrow=arrow(length = unit(0.1,"cm")), size=0.7, color = "grey") +
  geom_text_repel(data=arrows, aes(x=xend, y=yend), color = "black",
                  label=arrows$name, size=size/3, show.legend = FALSE)

pca_biplot 

filename <- paste0("./Figures/FigS4-S7-PCA/metabol/FigS7-PCA_biplot_FTICR_class.pdf")
ggsave(filename,units=c('in'),width=w,height=h,dpi=res,pca_biplot)

####statistical tests#####
##NMR###
#wilcoxan

df.s <- df %>%
  rownames_to_column(., var = "sampleID") %>%
  separate(., sampleID, into = c('unique','Condition',	'Time',	'Site',	'Label'), sep = "_", remove = T) %>%
  dplyr::rename ('Oxoisocaproate' = '2.Oxoisocaproate') %>%
  dplyr::rename('Hydroxybutyrate' = '3.Hydroxybutyrate')

modelList <- list()
for(i in 6:29) {
  fmla <- formula(paste(names(df.s) [i]," ~ Condition"))
  modelList[[i]]<- wilcox.test(fmla, data = df.s, paired = FALSE)
}
modelList

sink("./Output/PCA_meta/Wilcox-results.txt")
print(modelList)
sink()

#linear mixed effects
modelList <- list()
for(i in 6:29) {
  fmla <- formula(paste(names(df.s) [i]," ~ Condition"))
  lme<- lme(fmla,
                        random = list(Site = ~1),
                        data = df.s,
                        #weights =  varIdent(form = ~1|Condition),
                        na.action=na.omit)
  modelList[[i]] <- summary(lme)
}
modelList

sink("./Output/PCA_meta/NMR-LME.txt")
print(modelList)
sink()

#permanova
dist <- dist(df)
adonis(dist~Condition, df.s)

##FTICR###
#wilcoxan

df.f.s <- df.f %>%
  rownames_to_column(., var = "sampleID") %>%
  separate(., sampleID, into = c('Soil','Condition',	'Collar',	'Site',	'Time'), sep = "_", remove = T)

#permanova
dist <- dist(df.f)
adonis(dist~Condition, df.f.s)

#linear mixed effects


##FTICR classes###
#wilcoxan

df.c.s <- df.c %>%
  rownames_to_column(., var = "sampleID") %>%
  separate(., sampleID, into = c('Soil','Condition',	'Collar',	'Site',	'Time'), sep = "_", remove = T)

modelList <- list()
for(i in 6:14) {
  fmla <- formula(paste(names(df.c.s) [i]," ~ Condition"))
  modelList[[i]]<- wilcox.test(fmla, data = df.c.s, paired = FALSE)
}
modelList

sink("./Output/PCA_meta/Wilcox-results-FTICR-class.txt")
print(modelList)
sink()

#linear mixed effects
modelList <- list()
for(i in 6:14) {
  fmla <- formula(paste(names(df.c.s) [i]," ~ Condition"))
  lme<- lme(fmla,
            random = list(Site = ~1),
            data = df.c.s,
            #weights =  varIdent(form = ~1|Condition),
            na.action=na.omit)
  modelList[[i]] <- summary(lme)
}
modelList

sink("./Output/PCA_meta/FTICR-class-LME.txt")
print(modelList)
sink()

#permanova
dist <- dist(df.c)
adonis(dist~Condition, df.c.s)
