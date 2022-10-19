#---
#  title: "glm"
#author: "Linnea Honeker"
#date: "10/10/2022"
#---

# Linking gene expression and voc fluxes

library(tidyverse)
library(PCAtest)
library(MASS)
library(factoextra)
library(corrplot)

source('./RScripts/GeneFlux/extra_functions.R')

# PCA on full KO list--------------------------------------------------------------------
#####acetaet###############
###subset to only acetate cycling genes. KO_interest object created in glm.Rmd
KO_acetate <- KO_interest %>%
  filter(Compound == "Acetate") 

KO_acetate_list <- KO_acetate$KO
flux_list <- c("flux.acetate", "flux.acetone", "flux.diacetyl", "Condition", "Site")

idx <- match(c(KO_acetate_list,flux_list), names(merged))
idx <- idx[!is.na(idx)]


sub_acetate <- merged[,idx]

sub_acetate2 <- sub_acetate %>%
  mutate(tmt = ifelse(Condition == "drought", 1, 0)) %>%
  mutate(site = ifelse(Site == "Site3", 1, 0)) %>%
  dplyr::select(-Condition, -Site)


### test for significance of loadings of acetate correlation matrix - spearman using pca test ###
cor_acetate <- cor(sub_acetate2, method = c("spearman"))
mu <- rep(0,32)
cor_acetate <- as.matrix(cor_acetate)

ex025 <- mvrnorm(18, mu = mu, Sigma = cor_acetate)
result<-PCAtest(ex025, 1000, 1000, 0.05, varcorr=FALSE, counter=FALSE, plot=TRUE)

sig.load <- c("K00128", "K01067", "K01512", "K01895", "K00138", "K02535", "K18118", "K01739", 
              "K12339", "K19670", "K01438", "K00156",  "K01026", "K01040", "K01060",
            "K16263", "flux.acetate", "flux.acetone",  "tmt", "site")

### PCA on acetate #######
PCAdata.acetate = prcomp(na.omit(sub_acetate2), scale = TRUE)
PCAdata.acetate

pca_results <- get_pca_ind(PCAdata.acetate) #pca[["x"]]
var <- get_pca_var(PCAdata.acetate)

pca_coordinates <- as.data.frame(pca_results$coord[,c(1,2,3,4,5,6)])
colnames(pca_coordinates) <- c('PC1','PC2','PC3','PC4','PC5','PC6')

# bar plot showing cos2 - quality of represention for first 2 axes
fviz_cos2(PCAdata.acetate, choice = "var", axes = 1:3)

# correlation plost showing correlations between variables and first 2 axes
corrplot(var$cor[,1:3], is.corr=FALSE)


fviz_pca_var(PCAdata.acetate, axes = c(1, 2),# col.ind = PCAdata3$N.addition,
             #col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE,
             select.var = list(name=sig.load)
)
fviz_pca_ind(PCAdata.acetate,
             geom = c("point"),
             axes = c(1, 2),
             col.ind = as.factor(sub_acetate$Condition), # color by groups
             #palette = c("#00AFBB",  "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             #ellipse.type = "confidence",
             ellipse.level=0.95,
             legend.title = "Condition",
             repel = TRUE,
             select.var = list(name=sig.load)
)

fviz_pca_biplot(PCAdata.acetate,
             geom = c("point"),
             axes = c(1, 2),
             col.ind = as.factor(sub_acetate$Site), # color by groups
             #palette = c("#00AFBB",  "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             #ellipse.type = "confidence",
             ellipse.level=0.95,
             legend.title = "Condition",
             repel = TRUE,
             select.var = list(name=sig.load)
)

# PCA on original (subset) KO list
#####acetaet###############
###subset to only acetate cycling genes. KO_interest object created in glm.Rmd
KO_acetate <- KO_interest_subset %>%
  filter(VOC == "Acetate") 

KO_acetate_list <- KO_acetate$KO
flux_list <- c("flux.acetate", "flux.acetone", "flux.diacetyl", "Condition", "Site")

idx <- match(c(KO_acetate_list,flux_list), names(merged_subset))
idx <- idx[!is.na(idx)]


sub_acetate <- merged_subset[,idx]

sub_acetate2 <- sub_acetate %>%
  mutate(tmt = ifelse(Condition == "drought", 1, 0)) %>%
  dplyr::select(-Condition, -Site)


### test for significance of loadings of acetate correlation matrix - spearman using pca test ###
cor_acetate <- cor(sub_acetate2, method = c("spearman"))
mu <- rep(0,12)
cor_acetate <- as.matrix(cor_acetate)

ex025 <- mvrnorm(18, mu = mu, Sigma = cor_acetate)
result<-PCAtest(ex025, 1000, 1000, 0.05, varcorr=FALSE, counter=FALSE, plot=TRUE)

sig.load <- c("K00128", "K01067", "K01512", "K01895", "K00138", "K02535", "K18118", "K01739", 
              "K12339", "K19670", "K01438", "K00156",  "K01026", "K01040", "K01060",
              "K16263", "flux.acetate", "flux.acetone",  "tmt")

### PCA on acetate #######
PCAdata.acetate = prcomp(na.omit(sub_acetate2), scale = TRUE)
PCAdata.acetate

pca_results <- get_pca_ind(PCAdata.acetate) #pca[["x"]]
var <- get_pca_var(PCAdata.acetate)

pca_coordinates <- as.data.frame(pca_results$coord[,c(1,2,3,4,5,6)])
colnames(pca_coordinates) <- c('PC1','PC2','PC3','PC4','PC5','PC6')

# bar plot showing cos2 - quality of represention for first 2 axes
fviz_cos2(PCAdata.acetate, choice = "var", axes = 1:2)

# correlation plost showing correlations between variables and first 2 axes
corrplot(var$cor[,1:2], is.corr=FALSE)


fviz_pca_var(PCAdata.acetate, axes = c(1, 2),# col.ind = PCAdata3$N.addition,
             #col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE
)

fviz_pca_ind(PCAdata.acetate,
             geom = c("point"),
             axes = c(1, 2),
             col.ind = as.factor(sub_acetate$Condition), # color by groups
             #palette = c("#00AFBB",  "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             #ellipse.type = "confidence",
             ellipse.level=0.95,
             legend.title = "Condition",
             repel = TRUE
)

fviz_pca_biplot(PCAdata.acetate,
                geom = c("point"),
                axes = c(1, 2),
                col.ind = as.factor(sub_acetate$Condition), # color by groups
                #palette = c("#00AFBB",  "#FC4E07"),
                addEllipses = TRUE, # Concentration ellipses
                #ellipse.type = "confidence",
                ellipse.level=0.95,
                legend.title = "Condition",
                repel = TRUE
)




######Acetone########
#subset to only acetone cycling genes
KO_acetone <- KO_interest %>%
  filter(Compound == "Acetone") 

KO_acetone_list <- KO_acetone$KO
flux_list <- c("flux.acetate", "flux.acetone", "flux.diacetyl", "Condition")

idx <- match(c(KO_acetone_list,flux_list), names(merged))
idx <- idx[!is.na(idx)]


sub_acetone <- merged[,idx]

sub_acetone2 <- sub_acetone %>%
  mutate(tmt = ifelse(Condition == "drought", 1, 0)) %>%
  dplyr::select(-Condition)


### test for significance of loadings of acetone correlation matrix - spearman using pca test ###
cor_acetone <- cor(sub_acetone2, method = c("spearman"))
mu <- rep(0,11)
cor_acetate <- as.matrix(cor_acetate)

ex025 <- mvrnorm(18, mu = mu, Sigma = cor_acetone)
result<-PCAtest(ex025, 100, 100, 0.05, varcorr=FALSE, counter=FALSE, plot=TRUE)


###PCA on acetone ######
PCAdata.acetone = prcomp(na.omit(sub_acetone), scale = TRUE)
PCAdata.acetone

pca_results <- get_pca_ind(PCAdata.acetone) #pca[["x"]]
var <- get_pca_var(PCAdata.acetone)


pca_coordinates <- as.data.frame(pca_results$coord[,c(1,2,3,4,5,6)])
colnames(pca_coordinates) <- c('PC1','PC2','PC3','PC4','PC5','PC6')



fviz_cos2(PCAdata.acetone, choice = "var", axes = 1:2)
corrplot(var$cor[,1:8], is.corr=TRUE)


fviz_pca_var(PCAdata.acetone, axes = c(1, 2),# col.ind = PCAdata3$N.addition,
             #col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE#,
             #select.var = list(contrib = 20)
)
fviz_pca_ind(PCAdata.acetone,
             geom = c("point"),
             axes = c(1, 2),
             col.ind = as.factor(data$Condition), # color by groups
             #palette = c("#00AFBB",  "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             #ellipse.type = "confidence",
             ellipse.level=0.95,
             legend.title = "Day of Year",
             repel = TRUE
)

fviz_pca_biplot(PCAdata.acetone,
                geom = c("point"),
                axes = c(1, 2),
                col.ind = as.factor(sub_acetone$Condition), # color by groups
                #palette = c("#00AFBB",  "#FC4E07"),
                addEllipses = TRUE, # Concentration ellipses
                #ellipse.type = "confidence",
                ellipse.level=0.95,
                legend.title = "Day of Year",
                repel = TRUE
)


PCAdata.acetone = princomp(na.omit(sub.acetone), cor = TRUE, scores = TRUE)
PCAdata.acetone

fviz_pca_var(PCAdata.acetone, axes = c(1, 2),# col.ind = PCAdata3$N.addition,
             #col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)
fviz_pca_ind(PCAdata.acetone,
             geom = c("point"),
             axes = c(1, 2),
             col.ind = as.factor(data$Condition), # color by groups
             #palette = c("#00AFBB",  "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             #ellipse.type = "confidence",
             ellipse.level=0.95,
             legend.title = "Day of Year",
             repel = TRUE
)

fviz_pca_biplot(PCAdata.acetone,
                geom = c("point"),
                axes = c(1, 2),
                col.ind = as.factor(data$Condition), # color by groups
                #palette = c("#00AFBB",  "#FC4E07"),
                addEllipses = TRUE, # Concentration ellipses
                #ellipse.type = "confidence",
                ellipse.level=0.95,
                legend.title = "Day of Year",
                repel = TRUE
)



# All genes of interest
KO_list <- KO_interest$KO
flux_list <- c("log.flux.acetate", "log.flux.acetone", "log.flux.diacetyl", "Condition")

idx <- match(c(KO_list,flux_list), names(merged))
idx <- idx[!is.na(idx)]


sub_merged <- merged[,idx]

sub_merged <- sub_merged %>%
  mutate(tmt = ifelse(Condition == "drought", 1, 0)) %>%
  dplyr::select(-Condition)


PCAdata = prcomp(na.omit(sub_merged), scale = TRUE)
PCAdata

pca_results <- get_pca_ind(PCAdata) #pca[["x"]]
var <- get_pca_var(PCAdata)


pca_coordinates <- as.data.frame(pca_results$coord[,c(1,2,3,4,5,6)])
colnames(pca_coordinates) <- c('PC1','PC2','PC3','PC4','PC5','PC6')



fviz_cos2(PCAdata, choice = "var", axes = 1:2)
corrplot(var$cos2[,1:8], is.corr=FALSE)


fviz_pca_var(PCAdata, axes = c(1, 2),# col.ind = PCAdata3$N.addition,
             #col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE,
             select.var = list(contrib = 20)
)
fviz_pca_ind(PCAdata,
             geom = c("point"),
             axes = c(1, 2),
             col.ind = as.factor(data$Condition), # color by groups
             #palette = c("#00AFBB",  "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             #ellipse.type = "confidence",
             ellipse.level=0.95,
             legend.title = "Day of Year",
             repel = TRUE
)

fviz_pca_biplot(PCAdata,
                geom = c("point"),
                axes = c(1, 2),
                col.ind = as.factor(data$Condition), # color by groups
                #palette = c("#00AFBB",  "#FC4E07"),
                addEllipses = TRUE, # Concentration ellipses
                #ellipse.type = "confidence",
                ellipse.level=0.95,
                legend.title = "Day of Year",
                repel = TRUE,
                select.var= list(contrib=20)
)

# PCA on original (subset) KO list
#####all VOCs###############
###subset to only acetate cycling genes. KO_interest object created in glm.Rmd
KO <- KO_interest_subset %>%
  filter(VOC == "Acetate" |
           VOC == "Acetone" |
           VOC == "Diacetyl")

KO_list <- KO$KO
flux_list <- c("flux.acetate", "flux.acetone", "flux.diacetyl", "Condition", "Site")

idx <- match(c(KO_list,flux_list), names(merged_subset))
idx <- idx[!is.na(idx)]


sub <- merged_subset[,idx]

sub2<- sub %>%
  mutate(tmt = ifelse(Condition == "drought", 1, 0)) %>%
  dplyr::select(-Condition, -Site)


### test for significance of loadings of acetate correlation matrix - spearman using pca test ###
cor <- cor(sub2, method = c("spearman"))
mu <- rep(0,21)
cor<- as.matrix(cor)

ex025 <- mvrnorm(18, mu = mu, Sigma = cor)
result<-PCAtest(ex025, 1000, 1000, 0.05, varcorr=FALSE, counter=FALSE, plot=TRUE)

sig.load <- c("K00128", "K01067", "K01512", "K01895", "K00138", "K02535", "K18118", "K01739", 
              "K12339", "K19670", "K01438", "K00156",  "K01026", "K01040", "K01060",
              "K16263", "flux.acetate", "flux.acetone",  "tmt")

### PCA on acetate #######
PCAdata = prcomp(na.omit(sub2), scale = TRUE)
PCAdata

pca_results <- get_pca_ind(PCAdata) #pca[["x"]]
var <- get_pca_var(PCAdata)

pca_coordinates <- as.data.frame(pca_results$coord[,c(1,2,3,4,5,6)])
colnames(pca_coordinates) <- c('PC1','PC2','PC3','PC4','PC5','PC6')

# bar plot showing cos2 - quality of represention for first 2 axes
fviz_cos2(PCAdata, choice = "var", axes = 1:2)

# correlation plost showing correlations between variables and first 2 axes
corrplot(var$cor[,1:3], is.corr=FALSE)


fviz_pca_var(PCAdata, axes = c(1, 2),# col.ind = PCAdata3$N.addition,
             #col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE
)

fviz_pca_ind(PCAdata,
             geom = c("point"),
             axes = c(1, 2),
             col.ind = as.factor(sub_acetate$Condition), # color by groups
             #palette = c("#00AFBB",  "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             #ellipse.type = "confidence",
             ellipse.level=0.95,
             legend.title = "Condition",
             repel = TRUE
)

fviz_pca_biplot(PCAdata,
                geom = c("point"),
                axes = c(1, 2),
                col.ind = as.factor(sub_acetate$Condition), # color by groups
                #palette = c("#00AFBB",  "#FC4E07"),
                addEllipses = TRUE, # Concentration ellipses
                #ellipse.type = "confidence",
                ellipse.level=0.95,
                legend.title = "Condition",
                repel = TRUE
)

