##Linnea Honeker/
##2021-9-10
##R analysis of soil pyruvate metaT taxonomy data

library(data.table)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(ggplot2)
library(ape)
library(magrittr)
library(dplyr)
library(vegan)
library (DESeq2)
library(devtools)
library(RColorBrewer)
library(microbiome)
library(knitr)
library(plyr)

#set colors
colors_very_short <- c("black", "red", "blue", "green")
colors_short <- c("black","gray",
                 "darkred", "darksalmon", "green4",
                 "greenyellow", "lightslateblue","orange",
                 "moccasin", "lightpink", "lightsteelblue1", "navy", "blue4",
                 "darkgray", "black","brown", "cornflowerblue","darkgoldenrod",
                 "brown3","white","aliceblue","aquamarine","azure",
                 "beige","bisque","blue2","blueviolet","cyan","darkblue",
                 "chocolate","coral","coral4","aquamarine3","darkcyan","deeppink",
                 "darkred","darkslateblue")

colors_long <- c("black","gray","coral4","coral",
                "chartreuse3", "darkseagreen1","blue","lightblue",
                "yellow4", "yellow","darkorchid", "plum2",
                "darkred", "darksalmon", "green4",
                "greenyellow", "orange",
                "moccasin", "hotpink4", "lightpink", "lightblue4", "lightcyan3",  "lightslateblue", "lightsteelblue1", "navy", "blue4",
                "darkgray", "black","brown", "cornflowerblue","darkgoldenrod",
                "brown3","white","aliceblue","aquamarine","azure",
                "beige","bisque","blue2","blueviolet","cyan","darkblue",
                "chocolate","aquamarine3","darkcyan","deeppink",
                "darkred","darkslateblue")

#import biom table
otu=read.table('./Data/Taxonomy/MetaG/feature_table.csv',sep = ",", header = TRUE)

rownames(otu) <- otu[,1]
otu$Feature <- NULL

otu <- as.matrix(otu) 
is.matrix(otu)

otu <- otu_table(otu, taxa_are_rows = TRUE)

#import sample metadata
map = read.table('./Data/Taxonomy/MetaG/metadata.csv',sep = ",", header = TRUE)
row.names(map) <- map[,1]
map$SampleID <- NULL
map <- as.data.frame(map)
map <- sample_data(map)

#import taxnomy
taxa.frame = read.csv(file="./Data/Taxonomy/MetaG/metag_taxonomy.csv", header = TRUE, sep=",", 
                      blank.lines.skip = TRUE)
taxmat = as.matrix(taxa.frame)
rownames(taxmat) <- taxmat[,1]
taxmat <- taxmat[,-1]
tax=tax_table(taxmat)


#merge otu-table, sample-metadata, taxonomy
SP= merge_phyloseq(otu, tax,map)
SP

taxa(SP)

#filter unknowns
SP = subset_taxa(SP, Phylum != "unknown")


#filter to leave Bacteria Kingdom
SPb = subset_taxa(SP, Kingdom == "Bacteria" |
                    Kingdom == "Archaea")
SPb

#filter to fungal kingdom
SPf = subset_taxa(SP, Kingdom == "Eukaryota")

#taxonomy--------------------------------------------------------
#prune out phyla below 2% in each sample for all experiments
SP.rel.abun <- transform_sample_counts(SP, function(x) x/sum(x)) #convert counts to relative abundances
k.glom <- tax_glom(SP.rel.abun, taxrank = "Phylum") #agglomerate taxa to Phylum-level
k.dat <- psmelt(k.glom)
k.dat$Kingdom <- as.character(k.dat$Kingdom) #convert to character vector from factor


#prune out phyla below 2% in each sample for all experiments for bacteria/archaea
SPb.rel.abun <- transform_sample_counts(SPb, function(x) x/sum(x)) #convert counts to relative abundances
P.glom <- tax_glom(SPb.rel.abun, taxrank = "Phylum") #agglomerate taxa to Phylum-level
P.dat <- psmelt(P.glom)
P.dat$Phylum <- as.character(P.dat$Phylum) #convert to character vector from factor
P.dat$Phylum[P.dat$Abundance < 0.02] <- "Phylum < 2%"

SPb.rel.abun <- transform_sample_counts(SPb, function(x) x/sum(x)) #convert counts to relative abundances
C.glom <- tax_glom(SPb.rel.abun, taxrank = "Class") #agglomerate taxa to Phylum-level
C.dat <- psmelt(C.glom)
C.dat$Class <- as.character(C.dat$Class) #convert to character vector from factor
C.dat$Class[C.dat$Abundance < 0.02] <- "Class < 2%"

SPb.rel.abun <- transform_sample_counts(SPb, function(x) x/sum(x)) #convert counts to relative abundances
F.glom <- tax_glom(SPb.rel.abun, taxrank = "Family") #agglomerate taxa to Phylum-level
F.dat <- psmelt(F.glom)
F.dat$Family <- as.character(F.dat$Family) #convert to character vector from factor
F.dat$Family[F.dat$Abundance < 0.035] <- "Family < 3.5%"

SPb.rel.abun <- transform_sample_counts(SPb, function(x) x/sum(x)) #convert counts to relative abundances
G.glom <- tax_glom(SPb.rel.abun, taxrank = "Genus") #agglomerate taxa to Phylum-level
G.dat <- psmelt(G.glom)
G.dat$Genus <- as.character(G.dat$Genus) #convert to character vector from factor
G.dat$Genus[G.dat$Abundance < 0.035] <- "Genus < 3.5%"

#prune out phyla below 2% in each sample for all experiments for eukaryota
SPf.rel.abun <- transform_sample_counts(SPf, function(x) x/sum(x)) #convert counts to relative abundances
Pf.glom <- tax_glom(SPf.rel.abun, taxrank = "Phylum") #agglomerate taxa to Phylum-level
Pf.dat <- psmelt(Pf.glom)
Pf.dat$Phylum <- as.character(Pf.dat$Phylum) #convert to character vector from factor
Pf.dat$Phylum[Pf.dat$Abundance < 0.02] <- "Phylum < 2%"

SPf.rel.abun <- transform_sample_counts(SPf, function(x) x/sum(x)) #convert counts to relative abundances
Cf.glom <- tax_glom(SPf.rel.abun, taxrank = "Class") #agglomerate taxa to Phylum-level
Cf.dat <- psmelt(Cf.glom)
Cf.dat$Class <- as.character(Cf.dat$Class) #convert to character vector from factor
Cf.dat$Class[Cf.dat$Abundance < 0.02] <- "Class < 2%"

SPf.rel.abun <- transform_sample_counts(SPf, function(x) x/sum(x)) #convert counts to relative abundances
Ff.glom <- tax_glom(SPf.rel.abun, taxrank = "Family") #agglomerate taxa to Phylum-level
Ff.dat <- psmelt(Ff.glom)
Ff.dat$Family <- as.character(Ff.dat$Family) #convert to character vector from factor
Ff.dat$Family[Ff.dat$Abundance < 0.035] <- "Family < 3.5%"

SPf.rel.abun <- transform_sample_counts(SPf, function(x) x/sum(x)) #convert counts to relative abundances
Gf.glom <- tax_glom(SPf.rel.abun, taxrank = "Genus") #agglomerate taxa to Phylum-level
Gf.dat <- psmelt(Gf.glom)
Gf.dat$Genus <- as.character(Gf.dat$Genus) #convert to character vector from factor
Gf.dat$Genus[Gf.dat$Abundance < 0.035] <- "Genus < 3.5%"

#####################taxonomy plots- condition and time for eukaryota#########################

#######plot Phylum######################################
Pf.dat$Condition <- factor(Pf.dat$Condition, c("PreDrought", "Drought"))
Pf.dat$Time <- factor(Pf.dat$Time, c("0hr", "6hr", "48hr"))


p <- ggplot(Pf.dat, aes(x=Time, y=Abundance, fill=Phylum)) 

p_barplot <- p + geom_bar(aes(), stat = "identity", position = "fill") +
  facet_wrap(~ Condition, scales="free") +
  scale_fill_manual(values=colors_short) +
  theme_bw() +
  theme(text = element_text(size = 11),
        axis.text.x=element_text(hjust =0.5, vjust=0.2),
        legend.title = element_blank()
  )

filename <- paste0("./Figures/FigS3-Taxonomy/FigS3-metaG-Phylum.png")
ggsave(filename,width=5,height=5,dpi=1000,p_barplot)

#######plot Classes######################################
C.dat$Condition <- factor(C.dat$Condition, c("PreDrought", "Drought"))
C.dat$Time <- factor(C.dat$Time, c("0hr", "6hr", "48hr"))


p <- ggplot(C.dat, aes(x=Time, y=Abundance, fill=Class)) 

c_barplot <- p + geom_bar(aes(), stat = "identity", position = "fill") +
  facet_wrap(~ Condition, scales="free") +
  scale_fill_manual(values=colors_long) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x=element_text(size=14, angle=90,hjust =1, vjust=0.2),
        #legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
  )

filename <- paste0("./Figures/FigS3-Taxonomy/metaG-class.png")
ggsave(filename,width=6,height=6,dpi=1000,c_barplot)

#######plot family######################################
F.dat$Condition <- factor(F.dat$Condition, c("PreDrought", "Drought"))
F.dat$Time <- factor(F.dat$Time, c("0hr", "6hr", "48hr"))


p <- ggplot(F.dat, aes(x=Time, y=Abundance, fill=Family)) 

f_barplot <- p + geom_bar(aes(), stat = "identity", position = "fill") +
  facet_wrap(~ Condition, scales="free") +
  scale_fill_manual(values=colors_long) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x=element_text(size=14, angle=90,hjust =1, vjust=0.2),
        #legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
  )
filename <- paste0("./Figures/FigS3-Taxonomy/metaG-family.png")
ggsave(filename,width=6,height=6,dpi=1000,f_barplot)

#######plot Phylum######################################
P.dat$Condition <- factor(P.dat$Condition, c("PreDrought", "Drought"))
P.dat$Time <- factor(P.dat$Time, c("0hr", "6hr", "48hr"))


p <- ggplot(P.dat, aes(x=Time, y=Abundance, fill=Phylum)) 

p_barplot <- p + geom_bar(aes(), stat = "identity", position = "fill") +
  facet_wrap(~ Condition, scales="free") +
  scale_fill_manual(values=colors_short) +
  theme_bw() +
  theme(text = element_text(size = 11),
        axis.text.x=element_text(hjust =0.5, vjust=0.2),
        legend.title = element_blank()
  )

filename <- paste0("./Figures/FigS3-Taxonomy/FigS3-metaG-Phylum.png")
ggsave(filename,width=5,height=5,dpi=1000,p_barplot)

#######plot Classes######################################
C.dat$Condition <- factor(C.dat$Condition, c("PreDrought", "Drought"))
C.dat$Time <- factor(C.dat$Time, c("0hr", "6hr", "48hr"))


p <- ggplot(C.dat, aes(x=Time, y=Abundance, fill=Class)) 

c_barplot <- p + geom_bar(aes(), stat = "identity", position = "fill") +
  facet_wrap(~ Condition, scales="free") +
  scale_fill_manual(values=colors_long) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x=element_text(size=14, angle=90,hjust =1, vjust=0.2),
        #legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
  )

filename <- paste0("./Figures/FigS3-Taxonomy/metaG-class.png")
ggsave(filename,width=6,height=6,dpi=1000,c_barplot)

#######plot family######################################
F.dat$Condition <- factor(F.dat$Condition, c("PreDrought", "Drought"))
F.dat$Time <- factor(F.dat$Time, c("0hr", "6hr", "48hr"))


p <- ggplot(F.dat, aes(x=Time, y=Abundance, fill=Family)) 

f_barplot <- p + geom_bar(aes(), stat = "identity", position = "fill") +
  facet_wrap(~ Condition, scales="free") +
  scale_fill_manual(values=colors_long) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x=element_text(size=14, angle=90,hjust =1, vjust=0.2),
        #legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
  )
filename <- paste0("./Figures/FigS3-Taxonomy/metaG-family.png")
ggsave(filename,width=6,height=6,dpi=1000,f_barplot)

#####################taxonomy plots- condition and site#########################
#######plot kingdom######################################
k.dat$Condition <- factor(k.dat$Condition, c("PreDrought", "Drought"))

p <- ggplot(k.dat, aes(x=Site, y=Abundance, fill=Kingdom)) 

K_barplot <- p + geom_bar(aes(), stat = "identity", position = "fill") +
  facet_wrap(~ Condition, scales="free") +
  scale_fill_manual(values=colors_very_short) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x=element_text(size=14, angle=90,hjust =1, vjust=0.2),
        #legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
  )

filename <- paste0("./Figures/FigS3-Taxonomy/metaG-Kingdom-site.png")
ggsave(filename,width=6,height=6,dpi=1000,K_barplot)

#######plot Phylum######################################
P.dat$Condition <- factor(P.dat$Condition, c("PreDrought", "Drought"))


p <- ggplot(P.dat, aes(x=Site, y=Abundance, fill=Phylum)) 

p_barplot <- p + geom_bar(aes(), stat = "identity", position = "fill") +
  facet_wrap(~ Condition, scales="free") +
  scale_fill_manual(values=colors_short) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x=element_text(size=14, angle=90,hjust =1, vjust=0.2),
        #legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
  )

filename <- paste0("./Figures/FigS3-Taxonomy/metaG-Phylum-site.png")
ggsave(filename,width=6,height=6,dpi=1000,p_barplot)

#######plot Classes######################################
C.dat$Condition <- factor(C.dat$Condition, c("PreDrought", "Drought"))


p <- ggplot(C.dat, aes(x=Site, y=Abundance, fill=Class)) 

c_barplot <- p + geom_bar(aes(), stat = "identity", position = "fill") +
  facet_wrap(~ Condition, scales="free") +
  scale_fill_manual(values=colors_long) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x=element_text(size=14, angle=90,hjust =1, vjust=0.2),
        #legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
  )

filename <- paste0("./Figures/FigS3-Taxonomy/metaG-class-site.png")
ggsave(filename,width=6,height=6,dpi=1000,c_barplot)

#######plot family######################################
F.dat$Condition <- factor(F.dat$Condition, c("PreDrought", "Drought"))


p <- ggplot(F.dat, aes(x=Site, y=Abundance, fill=Family)) 

f_barplot <- p + geom_bar(aes(), stat = "identity", position = "fill") +
  facet_wrap(~ Condition, scales="free") +
  scale_fill_manual(values=colors_long) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x=element_text(size=14, angle=90,hjust =1, vjust=0.2),
        #legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
  )

filename <- paste0("./Figures/FigS3-Taxonomy/metaG-family-site.png")
ggsave(filename,width=6,height=6,dpi=1000,f_barplot)


#######beta diversity ########
SPb.ra = rarefy_even_depth(SPb,(min(sample_sums(SPb))))

sample_data(SPb.ra)$Condition <- factor(sample_data(SPb.ra)$Condition, levels = c("PreDrought","Drought"))
sample_data(SPb.ra)$Time <- factor(sample_data(SPb.ra)$Time, levels = c("0hr","6hr", "48hr"))
sample_data(SPb.ra)$Condition_time <- factor(sample_data(SPb.ra)$Condition_time, 
                                             levels = c("PreDrought_0hr","PreDrought_6hr", "PreDrought_48hr",
                                                        "Drought_0hr", "Drought_6hr", "Drought_48hr"))

my.ord <- ordinate(
  physeq = SPb.ra,
  method = "PCoA",
  distance = "bray")
  #weighted=TRUE)

#condition and time
ord.plot <- (
  plot_ordination(
    physeq = SPb.ra,
    ordination = my.ord,
    color = "Condition",
    #type = "split",
    shape = "Time",
    title = "PCoA of soil pyruate  (bray-curtis)"
  ) +
    geom_point(aes(color = Condition), size = 4) )
  
filename <- paste0("./Figures/FigS3-Taxonomy/metaG-PcoA-bray.png")
ggsave(filename,width=10,height=5,dpi=1000,ord.plot)

#condition and site
ord.plot <- (
  plot_ordination(
    physeq = SPb.ra,
    ordination = my.ord,
    color = "Condition",
    #type = "split",
    shape = "Site",
    title = "PCoA of soil pyruate  (bray-curtis)"
  ) +
    geom_point(aes(color = Condition), size = 4) )#+

filename <- paste0("./Figures/FigS3-Taxonomy/metaG-PcoA-bray-site.png")
ggsave(filename,width=10,height=5,dpi=1000,ord.plot)

#######alpha diversity#######
#condition and time
#define and plot alpha diversity measurements on unrarefied samples
my_colors_st = c("#B6D7A8", "#6AA84F", "#274E13", "#EA9999", "#CC0000", "#660000")
alpha_meas = c("Observed", "Shannon")
p <- plot_richness(SPb.ra, "Condition_time",scales = "free_y", measures="Observed")

#Add ggplot2 layer-boxplots
p1 <- p + scale_color_manual(values = my_colors_st) +
  geom_boxplot(data=p$data, aes(x=Condition_time, y=value, fill=Condition_time), 
               alpha=20, size = 0.5) +
  scale_fill_manual(values = my_colors_st) +
  scale_x_discrete(name = element_blank(), labels = c("0hr", "6hr", "48hr", "0hr", "6hr", "48hr")) +
  theme_bw() +
  theme(text = element_text(size = 11),
       legend.position = "none",
        axis.text.x=element_text(size=11,hjust =0.5, vjust=0.8))
filename <- paste0("./Figures/FigS3-Taxonomy/metaG-alpha-obderved.png")
ggsave(filename,width=3,height=4,dpi=1000,p1)

#condition and site
#define and plot alpha diversity measurements on unrarefied samples
p <- plot_richness(SPb.ra, "Site", measures="Observed")

#Add ggplot2 layer-boxplots
p1 <- p + scale_color_manual(values = c("orange", "magenta", "turquoise")) +
  geom_boxplot(data=p$data, aes(x=Site, y=value, fill=Site), 
               alpha=20, size = 0.5) +
  scale_fill_manual(values = c("orange", "magenta", "turquoise")) +
  scale_x_discrete(name = element_blank()) +
  theme_bw() +
  theme(text = element_text(size = 11),
        legend.position = "none",
        axis.text.x=element_text(size=11,hjust =0.5, vjust=0.8))
filename <- paste0("./Figures/FigS3-Taxonomy/metaG-alpha-obderved-site.png")
ggsave(filename,width=3,height=4,dpi=1000,p1)

#determine significance in alpha diversity
alpha = estimate_richness(SPb.ra)

pairwise.wilcox.test(alpha$Observed, sample_data(SPb)$Condition_time, p.adjust.method = "none")


#data:  alpha$Observed and sample_data(SPb)$Condition_time 
#
#Drought_0hr Drought_48hr Drought_6hr PreDrought_0hr PreDrought_48hr
#Drought_48hr    0.937       -            -           -              -              
#  Drought_6hr     0.394       0.699        -           -              -              
#  PreDrought_0hr  0.589       0.818        0.394       -              -              
#  PreDrought_48hr 0.310       0.699        0.818       1.000          -              
#  PreDrought_6hr  0.240       0.394        0.149       0.093          0.818          


#P value adjustment method: none 

pairwise.wilcox.test(alpha$Observed, sample_data(SPb)$Site, p.adjust.method = "none")
#data:  alpha$Observed and sample_data(SPb)$Condition 

#Drought
#PreDrought 0.05   

#P value adjustment method: none 
