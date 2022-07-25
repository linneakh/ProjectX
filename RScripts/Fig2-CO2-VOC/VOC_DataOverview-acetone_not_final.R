# Visualize Pyruvate voc Fluxes

# Johannes Ingrisch, modified by Linnea Honeker
# 2021-08-26, modified 4/1/22

# load packages
library(tidyverse)
library(plotly)
library(ggthemes)
library(patchwork)
library(dplyr)
library(data.table)
library(stringr)
library(flux)


theme_custom <- theme_few()

theme_set(theme_custom)

colors = c("#2ECC71" ,"#F39C12")

#Load 12C/(12C + 13C) data
data <- read.csv("./Data/CO2-VOCs/voc_raw/acetone-for-ggplot.csv") 
#####C13/(C12 + C13) #######
#re order factors for condition

data$Condition <- factor(data$Condition, c(
  "pre-drought", "drought"))

# filter out times outside of -12 hour to 48 hour

data_sub <- data %>%
  filter(Time > -0.3 & Time <2 )

# Isotope Sig acetone
filename=paste("./Figures/CO2-VOCs/acetone.png", sep = "")
png(filename ,width=3, height=4, unit='in', res = 1000)

data_sub %>%
  filter(Type == "sample") %>% 
  mutate(Time_hours = Time*24) %>%
  mutate(Site = case_when(
    startsWith(Sample, "S1") ~ "Site1",
    startsWith(Sample, "S2") ~ "Site2",
    startsWith(Sample, "S3") ~ "Site3"
  )) %>%
  ggplot(aes(x = Time_hours, y = Flux, color = Condition)) +
  facet_grid(Site ~ Pyruv) + 
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(aes(group = Condition), span = 0.5) +
  scale_color_manual(values = colors,
                     breaks = c("pre-drought", "drought"),
                     labels = c("Pre Drought", "Drought")) + 
  labs(x = "Hours after labeling", y = "13C/(12C + 13C) flux") +
  ylim(-150, 300) +
  theme(text = element_text(size = 12,
                           family = "Arial",
                           color = "black"),
    legend.position = "bottom",
    legend.title = element_blank())
dev.off()




# Plotwise
data_sub %>%
  filter(Type == "sample") %>% 
  mutate(Time_hours = Time*24) %>%
  filter(Pyruv == "C1") %>%
  ggplot(aes(x = Time_hours, y = Flux, color = Condition)) +
  facet_wrap(. ~ Sample, ncol = 5) + 
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  geom_point(alpha = 0.3) +
  geom_smooth(aes(group =  Condition), span = 0.5) +
  scale_color_manual(values = colors,
                     breaks = c("pre-drought", "drought"),
                     labels = c("Pre Drought", "Drought")) + 
  labs(x = "Hours after labeling", y = "13C/(12C + 13C) flux", subtitle = "C1 pyurvate") +
  theme(legend.position = "bottom")
ggsave("./Figures/CO2-VOCs/01_Pyruvate_acetone_C1_Plotwise.png")

data_sub %>%
  filter(Type == "sample") %>% 
  mutate(Time_hours = Time*24) %>%
  filter(Pyruv == "C2") %>%
  ggplot(aes(x = Time_hours, y = Flux, color = Condition)) +
  facet_wrap(. ~ Sample, ncol = 5) + 
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  geom_point(alpha = 0.3) +
  geom_smooth(aes(group = Condition), span = 0.5) +
  scale_color_manual(values = colors,
                     breaks = c("pre-drought", "drought"),
                     labels = c("Pre Drought", "Drought")) + 
  labs(x = "Hours after labeling", y = "13C/(12C + 13C) flux", subtitle = "C2 pyurvate") +
  theme(legend.position = "bottom")
ggsave("./Figures/CO2-VOCs/01_Pyruvate_acetone_C2_Plotwise.png")

# Isotope Sig - boxplots
filename=paste("./Figures/CO2-VOCs/aacetone-boxplot-site.png", sep = "")
png(filename ,width=4, height=4, unit='in', res = 1000)

data_sub %>%
  mutate(Site = case_when(
    startsWith(Sample, "S1") ~ "Site1",
    startsWith(Sample, "S2") ~ "Site2",
    startsWith(Sample, "S3") ~ "Site3")) %>%
  #filter(Type == "sample") %>% 
  mutate(label_state = ifelse(Time >0 & Time < 2, "post", "pre")) %>% 
  mutate(Time_hours = Time*24) %>%
  filter(label_state == "post") %>%
  ggplot(aes(x = Condition, y = Flux, fill = Condition)) +
  facet_grid(Pyruv ~ Site) + 
  #facet_grid(~ Pyruv) +
   geom_boxplot() +
  scale_x_discrete(labels = c("Pre Drought", "Drought")) +
  scale_fill_manual(values = colors,
                     breaks = c("pre-drought", "drought"),
                     labels = c("Pre Drought", "Drought")) + 
  labs(y = "13C/(12C+13C) flux" ) +
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        legend.title = element_blank())
dev.off()

# Isotope Sig - boxplots - C1
filename=paste("./Figures/CO2-VOCs/aacetone-boxplot-C1.png", sep = "")
png(filename ,width=2, height=4, unit='in', res = 1000)

data_sub %>%
  filter(Type == "sample") %>% 
  mutate(label_state = ifelse(Time >0 & Time < 2, "post", "pre")) %>% 
  mutate(Time_hours = Time*24) %>%
  filter(label_state == "post") %>%
  filter(Pyruv == "C1") %>%
  ggplot(aes(x = Condition, y = Flux, fill = Condition)) +
  facet_grid(. ~ Pyruv) + 
  geom_boxplot() +
  scale_x_discrete(labels = c("Pre Drought", "Drought")) +
  scale_fill_manual(values = colors,
                    breaks = c("pre-drought", "drought"),
                    labels = c("Pre Drought", "Drought")) + 
  ylim(-150, 300) +
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank())
dev.off()

# Isotope Sig - boxplots - C2
filename=paste("./Figures/CO2-VOCs/aacetone-boxplot-C2.png", sep = "")
png(filename ,width=2, height=4, unit='in', res = 1000)

data_sub %>%
  filter(Type == "sample") %>% 
  mutate(label_state = ifelse(Time >0 & Time < 2, "post", "pre")) %>% 
  mutate(Time_hours = Time*24) %>%
  filter(label_state == "post") %>%
  filter(Pyruv == "C2") %>%
  ggplot(aes(x = Condition, y = Flux, fill = Condition)) +
  facet_grid(. ~ Pyruv) + 
  geom_boxplot() +
  scale_x_discrete(labels = c("Pre Drought", "Drought")) +
  scale_fill_manual(values = colors,
                    breaks = c("pre-drought", "drought"),
                    labels = c("Pre Drought", "Drought")) + 
  ylim(-150, 300) +
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank())
dev.off()

####statistical analysis
#Flux - Condition, in C1 and C2 separated
data_c1 <- data_sub %>%
 filter(Pyruv == "C1")

data_c2 <- data_sub %>%
  filter(Pyruv == "C2")


wilcox.test(Flux ~ Condition, data = data_c1)
#data:  Flux by Condition
#W = 98902, p-value = 0.592
#alternative hypothesis: true location shift is not equal to 0


wilcox.test(Flux ~ Condition, data = data_c2)
#data:  Flux by Condition
#W = 192527, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0

########C13 data###############
#Load 13C data
data.13.P <- read.csv("./Data/CO2-VOCs/voc_raw/acetone-pre-drought-13C.csv") 
data.13.D <- read.csv("./Data/CO2-VOCs/voc_raw/acetone-drought-13C.csv") 

# process pre-drought dataset
# gather dataframe 
data.13.P.long <- data.13.P %>%
  pivot_longer(cols = contains("P"),
              names_to = c("Label","Time"),
              names_sep = "_",
              values_to= "Time_Flux") %>%
  drop_na(Time_Flux) %>%
  mutate(Time = ifelse(is.na(Time), "flux", Time))

# since flux and time are now in single column, with each pairs of rows containing info for single sample/timepoint,
# here we will loop through each pair of rows and spread the data so flux and time become separate rows

df.p.0 = data.frame(matrix (NA, ncol = 3, nrow = 1))

for (i in 1:nrow(data.13.P.long)) {
  if (isTRUE(data.13.P.long[i,1] == data.13.P.long[(i+1),1]) == TRUE) 
  {df.p.0[i,] <- spread(data.13.P.long[i:(i+1),], Time, Time_Flux)[1,]
  print (i)}
}

# erase rows with NA and name columns
df.p <- df.p.0 %>%
  rename(Label = X1, Flux = X2, Time = X3) %>%
  filter(Label != "NA")

#Add column for pyruvate label and add column for pre-drought
df.p.f <- df.p %>% 
  mutate_at("Label", str_replace, ".C", "_C") %>%
  separate(Label, into = c("Label", "Pyr"), sep = "_") %>%
  mutate(Condition = "Pre-Drought")

### process drought dataset
# gather dataframe 
data.13.D.long <- data.13.D %>%
  pivot_longer(cols = contains("P"),
               names_to = c("Label","Time"),
               names_sep = "_",
               values_to= "Time_Flux") %>%
  drop_na(Time_Flux) %>%
  mutate(Time = ifelse(is.na(Time), "flux", Time))

data.13.D.long$Label
# since flux and time are now in single column, with each pairs of rows containing info for single sample/timepoint,
# here we will loop through each pair of rows and spread the data so flux and time become separate rows
df.d.0 = data.frame(matrix (NA, ncol = 3, nrow = 1))

for (i in 1:nrow(data.13.D.long)) {
  if (isTRUE(data.13.D.long[i,1] == data.13.D.long[(i+1),1]) == TRUE) 
  {df.d.0[i,] <- spread(data.13.D.long[i:(i+1),], Time, Time_Flux)[1,]
  print (i)}
}

# erase rows with NA and name columns
df.d <- df.d.0 %>%
  rename(Label = X1, Flux = X2, Time = X3) %>%
  filter(Label != "NA")

#Add column for pyruvate label and add column for pre-drought
df.d.f <- df.d %>% 
  mutate_at("Label", str_replace, ".C", "_C") %>%
  mutate_at("Label", str_replace, ".W", "_W") %>%
  separate(Label, into = c("Label", "Pyr"), sep = "_") %>%
  mutate(Condition = "Drought")

df.acetone <- rbind(df.p.f, df.d.f)

# calculate "pre" and "water" labeling 13C flux, so that we can calcualte 13C flux from pyruvate
data_sub_pre <- df.acetone %>%
  filter(Time <= 0)

NP_pre <- mean(data_sub_pre$Flux)

data_sub_water <- df.acetone %>%
  filter(Pyr == "W")

NP_water <- mean(data_sub_water$Flux)

NP <- mean(NP_pre, NP_water)

# filter time points to between 0 and 48 hrs
df.acetone.f <- df.acetone %>%
  mutate(Time.h = Time * 48) %>%
  select(-Time) %>%
  filter(Time.h < 48 &
           Time.h > 0)

df.acetone.f$Time.h
df.acetone.f$Pyr
df.acetone.f$Condition
df.acetone.f$Label
df.acetone.f$Flux


###calculate area under the curve for C1 and C2 during pre-drought and drought
#Calculate AUC separately for each collar
df.acetone.f.d <- df.acetone.f %>%
  filter(Condition == "Drought")
df.acetone.f.p <- df.acetone.f %>%
  filter(Condition == "Pre-Drought") 



#Define samples that are either C1 or C2 labeled
Label.C1a = c("S1P1", "S1P3", "S1P5", "S2P1", "S2P3", "S2P7", "S3P1", "S3P3", "S3P5")
Label.C1b = c("S1P1", "S1P3", "S1P5", "S2P1", "S2P3", "S2P5", "S3P1", "S3P3", "S3P5")
Label.C2a = c("S1P2", "S1P4", "S1P6", "S2P2", "S2P4", "S2P6", "S3P2", "S3P4", "S3P6")
Label.C2b = c("S1P2", "S1P4", "S1P6", "S2P4", "S2P6", "S2P8", "S3P2", "S3P4", "S3P6")
Label.W = c("S1P2", "S1P4")

#create empty lists for each category
AUC_list_PD_C1 <- list()
AUC_list_PD_C2 <- list()
AUC_list_D_C1 <- list()
AUC_list_D_C2 <- list()
AUC_list_D_W <- list()


for (i in Label.C1a) {
  x <- df.acetone.f.p %>%
    filter(Label == i &
             Pyr != "W")
  i.AUC <- auc(x$Time.h, x$Flux) 
  AUC_list_PD_C1[[i]] <- i.AUC
}


for (i in Label.C2b) {
  x <- df.acetone.f.p %>%
    filter(Label == i &
      Pyr != "W")
  i.AUC <- auc(x$Time.h, x$Flux) 
  AUC_list_PD_C2[[i]] <- i.AUC
}


for (i in Label.C1b) {
  x <- df.acetone.f.d %>%
    filter(Label == i &
      Pyr != "W")
  i.AUC <- auc(x$Time.h, x$Flux) 
  AUC_list_D_C1[[i]] <- i.AUC
}


for (i in Label.C2a) {
  x <- df.acetone.f.d %>%
    filter(Label == i &
             Pyr != "W")
  i.AUC <- auc(x$Time.h, x$Flux) 
  AUC_list_D_C2[[i]] <- i.AUC
}

for (i in Label.W) {
  x <- df.acetone.f.d %>%
    filter(Label == i &
             Pyr == "W")
  i.AUC <- auc(x$Time.h, x$Flux) 
  AUC_list_D_W[[i]] <- i.AUC
}

#Create a data frame with area under the curves for each sample for statistical analysis and plotting
#start by converting lists to data frames, and renaming columns in order to merge

AUC.D.C1 <- as.data.frame(do.call(rbind,AUC_list_D_C1))
AUC.D.C1$Label <- rownames(AUC.D.C1)
AUC.D.C1$AUC <- AUC.D.C1$V1
AUC.D.C1$Pyr <- "C1"
AUC.D.C1$Cond <- "Drought"
AUC.D.C1$V1 <- NULL

AUC.D.C2 <- as.data.frame(do.call(rbind,AUC_list_D_C2))
AUC.D.C2$Label <- rownames(AUC.D.C2)
AUC.D.C2$AUC <- AUC.D.C2$V1
AUC.D.C2$Pyr <- "C2"
AUC.D.C2$Cond <- "Drought"
AUC.D.C2$V1 <- NULL

AUC.PD.C1 <- as.data.frame(do.call(rbind, AUC_list_PD_C1))
AUC.PD.C1$Label <- rownames(AUC.PD.C1)
AUC.PD.C1$AUC <- AUC.PD.C1$V1
AUC.PD.C1$Pyr <- "C1"
AUC.PD.C1$Cond <- "Pre_Drought"
AUC.PD.C1$V1 <- NULL

AUC.PD.C2 <- as.data.frame(do.call(rbind, AUC_list_PD_C2))
AUC.PD.C2$Label <- rownames(AUC.PD.C2)
AUC.PD.C2$AUC <- AUC.PD.C2$V1
AUC.PD.C2$Pyr <- "C2"
AUC.PD.C2$Cond <- "Pre_Drought"
AUC.PD.C2$V1 <- NULL

AUC.D.W <- as.data.frame(do.call(rbind, AUC_list_D_W))
AUC.D.W$Label <- rownames(AUC.D.W)
AUC.D.W$AUC <- AUC.D.W$V1
AUC.D.W$Pyr <- "W"
AUC.D.W$Cond <- "Drought"
AUC.D.W$V1 <- NULL

AUC.W.avg <- mean(AUC.D.W$AUC)

# merge dataframes, AUC is in umol/m2
AUC.C1 <- rbind(AUC.D.C1, AUC.PD.C1)
AUC.C2 <- rbind(AUC.D.C2, AUC.PD.C2)
AUC.0 <- rbind(AUC.C1, AUC.C2)
AUC.W <- rbind(AUC.C1, AUC.C2, AUC.D.W)


# calculate total mol of C (note, AUC is in mmol/m2, diameter of collar = 20 cm which means area is 314.16 cm2) 
AUC <- AUC.0 %>%
  mutate(AUCa = AUC - AUC.W.avg) %>%
  mutate(umol = ((AUCa * 314.16) / 10000) * 1000) %>%
  mutate(mol = umol / 1000000) %>%
  mutate(mass.13C.g = mol * 13) #add column with total mass of 13C (molar mass of 13C is 13 g/mol)

# calculate % of 13C pyruvate in 13C-acetone (note, 0.0011228385 mol of 13C-pyruvate added,
# 0.0034067682 mol of C, but only one carbon atom labeled, so 0.0011344894 mol of labeled C)
mass.13C.pyr.g = 0.0011228385 * 13

AUC <- AUC %>%
  mutate(perc.mol = (mol/0.0011228385) *100) %>%
  mutate(perc.mass = ((mass.13C.g / mass.13C.pyr.g) * 100))

# calculate average aucs and nmol of C
AUC_summary <- AUC %>%
  group_by(Cond, Pyr) %>%
  summarise_at(vars(c(AUC, umol, perc.mass)), list(name = mean))


# plot cummalative flux of 13C pyruvate
plot <- AUC %>%
  mutate(Cond=factor(Cond, levels = c("Pre_Drought", "Drought"))) %>%
  ggplot(aes(x= Cond, y = umol, fill = Cond)) +
  geom_boxplot() +
  facet_grid(.~Pyr) +
  scale_fill_manual(values = colors,
                    breaks = c("Pre_Drought", "Drought"),
                    labels = c("Pre Drought", "Drought")) +
  labs(y = expression("Cummulative efflux of "^13*"C-acetone (umol)")) +
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        legend.title = element_blank())
plot
ggsave("./Figures/CO2-VOCs/total-C-acetone-flux.png",units=c('in'),width=6, height = 4,dpi=300,plot)

# plot cummalative percentage of of 13C pyruvate
plot <- AUC %>%
  mutate(Cond=factor(Cond, levels = c("Pre_Drought", "Drought"))) %>%
  ggplot(aes(x= Cond, y = perc.mass, fill = Cond)) +
  geom_boxplot() +
  facet_grid(.~Pyr) +
  scale_fill_manual(values = colors,
                    breaks = c("Pre_Drought", "Drought"),
                    labels = c("Pre Drought", "Drought")) +
  labs(y = expression("Cummulative efflux of "^13*"C-acetone (% of 13C from pyruvate)")) +
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        legend.title = element_blank())
plot
ggsave("./Figures/CO2-VOCs/C-acetone-perc-flux.png",units=c('in'),width=6, height = 4,dpi=300,plot)

# plot cummalative flux of 13C pyruvate for just C2, because only C2 showed enrichment
plot <- AUC %>%
  mutate(Cond=factor(Cond, levels = c("Pre_Drought", "Drought"))) %>%
  filter(Pyr == "C2") %>%
  ggplot(aes(x= Cond, y = umol, fill = Cond)) +
  geom_boxplot() + 
  scale_fill_manual(values = colors,
                    breaks = c("Pre_Drought", "Drought"),
                    labels = c("Pre Drought", "Drought")) +
  labs(y = expression(""^13*"C-acetone "*mu*"mol")) +
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        legend.title = element_blank())
plot
ggsave("./Figures/CO2-VOCs/total-C-acetone-flux-C2.png",units=c('in'),width=4, height = 4,dpi=300,plot)

# plot cummalative percentage of of 13C pyruvate
plot <- AUC %>%
  mutate(Cond=factor(Cond, levels = c("Pre_Drought", "Drought"))) %>%
  filter(Pyr == "C2") %>%
  ggplot(aes(x= Cond, y = perc.mass, fill = Cond)) +
  geom_boxplot() +
  scale_fill_manual(values = colors,
                    breaks = c("Pre_Drought", "Drought"),
                    labels = c("Pre Drought", "Drought")) +
  labs(y = expression("Cummulative efflux of "^13*"C-acetone (% of 13C from pyruvate)")) +
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        legend.title = element_blank())
plot
ggsave("./Figures/CO2-VOCs/C-acetone-perc-flux-C2.png",units=c('in'),width=4, height = 4,dpi=300,plot)



