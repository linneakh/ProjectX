
# Visualize Pyruvate Co2 Fluxes

# Johannes Ingrisch, modified by Linnea Honeker
# 2021-08-26

# load packages
library(tidyverse)
library(plotly)
library(ggthemes)
library(patchwork)
library(flux)
library(dplyr)
library(data.table)


#set themes and colors
theme_custom <- theme_few()
theme_set(theme_custom)
colors = c("#2ECC71" ,"#F39C12")

#load data
data <- read.csv("./Data/CO2-VOCs/01_raw/S3_PyruvateLabeling_qc.csv") %>% 
  select(-X)

#add timestanp, however, this hasn't been working and not needed
data <- data %>% 
  mutate(Timestamp = lubridate::ymd_hms(Timestamp_MST, tz = "America/Phoenix"), 
         Timestamp_Label_MST = lubridate::ymd_hms(Timestamp_Label_MST, tz = "America/Phoenix")) %>% 
  select(-Timestamp_MST)

# Quickview
data <- data %>% 
  select(Label, Condition, Pyruv, Trel, Flux, D13C, AtF, Qc_iso, Qc_flux)

# quality control filtering
data_sub <- data %>% 
  filter(Qc_flux == 1, Qc_iso == 1)  %>% 
  select(-Qc_iso, -Qc_flux) %>% 
  filter(Trel > -10, Trel < 50, Flux > 0, D13C > -100)

##### calculate average AtF, below are equations from Johnny
# atF = Flux13/(flux13 + Flux12);
# d13C = (((flux13 / flux12)/Rvpdb)-1) * 1000
# F_13CO2 = Flux * AtF
# F_12CO2 = Flux - F_13CO2 = Flux * (1-AtF)
#  Rvpdb = 0.0111802 (from Werner and Brand, 2001)
# solving for AtF, Atf = d13/11.1802

# alculate AtF natural abundance for pre labeling time points
data_sub_pre <- data_sub %>%
  filter(Trel <= 0)

#calculate mean AtF natural abundance across all pre-labeling time points
AtFna_mean = mean(data_sub_pre$AtF)

# filter to all post-labeling time points
data_sub_post <- data_sub %>%
  filter(Trel > 0)

# Calculate AtF, flux13, flux12 (flux units: umol m-2 s-1)
data_sub_13C_flux <- data_sub_post %>%
  mutate(fCO2_13 = AtF * Flux) %>%
  mutate(AtFe = AtF - AtFna_mean) %>% #calculate AtF from 13C pyruvate
  mutate(fCO2_13_py = AtFe * Flux) #calculate 13C flux attributed to 13C pyruvate

####Calculations of AUCs using all data combined######
#Calculate area under the curve for drought 13C-CO2, units are umol/m2
data_sub_13C_flux_drought <- data_sub_13C_flux %>%
  filter(Condition == "drought")
auc(data_sub_13C_flux_drought$Trel, data_sub_13C_flux_drought$fCO2_13_py) #2.002319

#Calculate area under the curve for pre-drought 13C-CO2
data_sub_13C_flux_predrought <- data_sub_13C_flux %>%
  filter(Condition == "pre_drought")
auc(data_sub_13C_flux_predrought$Trel, data_sub_13C_flux_predrought$fCO2_13_py) #4.130672

#Calculate area under the curve for drought C1 13C-CO2
data_sub_13C_flux_droughtC1 <- data_sub_13C_flux %>%
  filter(Condition == "drought" &
           Pyruv == "C1")
CO2_C1_D <- auc(data_sub_13C_flux_droughtC1$Trel, data_sub_13C_flux_droughtC1$fCO2_13_py) #2.891094

#Calculate area under the curve for drought C2 13C-CO2, units are umol/m2
data_sub_13C_flux_droughtC2 <- data_sub_13C_flux %>%
  filter(Condition == "drought" &
           Pyruv == "C2")
CO2_C2_D <- auc(data_sub_13C_flux_droughtC2$Trel, data_sub_13C_flux_droughtC2$fCO2_13_py) #1.219856

#Calculate area under the curve for predrought C1 13C-CO2
data_sub_13C_flux_predroughtC1 <- data_sub_13C_flux %>%
  filter(Condition == "pre_drought" &
           Pyruv == "C1")
CO2_C1_PD <- auc(data_sub_13C_flux_predroughtC1$Trel, data_sub_13C_flux_predroughtC1$fCO2_13_py) #5.985825

#Calculate area under the curve for predrought C2 13C-CO2, units are umol/m2
data_sub_13C_flux_predroughtC2 <- data_sub_13C_flux %>%
  filter(Condition == "pre_drought" &
           Pyruv == "C2")
CO2_C2_PD <- auc(data_sub_13C_flux_predroughtC2$Trel, data_sub_13C_flux_predroughtC2$fCO2_13_py) #1.870699

# C1 - C2 diff during pre-drought (biosynthesis) =
diff_D <- CO2_C1_D - CO2_C2_D #1.67 umol/m2
diff_PD <- CO2_C1_PD - CO2_C2_PD #4.12 umol/m2

########Calculate AUC separately for each collar, then average######################
#separate out drought and pre-drought
data_sub_13C_flux_drought <- data_sub_13C_flux %>%
  filter(Condition == "drought") %>%
  mutate(Trel.s = Trel * 360) %>%
  select(-Trel)

data_sub_13C_flux_predrought <- data_sub_13C_flux %>%
  filter(Condition == "pre_drought") %>%
  mutate(Trel.s = Trel * 360) %>%
  select(-Trel)


#Define samples that are either C1 or C2 labeled
Label.C1a = c("S1.P1", "S1.P3", "S1.P5", "S2.P1", "S2.P3", "S2.P7", "S3.P1", "S3.P3", "S3.P5")
Label.C1b = c("S1.P1", "S1.P3", "S1.P5", "S2.P1", "S2.P3", "S2.P5", "S3.P1", "S3.P3", "S3.P5")
Label.C2a = c("S1.P2", "S1.P4", "S1.P6", "S2.P2", "S2.P4", "S2.P6", "S3.P2", "S3.P4", "S3.P6")
Label.C2b = c("S1.P2", "S1.P4", "S1.P6",  "S2.P4", "S2.P6","S2.P8", "S3.P2", "S3.P4", "S3.P6")

#create empty lists for each category
AUC_list_PD_C1 <- list()
AUC_list_PD_C2 <- list()
AUC_list_D_C1 <- list()
AUC_list_D_C2 <- list()

for (i in Label.C1a) {
  x <- data_sub_13C_flux_predrought %>%
  filter(Label == i)
  i.AUC <- auc(x$Trel.s, x$fCO2_13_py) 
  AUC_list_PD_C1[[i]] <- i.AUC
}
  
for (i in Label.C2b) {
  x <- data_sub_13C_flux_predrought %>%
    filter(Label == i)
  i.AUC <- auc(x$Trel.s, x$fCO2_13_py) 
  AUC_list_PD_C2[[i]] <- i.AUC
}

for (i in Label.C1b) {
  x <- data_sub_13C_flux_drought %>%
    filter(Label == i)
  i.AUC <- auc(x$Trel.s, x$fCO2_13_py) 
  AUC_list_D_C1[[i]] <- i.AUC
}

for (i in Label.C2a) {
  x <- data_sub_13C_flux_drought %>%
    filter(Label == i)
  i.AUC <- auc(x$Trel.s, x$fCO2_13_py) 
  AUC_list_D_C2[[i]] <- i.AUC
}

# Create a data frame with area under the curves for each sample for statistical analysis and plotting
# start by converting lists to data frames, and renaming columns in order to merge

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

# merge dataframes, AUC is in umol/cm2
AUC.C1 <- rbind(AUC.D.C1, AUC.PD.C1)
AUC.C2 <- rbind(AUC.D.C2, AUC.PD.C2)
AUC.0 <- rbind(AUC.C1, AUC.C2)


# calculate total mmol of C (note, AUC is in umol/m2, diameter of collar = 20 cm) 
# mmol: multiply AUC by area of collar (314.16cm2) to get umol, then divide by 1000 to convert from umol to mmol
# mol: divide mmol by 1000 to convert to mol
AUC <- AUC.0 %>%
  mutate(umol = ((AUC * 314.16) /10000)) %>%
  mutate(mol = umol / 1000000) %>%
  mutate(mass.13C.g = mol * 13) #add column with total mass of 13C (molar mass of 13C is 13 g/mol)

# calculate % of 13C pyruvate in 13C-CO2 (note, 0.10g of 13C-pyruvate added, which = 0.0011228385 mol(molar mass of 13C-pyruvate = 89.06g/mol),
# Since only one carbon atom labeled, 13C is also 0.0011344894 mol, and total mass of 13C can by calculated)
mass.13C.pyr.g = 0.0011228385 * 13

AUC <- AUC %>%
  mutate(perc.mass = (mass.13C.g/mass.13C.pyr.g) *100) %>%
  mutate(perc.mol = (mol/ 0.0011228385) * 100)

# calculate average aucs and nmol of C
AUC_summary <- AUC %>%
  group_by(Cond, Pyr) %>%
  summarise_at(vars(c(AUC, umol, mass.13C.g, perc.mass)), list(name = mean))

# total C for biosynthesis during drought
AUC_summary_wide <- AUC_summary %>%
  select(-c(AUC_name, perc.mass_name, umol_name)) %>%
  pivot_wider(names_from = c(Pyr, Cond), values_from = mass.13C.g_name) %>%
  mutate(Biosynthesis.D = C1_Drought - C2_Drought) %>%
  mutate(Biosynthesis.PD = C1_Pre_Drought - C2_Pre_Drought) %>%
  mutate(TCA_cycle.D = C1_Drought - Biosynthesis.D) %>%
  mutate(TCA_cycle.PD = C1_Pre_Drought - Biosynthesis.PD) 
 

# plot mmol of 13C pyruvate
plot <- AUC %>%
  mutate(Cond=factor(Cond, levels = c("Pre_Drought", "Drought"))) %>%
  ggplot(aes(x= Cond, y = umol, fill = Cond)) +
  geom_boxplot() +
  facet_grid(.~Pyr) +
  scale_fill_manual(values = colors,
                    breaks = c("Pre_Drought", "Drought"),
                    labels = c("Pre Drought", "Drought")) +
  labs(y = expression("Cummulative efflux of "^13*"C-CO"[2]*" (umol)")) +
  theme(text = element_text(size = 12,
                          family = "Arial",
                          color = "black"),
      legend.position = "bottom",
      axis.title.x = element_blank(),
      legend.title = element_blank())
plot
ggsave("./Figures/CO2-VOCs/total-C-CO2-flux.png",units=c('in'),width=6, height = 4,dpi=300,plot)

# plot % of 13C pyruvate
plot <- AUC %>%
  mutate(Cond=factor(Cond, levels = c("Pre_Drought", "Drought"))) %>%
  ggplot(aes(x= Cond, y = perc, fill = Cond)) +
  geom_boxplot() +
  facet_grid(.~Pyr) +
  scale_fill_manual(values = colors,
                    breaks = c("Pre_Drought", "Drought"),
                    labels = c("Pre Drought", "Drought")) +
  labs(y = expression("Cummulative efflux of "^13*"C-CO"[2]*" (% of 13C injected)")) +
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        legend.title = element_blank())
plot
ggsave("./Figures/CO2-VOCs/total-C-CO2-per-flux.png",units=c('in'),width=6, height = 4,dpi=300,plot)

# calculate C1/C2
AUC_list_PD_C1_avg/AUC_list_PD_C2_avg # 2.87
AUC_list_D_C1_avg/AUC_list_D_C2_avg # 2.56


