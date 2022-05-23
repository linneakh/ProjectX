
# Visualize Pyruvate Co2 Fluxes

# Johannes Ingrisch, modified by Linnea Honeker
# 2021-08-26

# load packages
library(tidyverse)
library(plotly)
library(ggthemes)
library(patchwork)
library(flux)

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
  select(Label, Condition, Pyruv, Trel, Flux, D13C, Qc_iso, Qc_flux)

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

# separate out data from before and after, calculate AtF natural abundance for pre
vpdb = 0.01123720

data_sub_pre <- data_sub %>%
  filter(Trel <= 0) %>%
  mutate(R = (D13C*vpdb)/1000 + vpdb) %>%
  mutate(AtFna = R / (1 + R))

AtFna_mean = mean(data_sub_pre$AtFna)

data_sub_post <- data_sub %>%
  filter(Trel > 0)

# Calculate AtF, flux13, flux12
data_sub_13C_flux <- data_sub_post %>%
  mutate(R = (D13C*vpdb)/1000 + vpdb) %>%
  mutate(AtF = R / (1 + R)) %>%
  mutate(fCO2_13 = AtF * Flux) %>%
  mutate(AtFe = AtF - AtFna_mean) %>%
  mutate(fCO2_13_py = AtFe * Flux)

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
