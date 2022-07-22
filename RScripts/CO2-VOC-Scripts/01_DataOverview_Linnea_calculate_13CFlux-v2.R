
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



# calculate C1/C2
AUC_list_PD_C1_avg/AUC_list_PD_C2_avg # 2.87
AUC_list_D_C1_avg/AUC_list_D_C2_avg # 2.56


