
# Visualize Pyruvate Co2 Fluxes

# Johannes Ingrisch, modified by Linnea Honeker
# 2021-08-26

# load packages
library(tidyverse)
library(plotly)
library(ggthemes)
library(patchwork)
library(zoo)
library(dplyr)
library(dun.test)
library(DescTools)

theme_custom <- theme_few()

theme_set(theme_custom)

colors = c("#2ECC71" ,"#F39C12")
colors_pyruv = c("blue" ,"pink")


###C1/C2 ratios, separate by site and divide C1:C2 ratios of nearby locations (paired locations)
# modify data_sub table to add time intervals and calculate C1 and C2 ratios per site 
ratios_per_set <- data_sub_13C_flux %>%  #data_sun_13C_flux object created in Fig2-B-G-H-I-J-cumulative_fluxes.RMD
  mutate(label_state = ifelse(Trel >0, "post", "pre")) %>% 
  filter(label_state == "post") %>%
  filter(Pyruv != "W") %>%
  select(-label_state) %>%
  filter(Label != "S2.P2") %>%
  mutate(Site = case_when(
    startsWith(Label, "S1") ~ "Site1",
    startsWith(Label, "S2") ~ "Site2",
    startsWith(Label, "S3") ~ "Site3"
  )) %>%
  mutate(Set = case_when(
    endsWith(Label, "P1") ~ "Set1",
    endsWith(Label, "P2") ~ "Set1",
    endsWith(Label, "P3") ~ "Set2",
    endsWith(Label, "P4") ~ "Set2",
    endsWith(Label, "P5") ~ "Set3",
    endsWith(Label, "P6") ~ "Set3",
    endsWith(Label, "P7") ~ "Set4",
    endsWith(Label, "P8") ~ "Set4"
  )) %>%
  mutate(Time = case_when(
    (Trel <= 3) ~ "3",
    (Trel > 3 & Trel <= 6) ~ "6",
    (Trel > 6 & Trel <= 12) ~ "12",
    (Trel > 12 & Trel <= 18) ~ "18",
    (Trel > 18 & Trel <= 24) ~ "24",
    (Trel > 24 & Trel <= 30) ~ "30",
    (Trel > 30 & Trel <= 36) ~ "36",
    (Trel > 36 & Trel <= 42) ~ "42",
    (Trel > 42 & Trel <= 100) ~ "48"
  )) %>%
  mutate(Time=factor(Time, levels = c("3", "6", "12", "18", "24", "30", "36", "42", "48"))) %>%
  select(-c(Trel, D13C, Flux, AtF, AtFe, Flux.h, Label)) %>%
  group_by(Condition, Site, Set, Pyruv, Time) %>%
  summarise_all(funs(mean)) %>%
  pivot_wider(names_from = Pyruv, values_from = c(fCO2_13, fCO2_13_py)) %>%
  mutate(C1_C2_diff_fCO2_13_py = fCO2_13_py_C1 - fCO2_13_py_C2) %>%
  mutate(C1_C2_ratio_sum_fCO2_13_py = fCO2_13_py_C1/(fCO2_13_py_C1 + fCO2_13_py_C2))

#set Condition factor in order
ratios_per_set$Condition <- factor(ratios_per_set$Condition, c("pre_drought", "drought"))

# graph of all differences (flux13Conly)
ratios_per_set %>%
  ggplot(aes(x = Time, y = C1_C2_diff_fCO2_13_py, fill = Condition)) +
  geom_boxplot() +
  labs(x = "Time post pyruvate (h)", y= expression("C allocation to biosynthesis ("*mu*"mol m"^-2*"s"^-1*")")) +
  #facet_wrap(. ~ Site) + 
  scale_fill_manual(values = colors,
                    breaks = c("pre_drought", "drought"),
                    labels = c("Pre Drought", "Drought")) + 
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        legend.title = element_blank())
ggsave("Figures/Fig2-S1-CO2-VOCs/Fig2C-C1-C2-difference-13C-CO2-py-w-set4.png", width = 4, height = 4, dpi = 1000)

# graph of all ratios (flux 13C from pyruvate, C1/(C1+C2))
ratios_per_set %>%
  #filter(Set != "Set4") %>%
  ggplot(aes(x = Time, y = C1_C2_ratio_sum_fCO2_13_py, fill = Condition)) +
  geom_boxplot() +
  labs(x = "Time post pyruvate (h)", y= expression("C allocation to biosynthesis (proportion)")) +
  #facet_wrap(~Site) + 
  scale_fill_manual(values = colors,
                    breaks = c("pre_drought", "drought"),
                    labels = c("Pre Drought", "Drought")) + 
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        legend.title = element_blank())
ggsave("Figures/Fig2-S1-CO2-VOCs/FigS1-C1-C2-ratio-sum-13C-CO2-py-w-set4.png", width = 5, height = 4, dpi = 1000)

####statistical analysis
# ratios for differences in C1 - C2 deltaC-CO2r - new calculations for 13C-CO2 from pyruvate
#t-test
ratios_per_set_3hr<- ratios_per_set %>%
  filter(Time == "3" )
t.test(C1_C2_diff_fCO2_13_py ~ Condition, data = ratios_per_set_3hr )

ratios_per_set_6hr<- ratios_per_set %>%
  filter(Time == "6" )
t.test(C1_C2_diff_fCO2_13_py ~ Condition, data = ratios_per_set_6hr ) #p = 0.004931

ratios_per_set_12hr<- ratios_per_set %>%
  filter(Time == "12" )
t.test(C1_C2_diff_fCO2_13_py ~ Condition, data = ratios_per_set_12hr ) #p=0.001838

ratios_per_set_18hr<- ratios_per_set %>%
  filter(Time == "18" )
t.test(C1_C2_diff_fCO2_13_py ~ Condition, data = ratios_per_set_18hr ) #p=0.0004503

ratios_per_set_24hr<- ratios_per_set %>%
  filter(Time == "24" )
t.test(C1_C2_diff_fCO2_13_py ~ Condition, data = ratios_per_set_24hr ) #0.03411

ratios_per_set_30hr<- ratios_per_set %>%
  filter(Time == "30" )
t.test(C1_C2_diff_fCO2_13_py ~ Condition, data = ratios_per_set_30hr )

ratios_per_set_36hr<- ratios_per_set %>%
  filter(Time == "36" )
t.test(C1_C2_diff_fCO2_13_py ~ Condition, data = ratios_per_set_36hr )

ratios_per_set_42hr<- ratios_per_set %>%
  filter(Time == "42" )
t.test(C1_C2_diff_fCO2_13_py ~ Condition, data = ratios_per_set_42hr )

ratios_per_set_48hr<- ratios_per_set %>%
  filter(Time == "48" )
t.test(C1_C2_diff_fCO2_13_py ~ Condition, data = ratios_per_set_48hr )

#wilcoxan test
ratios_per_set_3hr<- ratios_per_set %>%
  filter(Time == "3" )
wilcox.test(C1_C2_diff_fCO2_13_py ~ Condition, data = ratios_per_set_3hr ) # 0.03595

ratios_per_set_6hr<- ratios_per_set %>%
  filter(Time == "6" )
wilcox.test(C1_C2_diff_fCO2_13_py ~ Condition, data = ratios_per_set_6hr ) # 0.007898

ratios_per_set_12hr<- ratios_per_set %>%
  filter(Time == "12" )
wilcox.test(C1_C2_diff_fCO2_13_py ~ Condition, data = ratios_per_set_12hr ) #p=0.001563

ratios_per_set_18hr<- ratios_per_set %>%
  filter(Time == "18" )
wilcox.test(C1_C2_diff_fCO2_13_py ~ Condition, data = ratios_per_set_18hr ) #p=0.002468

ratios_per_set_24hr<- ratios_per_set %>%
  filter(Time == "24" )
wilcox.test(C1_C2_diff_fCO2_13_py ~ Condition, data = ratios_per_set_24hr ) #0.02065

ratios_per_set_30hr<- ratios_per_set %>%
  filter(Time == "30" )
wilcox.test(C1_C2_diff_fCO2_13_py ~ Condition, data = ratios_per_set_30hr )

ratios_per_set_36hr<- ratios_per_set %>%
  filter(Time == "36" )
wilcox.test(C1_C2_diff_fCO2_13_py ~ Condition, data = ratios_per_set_36hr )

ratios_per_set_42hr<- ratios_per_set %>%
  filter(Time == "42" )
wilcox.test(C1_C2_diff_fCO2_13_py ~ Condition, data = ratios_per_set_42hr )

ratios_per_set_48hr<- ratios_per_set %>%
  filter(Time == "48" )
wilcox.test(C1_C2_diff_fCO2_13_py ~ Condition, data = ratios_per_set_48hr )

# C1/(1 + C2) ratios of 13C-CO2-from-pyruvate- statistical tests
#t-test
ratios_per_set_3hr<- ratios_per_set %>%
  filter(Time == "3" )
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_3hr )

ratios_per_set_6hr<- ratios_per_set %>%
  filter(Time == "6" )
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_6hr ) #p = 0.004931

ratios_per_set_12hr<- ratios_per_set %>%
  filter(Time == "12" )
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_12hr ) #p=0.001838

ratios_per_set_18hr<- ratios_per_set %>%
  filter(Time == "18" )
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_18hr ) #p=0.0004503

ratios_per_set_24hr<- ratios_per_set %>%
  filter(Time == "24" )
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_24hr ) #0.03411

ratios_per_set_30hr<- ratios_per_set %>%
  filter(Time == "30" )
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_30hr )

ratios_per_set_36hr<- ratios_per_set %>%
  filter(Time == "36" )
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_36hr )

ratios_per_set_42hr<- ratios_per_set %>%
  filter(Time == "42" )
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_42hr )

ratios_per_set_48hr<- ratios_per_set %>%
  filter(Time == "48" )
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_48hr )

#wilcoxan test
ratios_per_set_3hr<- ratios_per_set %>%
  filter(Time == "3" )
wilcox.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_3hr ) # 0.03595

ratios_per_set_6hr<- ratios_per_set %>%
  filter(Time == "6" )
wilcox.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_6hr ) # 0.007898

ratios_per_set_12hr<- ratios_per_set %>%
  filter(Time == "12" )
wilcox.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_12hr ) #p=0.001563

ratios_per_set_18hr<- ratios_per_set %>%
  filter(Time == "18" )
wilcox.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_18hr ) #p=0.002468

ratios_per_set_24hr<- ratios_per_set %>%
  filter(Time == "24" )
wilcox.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_24hr ) #0.02065

ratios_per_set_30hr<- ratios_per_set %>%
  filter(Time == "30" )
wilcox.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_30hr )

ratios_per_set_36hr<- ratios_per_set %>%
  filter(Time == "36" )
wilcox.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_36hr )

ratios_per_set_42hr<- ratios_per_set %>%
  filter(Time == "42" )
wilcox.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_42hr )

ratios_per_set_48hr<- ratios_per_set %>%
  filter(Time == "48" )
wilcox.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_48hr )

# C1/(1 + C2) ratios of 13C-CO2-from-pyruvate- statistical tests by site
#t-test Site 1
ratios_per_set_3hr_Site1<- ratios_per_set %>%
  filter(Time == "3" &
           Site == "Site1")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_3hr_Site1)

ratios_per_set_6hr_Site1<- ratios_per_set %>%
  filter(Time == "6" &
           Site == "Site1")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_6hr_Site1)

ratios_per_set_12hr_Site1<- ratios_per_set %>%
  filter(Time == "12" &
           Site == "Site1")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_12hr_Site1)

ratios_per_set_18hr_Site1<- ratios_per_set %>%
  filter(Time == "18" &
           Site == "Site1")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_18hr_Site1)

ratios_per_set_24hr_Site1<- ratios_per_set %>%
  filter(Time == "24" &
           Site == "Site1")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_24hr_Site1)

ratios_per_set_30hr_Site1<- ratios_per_set %>%
  filter(Time == "30" &
           Site == "Site1")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_30hr_Site1)

ratios_per_set_36hr_Site1<- ratios_per_set %>%
  filter(Time == "36" &
           Site == "Site1")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_36hr_Site1)

ratios_per_set_42hr_Site1<- ratios_per_set %>%
  filter(Time == "42" &
           Site == "Site1")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_42hr_Site1)

ratios_per_set_48hr_Site1<- ratios_per_set %>%
  filter(Time == "48" &
           Site == "Site1")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_48hr_Site1)

ratios_per_set_3hr_Site3<- ratios_per_set %>%
  filter(Time == "3" &
           Site == "Site3")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_3hr_Site3)

ratios_per_set_6hr_Site3<- ratios_per_set %>%
  filter(Time == "6" &
           Site == "Site3")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_6hr_Site3)

ratios_per_set_12hr_Site3<- ratios_per_set %>%
  filter(Time == "12" &
           Site == "Site3")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_12hr_Site3)

ratios_per_set_18hr_Site3<- ratios_per_set %>%
  filter(Time == "18" &
           Site == "Site3")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_18hr_Site3)

ratios_per_set_24hr_Site3<- ratios_per_set %>%
  filter(Time == "24" &
           Site == "Site3")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_24hr_Site3)

ratios_per_set_30hr_Site3<- ratios_per_set %>%
  filter(Time == "30" &
           Site == "Site3")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Condition, data = ratios_per_set_30hr_Site3)

## C1/C1+C2 are sites different
ratios_per_set_3hr_D<- ratios_per_set %>%
  filter(Time == "3" &
           Condition == "drought" &
           Site != "Site2")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Site, data = ratios_per_set_3hr_D)

ratios_per_set_6hr_D<- ratios_per_set %>%
  filter(Time == "6" &
           Condition == "drought" &
           Site != "Site2")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Site, data = ratios_per_set_6hr_D)

ratios_per_set_12hr_D<- ratios_per_set %>%
  filter(Time == "12" &
           Condition == "drought" &
           Site != "Site2")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Site, data = ratios_per_set_12hr_D)

ratios_per_set_18hr_D<- ratios_per_set %>%
  filter(Time == "18" &
           Condition == "drought" &
           Site != "Site2")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Site, data = ratios_per_set_18hr_D)

ratios_per_set_24hr_D<- ratios_per_set %>%
  filter(Time == "24" &
           Condition == "drought" &
           Site != "Site2")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Site, data = ratios_per_set_24hr_D)

ratios_per_set_30hr_D<- ratios_per_set %>%
  filter(Time == "30" &
           Condition == "drought" &
           Site != "Site2")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Site, data = ratios_per_set_30hr_D)

ratios_per_set_36hr_D<- ratios_per_set %>%
  filter(Time == "36" &
           Condition == "drought" &
           Site != "Site2")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Site, data = ratios_per_set_36hr_D)

ratios_per_set_42hr_D<- ratios_per_set %>%
  filter(Time == "42" &
           Condition == "drought" &
           Site != "Site2")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Site, data = ratios_per_set_42hr_D)

ratios_per_set_48hr_D<- ratios_per_set %>%
  filter(Time == "48" &
           Condition == "drought" &
           Site != "Site2")
t.test(C1_C2_ratio_sum_fCO2_13_py ~ Site, data = ratios_per_set_48hr_D)

## C1/C1+C2 are sites different wilcoxan test
ratios_per_set_3hr_D<- ratios_per_set %>%
  filter(Time == "3" &
           Condition == "drought" &
           Site != "Site2")
wilcox.test(C1_C2_ratio_sum_fCO2_13_py ~ Site, data = ratios_per_set_3hr_D)

ratios_per_set_6hr_D<- ratios_per_set %>%
  filter(Time == "6" &
           Condition == "drought" &
           Site != "Site2")
wilcox.test(C1_C2_ratio_sum_fCO2_13_py ~ Site, data = ratios_per_set_6hr_D)

ratios_per_set_12hr_D<- ratios_per_set %>%
  filter(Time == "12" &
           Condition == "drought" &
           Site != "Site2")
wilcox.test(C1_C2_ratio_sum_fCO2_13_py ~ Site, data = ratios_per_set_12hr_D)

ratios_per_set_18hr_D<- ratios_per_set %>%
  filter(Time == "18" &
           Condition == "drought" &
           Site != "Site2")
wilcox.test(C1_C2_ratio_sum_fCO2_13_py ~ Site, data = ratios_per_set_18hr_D)

ratios_per_set_24hr_D<- ratios_per_set %>%
  filter(Time == "24" &
           Condition == "drought" &
           Site != "Site2")
wilcox.test(C1_C2_ratio_sum_fCO2_13_py ~ Site, data = ratios_per_set_24hr_D)

ratios_per_set_30hr_D<- ratios_per_set %>%
  filter(Time == "30" &
           Condition == "drought" &
           Site != "Site2")
wilcox.test(C1_C2_ratio_sum_fCO2_13_py ~ Site, data = ratios_per_set_30hr_D)

ratios_per_set_36hr_D<- ratios_per_set %>%
  filter(Time == "36" &
           Condition == "drought" &
           Site != "Site2")
wilcox.test(C1_C2_ratio_sum_fCO2_13_py ~ Site, data = ratios_per_set_36hr_D)

ratios_per_set_42hr_D<- ratios_per_set %>%
  filter(Time == "42" &
           Condition == "drought" &
           Site != "Site2")
wilcox.test(C1_C2_ratio_sum_fCO2_13_py ~ Site, data = ratios_per_set_42hr_D)

ratios_per_set_48hr_D<- ratios_per_set %>%
  filter(Time == "48" &
           Condition == "drought" &
           Site != "Site2")
wilcox.test(C1_C2_ratio_sum_fCO2_13_py ~ Site, data = ratios_per_set_48hr_D)

