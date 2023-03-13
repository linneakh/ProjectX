
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

ratios_per_set <- data_sub_13C_flux %>%  #data_sub_13C_flux object created in Fig2-B-G-H-I-J-cumulative_fluxes.RMD
  mutate(label_state = ifelse(Trel >0, "post", "pre")) %>% 
  filter(label_state == "post") %>%
  filter(Pyruv != "W") %>%
  dplyr::select(-label_state) %>%
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
  dplyr::select(-c(Trel, D13C, Flux, AtF, AtFe, Flux.h, Label)) %>%
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
  labs(x = "Time post pyruvate (h)", y= expression("C allocation to biosynthesis ("*mu*"mol m"^-2*"h"^-1*")")) +
  #facet_wrap(. ~ Site) + 
  scale_fill_manual(values = colors,
                    breaks = c("pre_drought", "drought"),
                    labels = c("Pre Drought", "Drought")) + 
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        legend.title = element_blank())
ggsave("Figures/Fig2-S1-CO2-VOCs/Fig2C-C1-C2-difference-13C-CO2-py-w-set4.png", width = 4, height = 4.1, dpi = 1000)

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
#ggsave("Figures/Fig2-S1-CO2-VOCs/FigS1-C1-C2-ratio-sum-13C-CO2-py-w-set4.png", width = 5, height = 4, dpi = 1000)

####statistical analysis
# ratios for differences in C1 - C2 deltaC-CO2r - new calculations for 13C-CO2 from pyruvate
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

#linear mixed effect
# subset to time = 3. (DF = 6, t-value = -2.40, p = 0.0533)
statdat.3 <- ratios_per_set %>%
  filter(Time == "3") 
lme.3 <- lme(C1_C2_diff_fCO2_13_py ~ Condition,
                  random = list(Site = ~1,
                                Set =  ~1),
                  data = statdat.3,
                  weights =  varIdent(form = ~1|Condition),
             na.action=na.omit
)
summary(lme.3)

# subset to time = 6. (DF=6, t-value=-3.52, p=0.0124)
statdat.6 <- ratios_per_set %>%
  filter(Time == "6") 
lme.6 <- lme(C1_C2_diff_fCO2_13_py ~ Condition,
             random = list(Site = ~1,
                           Set =  ~1),
             data = statdat.6,
             weights =  varIdent(form = ~1|Condition),
             na.action=na.omit
)
summary(lme.6)

# subset to time = 12. (t-value = -4.26, p = 0.0053)
statdat.12 <- ratios_per_set %>%
  filter(Time == "12") 
lme.12 <- lme(C1_C2_diff_fCO2_13_py ~ Condition,
             random = list(Site = ~1,
                           Set =  ~1),
             data = statdat.12,
             weights =  varIdent(form = ~1|Condition),
             na.action=na.omit
)
summary(lme.12)

# subset to time = 18. (t=-4.90, p = 0.0027)
statdat.18 <- ratios_per_set %>%
  filter(Time == "18") 
lme.18 <- lme(C1_C2_diff_fCO2_13_py ~ Condition,
             random = list(Site = ~1,
                           Set =  ~1),
             data = statdat.18,
             weights =  varIdent(form = ~1|Condition),
             na.action=na.omit
)
summary(lme.18)

# subset to time = 24. (t=-2.4, p=0.0533)
statdat.24 <- ratios_per_set %>%
  filter(Time == "24") 
lme.24 <- lme(C1_C2_diff_fCO2_13_py ~ Condition,
             random = list(Site = ~1,
                           Set =  ~1),
             data = statdat.3,
             weights =  varIdent(form = ~1|Condition),
             na.action=na.omit
)
summary(lme.24)

# subset to time = 30. (t=00.88, p = 0.4116)
statdat.30 <- ratios_per_set %>%
  filter(Time == "30") 
lme.30 <- lme(C1_C2_diff_fCO2_13_py ~ Condition,
             random = list(Site = ~1,
                           Set =  ~1),
             data = statdat.30,
             weights =  varIdent(form = ~1|Condition),
             na.action=na.omit
)
summary(lme.30)

# subset to time = 36. (t=-0.64, p=0.5506)
statdat.36 <- ratios_per_set %>%
  filter(Time == "36") 
lme.36 <- lme(C1_C2_diff_fCO2_13_py ~ Condition,
             random = list(Site = ~1,
                           Set =  ~1),
             data = statdat.36,
             weights =  varIdent(form = ~1|Condition),
             na.action=na.omit
)
summary(lme.36)

# subset to time = 42. (t=-0.15, p = 0.885)
statdat.42 <- ratios_per_set %>%
  filter(Time == "42") 
lme.42<- lme(C1_C2_diff_fCO2_13_py ~ Condition,
              random = list(Site = ~1,
                            Set =  ~1),
              data = statdat.42,
              weights =  varIdent(form = ~1|Condition),
              na.action=na.omit
)
summary(lme.42)

# subset to time = 48. (t=0.400, p=0.7056)
statdat.48 <- ratios_per_set %>%
  filter(Time == "48") 
lme.48 <- lme(C1_C2_diff_fCO2_13_py ~ Condition,
              random = list(Site = ~1,
                            Set =  ~1),
              data = statdat.48,
              weights =  varIdent(form = ~1|Condition),
              na.action=na.omit
)
summary(lme.48)

# C1/(1 + C2) ratios of 13C-CO2-from-pyruvate- statistical tests

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

#linear mixed effect
# subset to time = 3. (0.57789  0.5844)
statdat.3 <- ratios_per_set %>%
  filter(Time == "3") 
lme.3 <- lme(C1_C2_ratio_sum_fCO2_13_py ~ Condition,
             random = list(Site = ~1,
                           Set =  ~1),
             data = statdat.3,
             weights =  varIdent(form = ~1|Condition),
             na.action=na.omit
)
summary(lme.3)

# subset to time = 6. (-0.89921  0.4032)
statdat.6 <- ratios_per_set %>%
  filter(Time == "6") 
lme.6 <- lme(C1_C2_ratio_sum_fCO2_13_py ~ Condition,
             random = list(Site = ~1,
                           Set =  ~1),
             data = statdat.6,
             weights =  varIdent(form = ~1|Condition),
             na.action=na.omit
)
summary(lme.6)

# subset to time = 12. (-1.47846  0.1898)
statdat.12 <- ratios_per_set %>%
  filter(Time == "12") 
lme.12 <- lme(C1_C2_ratio_sum_fCO2_13_py ~ Condition,
              random = list(Site = ~1,
                            Set =  ~1),
              data = statdat.12,
              weights =  varIdent(form = ~1|Condition),
              na.action=na.omit
)
summary(lme.12)

# subset to time = 18. (-1.35718  0.2236)
statdat.18 <- ratios_per_set %>%
  filter(Time == "18") 
lme.18 <- lme(C1_C2_ratio_sum_fCO2_13_py ~ Condition,
              random = list(Site = ~1,
                            Set =  ~1),
              data = statdat.18,
              weights =  varIdent(form = ~1|Condition),
              na.action=na.omit
)
summary(lme.18)

# subset to time = 24. (0.57789  0.5844)
statdat.24 <- ratios_per_set %>%
  filter(Time == "24") 
lme.24 <- lme(C1_C2_ratio_sum_fCO2_13_py ~ Condition,
              random = list(Site = ~1,
                            Set =  ~1),
              data = statdat.3,
              weights =  varIdent(form = ~1|Condition),
              na.action=na.omit
)
summary(lme.24)

# subset to time = 30. (-0.623992  0.5556)
statdat.30 <- ratios_per_set %>%
  filter(Time == "30") 
lme.30 <- lme(C1_C2_ratio_sum_fCO2_13_py ~ Condition,
              random = list(Site = ~1,
                            Set =  ~1),
              data = statdat.30,
              weights =  varIdent(form = ~1|Condition),
              na.action=na.omit
)
summary(lme.30)

# subset to time = 36. (-0.203234   0.847)
statdat.36 <- ratios_per_set %>%
  filter(Time == "36") 
lme.36 <- lme(C1_C2_ratio_sum_fCO2_13_py ~ Condition,
              random = list(Site = ~1,
                            Set =  ~1),
              data = statdat.36,
              weights =  varIdent(form = ~1|Condition),
              na.action=na.omit
)
summary(lme.36)

# subset to time = 42. (0.503519  0.6360)
statdat.42 <- ratios_per_set %>%
  filter(Time == "42") 
lme.42<- lme(C1_C2_ratio_sum_fCO2_13_py ~ Condition,
             random = list(Site = ~1,
                           Set =  ~1),
             data = statdat.42,
             weights =  varIdent(form = ~1|Condition),
             na.action=na.omit
)
summary(lme.42)

# subset to time = 48. (1.147460  0.3031)
statdat.48 <- ratios_per_set %>%
  filter(Time == "48") 
lme.48 <- lme(C1_C2_ratio_sum_fCO2_13_py ~ Condition,
              random = list(Site = ~1,
                            Set =  ~1),
              data = statdat.48,
              weights =  varIdent(form = ~1|Condition),
              na.action=na.omit
)
summary(lme.48)

# C1/(1 + C2) ratios of 13C-CO2-from-pyruvate- statistical tests by site

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

