
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

theme_custom <- theme_few()

theme_set(theme_custom)

colors = c("#2ECC71" ,"#F39C12")
colors_pyruv = c("blue" ,"pink")


data <- read.csv("./Data/CO2-VOCs/01_raw/S3_PyruvateLabeling_qc.csv") %>% 
  select(-X)


data <- data %>% 
  mutate(Timestamp = lubridate::ymd_hms(Timestamp_MST, tz = "America/Phoenix"), 
         Timestamp_Label_MST = lubridate::ymd_hms(Timestamp_Label_MST, tz = "America/Phoenix")) %>% 
  select(-Timestamp_MST)


# Quickview


data <- data %>% 
  select(Label, Condition, Pyruv, Trel, Flux, D13C, Qc_iso, Qc_flux)

data$Condition <- factor(data$Condition, c(
  "pre_drought", "drought"))

data_sub <- data %>% 
  filter(Qc_flux == 1, Qc_iso == 1)  %>% 
  select(-Qc_iso, -Qc_flux) %>% 
  filter(Trel > -10, Trel < 50, Flux > 0, D13C > -100) %>%
  arrange(Trel)


# overall flux of CO2
# plot pre-drought vs drought faceted by pyruvate label
data_sub %>%
  gather(Variable, Value, -Label, -Condition, -Pyruv, -Trel) %>% 
  mutate(Condition=factor(Condition, levels = c("pre_drought", "drought"))) %>%
  filter(Variable == "Flux") %>% 
  mutate(Site = case_when(
    startsWith(Label, "S1") ~ "Site1",
    startsWith(Label, "S2") ~ "Site2",
    startsWith(Label, "S3") ~ "Site3"
  )) %>%
  ggplot(aes(x = Trel, y = Value, color = Condition)) +
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(span = 0.5) +
  facet_wrap(.~ Pyruv) +
  scale_color_manual(values = colors,
                     breaks = c("pre_drought", "drought"),
                     labels = c("Pre Drought", "Drought")) + 
  labs(x = "Time post pyruvate (hr)", y = bquote(delta^13 *CO[2])) +
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        legend.title = element_blank())
ggsave("Figures/CO2-VOCs/C1_comb_C2_flux-update-pyruv.png", width = 4, height = 3, dpi = 1000)

# plot pre-drought vs drought not faceted
data_sub %>% 
  gather(Variable, Value, -Label, -Condition, -Pyruv, -Trel) %>% 
  mutate(Condition=factor(Condition, levels = c("pre_drought", "drought"))) %>% 
  filter(Variable == "Flux") %>% 
  mutate(label_state = ifelse(Trel >0, "post", "pre")) %>% 
  ggplot(aes(x = Trel, y = Value, color = Condition)) +
  #facet_grid(. ~ Pyruv) + 
  geom_point(alpha = 0.3, size = 2) +
  geom_smooth(span = 0.5) +
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  scale_color_manual(values = colors,
                     breaks = c("pre_drought", "drought"),
                     labels = c("Pre Drought", "Drought")) + 
  labs(x = "Hours after labeling", y = bquote("Soil" ~ CO[2]~"efflux (" *mu*mol ~m^-2*s^-1*")")) +
  theme(text = element_text(size = 17,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        legend.title = element_blank()) +
  ylim(0,6)
dev.off()

# Isotope Sig
# plot pre-drought vs drought faceted by label
data_sub %>%
  gather(Variable, Value, -Label, -Condition, -Pyruv, -Trel) %>% 
  mutate(Condition=factor(Condition, levels = c("pre_drought", "drought"))) %>%
  filter(Variable == "D13C" & Pyruv != "water") %>% 
  mutate(label_state = ifelse(Trel >0, "post", "pre")) %>% 
  mutate(Site = case_when(
    startsWith(Label, "S1") ~ "Site1",
    startsWith(Label, "S2") ~ "Site2",
    startsWith(Label, "S3") ~ "Site3"
  )) %>%
  ggplot(aes(x = Trel, y = Value, color = Condition)) +
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(aes(group = interaction(label_state, Condition)), span = 0.5) +
  facet_wrap( ~ Pyruv) +
  scale_color_manual(values = colors,
                     breaks = c("pre_drought", "drought"),
                     labels = c("Pre Drought", "Drought")) + 
  labs(x = "Time post pyruvate (hr)", y = bquote(delta^13 *CO[2])) +
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        legend.title = element_blank())
ggsave("Figures/CO2-VOCs/C1_C2_delta_co2-update-cond.png", width = 4, height = 3, dpi = 1000)

# plot pre-drought vs drought faceted by site label
data_sub %>%
  gather(Variable, Value, -Label, -Condition, -Pyruv, -Trel) %>% 
  mutate(Condition=factor(Condition, levels = c("pre_drought", "drought"))) %>%
  filter(Variable == "D13C" & Pyruv != "water") %>% 
  mutate(label_state = ifelse(Trel >0, "post", "pre")) %>% 
  mutate(Site = case_when(
    startsWith(Label, "S1") ~ "Site1",
    startsWith(Label, "S2") ~ "Site2",
    startsWith(Label, "S3") ~ "Site3"
  )) %>%
  ggplot(aes(x = Trel, y = Value, color = Condition)) +
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(aes(group = interaction(label_state, Condition)), span = 0.5) +
  facet_wrap(Site ~ Pyruv, ncol = 2) +
  scale_color_manual(values = colors,
                     breaks = c("pre_drought", "drought"),
                    labels = c("Pre Drought", "Drought")) + 
  labs(x = "Time post pyruvate (hr)", y = bquote(delta^13 *CO[2])) +
  theme(text = element_text(size = 10,
                           family = "Arial",
                           color = "black"),
    legend.position = "bottom",
    legend.title = element_blank())
ggsave("Figures/CO2-VOCs/C1_C2_delta_co2-update-cond_by_site.png", width = 4, height = 3, dpi = 1000)

# plot C1 vs C2 faceted by condition and site 
data_sub %>%
  gather(Variable, Value, -Label, -Condition, -Pyruv, -Trel) %>% 
  mutate(Condition=factor(Condition, levels = c("pre_drought", "drought"))) %>%
  filter(Variable == "D13C" & Pyruv != "water") %>% 
  mutate(label_state = ifelse(Trel >0, "post", "pre")) %>% 
  mutate(Site = case_when(
    startsWith(Label, "S1") ~ "Site1",
    startsWith(Label, "S2") ~ "Site2",
    startsWith(Label, "S3") ~ "Site3"
  )) %>%
  ggplot(aes(x = Trel, y = Value, color = Pyruv)) +
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(aes(group = interaction(label_state, Pyruv)), span = 0.5) +
  facet_wrap(Condition ~ Site) +
  scale_color_manual(values = colors_pyruv,
                     breaks = c("C1", "C2")) +
  labs(x = "Time post pyruvate (hr)", y = bquote(delta^13 *CO[2])) +
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        legend.title = element_blank())
ggsave("Figures/CO2-VOCs/C1_C2_delta_co2-update-label_by_cond_by_site.png", width = 5, height = 4, dpi =1000)


# Plotwise
filename=paste("Figures/CO2-VOCs/C1_comb_C2_flux-update-plotwise.png", sep = "")
png(filename ,width=9, height=5, unit='in', res = 1000)

data_sub %>% 
  gather(Variable, Value, -Label, -Condition, -Pyruv, -Trel) %>% 
  mutate(Condition=factor(Condition, levels = c("pre_drought", "drought"))) %>% 
  filter(Variable == "Flux") %>% 
  mutate(label_state = ifelse(Trel >0, "post", "pre")) %>% 
  mutate(Site = case_when(
    startsWith(Label, "S1") ~ "Site1",
    startsWith(Label, "S2") ~ "Site2",
    startsWith(Label, "S3") ~ "Site3"
  )) %>%
  ggplot(aes(x = Trel, y = Value, color = Condition)) +
  facet_wrap(~ Label, ncol = 4) + 
  geom_point(alpha = 0.3, size = 2) +
  geom_smooth(span = 0.5) +
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  scale_color_manual(values = colors,
                     breaks = c("pre_drought", "drought"),
                     labels = c("Pre Drought", "Drought")) + 
  labs(x = "Hours after labeling", y = bquote("Soil" ~ CO[2]~"efflux (" *mu*mol ~m^-2*s^-1*")")) +
  theme(text = element_text(size = 17,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        legend.title = element_blank()) +
  ylim(0,6)
dev.off()


data_sub %>%
  gather(Variable, Value, -Label, -Condition, -Pyruv, -Trel) %>% 
  mutate(Condition=factor(Condition, levels = c("pre_drought", "drought"))) %>%
  filter(Variable == "D13C" & Pyruv == "C1") %>% 
  mutate(label_state = ifelse(Trel >0, "post", "pre")) %>% 
  mutate(Site = case_when(
    startsWith(Label, "S1") ~ "Site1",
    startsWith(Label, "S2") ~ "Site2",
    startsWith(Label, "S3") ~ "Site3"
  )) %>%
  ggplot(aes(x = Trel, y = Value, color = Condition)) +
  facet_wrap(.~ Label) + 
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  geom_point(alpha = 0.3) +
  geom_smooth(aes(group = interaction(label_state, Condition)), span = 0.5) +
  scale_color_manual(values = colors,
                     breaks = c("pre_drought", "drought"),
                     labels = c("Pre Drought", "Drought")) + 
  labs(x = "Hours after labeling", y = bquote(delta^13 *CO[2]), subtitle = "C1-Pyruvate") +
  theme(legend.position = "bottom")
ggsave("./Figures/CO2-VOCs/01_Pyruvate_SoilResp_C1_Plotwise_update.png")

data_sub %>% 
  gather(Variable, Value, -Label, -Condition, -Pyruv, -Trel) %>% 
  mutate(Condition=factor(Condition, levels = c("pre_drought", "drought"))) %>% 
  filter(Variable == "D13C" & Pyruv =="C2") %>% 
  mutate(label_state = ifelse(Trel >0, "post", "pre")) %>% 
  ggplot(aes(x = Trel, y = Value, color = Condition)) +
  facet_wrap(. ~ Label) + 
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  geom_point(alpha = 0.3) +
  geom_smooth(aes(group = interaction(label_state, Condition)), span = 0.5) +
  scale_color_manual(values = colors,
                     breaks = c("pre_drought", "drought"),
                     labels = c("Pre Drought", "Drought")) + 
  labs(x = "Hours after labeling", y = bquote(delta^13 *CO[2]), subtitle = "C2-Pyruvate") +
  theme(legend.position = "bottom")
ggsave("./Figures/CO2-VOCs/01_Pyruvate_SoilResp_C2_Plotwise_update.png")

# Isotope Sig - boxplots
filename=paste("Figures/CO2-VOCs/C1_C2_delta_co2-boxplot.png", sep = "")
png(filename ,width=4, height=4, unit='in', res = 1000)

data_sub %>%
  gather(Variable, Value, -Label, -Condition, -Pyruv, -Trel) %>% 
  mutate(Condition=factor(Condition, levels = c("pre_drought", "drought"))) %>%
  filter(Variable == "D13C") %>% 
  mutate(label_state = ifelse(Trel >0, "post", "pre")) %>% 
  filter(label_state == "post") %>%
  ggplot(aes(x = Condition, y = Value, fill = Condition)) +
  facet_grid(. ~ Pyruv) + 
  geom_boxplot() +
  scale_x_discrete(labels = c("Pre Drought", "Drought")) +
  scale_fill_manual(values = colors,
                     breaks = c("pre_drought", "drought"),
                     labels = c("Pre Drought", "Drought")) + 
  labs(x = "", y = bquote(delta^13 *CO[2])) +
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        legend.title = element_blank())
dev.off()


filename=paste("Figures/CO2-VOCs/C1_comb_C2_flux-barplot.png", sep = "")
png(filename ,width=2.5, height=4, unit='in', res = 1000)

data_sub %>% 
  gather(Variable, Value, -Label, -Condition, -Pyruv, -Trel) %>% 
  mutate(Condition=factor(Condition, levels = c("pre_drought", "drought"))) %>% 
  filter(Variable == "Flux") %>% 
  mutate(label_state = ifelse(Trel >0, "post", "pre")) %>% 
  filter(label_state == "post") %>%
  ggplot(aes(x = Condition, y = Value, fill = Condition)) +
  #facet_grid(. ~ Pyruv) + 
  geom_boxplot() +
  scale_x_discrete(labels = c("Pre Drought", "Drought")) +
  scale_fill_manual(values = colors,
                     breaks = c("pre_drought", "drought"),
                     labels = c("Pre Drought", "Drought")) + 
  labs(y = bquote("Soil" ~ CO[2]~"efflux (" *mu*mol ~m^-2*s^-1*")")) +
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        legend.title = element_blank()) +
  ylim(0,6)
dev.off()

###C1/C2 ratios, separate by site and divide C1:C2 ratios of nearby locations (paired locations)
# modify data_sub table to add time intervals and calculate C1 and C2 ratios per site 
ratios_per_set <- data_sub %>%
  mutate(label_state = ifelse(Trel >0, "post", "pre")) %>% 
  filter(label_state == "post") %>%
  filter(Pyruv != "water") %>%
  select(-label_state) %>%
  #filter(Label != "S2.P1") %>%
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
  select(-Trel) %>%
  group_by(Condition, Site, Set, Pyruv, Time) %>%
  summarise_all(funs(mean)) %>%
  pivot_wider(names_from = Pyruv, values_from = c(D13C, Flux)) %>%
  mutate(C1_C2_ratio = Flux_C1/Flux_C2) %>%
  mutate(C1_C2_diff = Flux_C1 - Flux_C2) %>%
  mutate(C1_C2_d13c_ratio = D13C_C1/D13C_C2) %>%
  mutate(C1_C2_d13c_diff = D13C_C1 - D13C_C2)

# graph of all differences
ratios_per_set %>%
  filter(Set != "Set4") %>%
  ggplot(aes(x = Time, y = C1_C2_d13c_diff, fill = Condition)) +
  geom_boxplot() +
  labs(x = "Time post pyruvate (hr)", y = bquote(delta^13 *CO[2] * (C1 - C2))) +
  #facet_wrap(. ~ Site) + 
  scale_fill_manual(values = colors,
                    breaks = c("pre_drought", "drought"),
                    labels = c("Pre Drought", "Drought")) + 
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        legend.title = element_blank())
ggsave("Figures/CO2-VOCs/C1-C2-difference-all-sites.png", width = 5, height = 3, dpi = 1000)



####statistical analysis
#Flux - Condition, in C1 and C2 separated
data_sub_c1 <- data_sub %>%
 filter(Pyruv == "C1")

data_sub_c2 <- data_sub %>%
  filter(Pyruv == "C2")


wilcox.test(Flux ~ Condition, data = data_sub_c1)

wilcox.test(Flux ~ Condition, data = data_sub_c2)

#Flux - Condition, in C1 and C2 combined

wilcox.test(Flux ~ Condition, data = data_sub)
#data:  Flux by Condition
#W = 786373, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0

t.test(Flux ~ Condition, data = data_sub)
#Welch Two Sample t-test
#
#data:  Flux by Condition
#t = 8.7988, df = 2184.3, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.5834676 0.9181427
#sample estimates:
#  mean in group pre_drought     mean in group drought 
#3.200523                  2.449717 

#Delta - Condition, in C1 and C2 separated
wilcox.test(D13C ~ Condition, data = data_sub_c1)
#data:  D13C by Condition
#W = 193599, p-value = 1.712e-11
#alternative hypothesis: true location shift is not equal to 0

wilcox.test(D13C ~ Condition, data = data_sub_c2)
#data:  D13C by Condition
#W = 178428, p-value = 0.0001037
#alternative hypothesis: true location shift is not equal to 0

#Flux - Condition, in C1 and C2 combined

wilcox.test(D13C ~ Condition, data = data_sub)
#data:  D13C by Condition
#W = 716268, p-value = 1.255e-08
#alternative hypothesis: true location shift is not equal to 0

t.test(D13C ~ Condition, data = data_sub)
#data:  D13C by Condition
#t = 11.216, df = 1828.7, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  1138.909 1621.638
#sample estimates:
#  mean in group pre_drought     mean in group drought 
#3176.336                  1796.062 

# ratios for differences in C1 - C2 deltaC-CO2r
#t-test
ratios_per_set_3hr<- ratios_per_set %>%
  filter(Time == "3" )
t.test(C1_C2_d13c_diff ~ Condition, data = ratios_per_set_3hr )

ratios_per_set_6hr<- ratios_per_set %>%
  filter(Time == "6" )
t.test(C1_C2_d13c_diff ~ Condition, data = ratios_per_set_6hr ) #p = 0.0176

ratios_per_set_12hr<- ratios_per_set %>%
  filter(Time == "12" )
t.test(C1_C2_d13c_diff ~ Condition, data = ratios_per_set_12hr ) #p=0.00418

ratios_per_set_18hr<- ratios_per_set %>%
  filter(Time == "18" )
t.test(C1_C2_d13c_diff ~ Condition, data = ratios_per_set_18hr ) #p=0.00158

ratios_per_set_24hr<- ratios_per_set %>%
  filter(Time == "24" )
t.test(C1_C2_d13c_diff ~ Condition, data = ratios_per_set_24hr ) #p=0.0371

ratios_per_set_30hr<- ratios_per_set %>%
  filter(Time == "30" )
t.test(C1_C2_d13c_diff ~ Condition, data = ratios_per_set_30hr )

ratios_per_set_36hr<- ratios_per_set %>%
  filter(Time == "36" )
t.test(C1_C2_d13c_diff ~ Condition, data = ratios_per_set_36hr )

ratios_per_set_42hr<- ratios_per_set %>%
  filter(Time == "42" )
t.test(C1_C2_d13c_diff ~ Condition, data = ratios_per_set_42hr )

ratios_per_set_48hr<- ratios_per_set %>%
  filter(Time == "48" )
t.test(C1_C2_d13c_diff ~ Condition, data = ratios_per_set_48hr )

#wilcoxan test
ratios_per_set_3hr<- ratios_per_set %>%
  filter(Time == "3" )
wilcox.test(C1_C2_d13c_diff ~ Condition, data = ratios_per_set_3hr )

ratios_per_set_6hr<- ratios_per_set %>%
  filter(Time == "6" )
wilcox.test(C1_C2_d13c_diff ~ Condition, data = ratios_per_set_6hr ) #p = 0.0176

ratios_per_set_12hr<- ratios_per_set %>%
  filter(Time == "12" )
wilcox.test(C1_C2_d13c_diff ~ Condition, data = ratios_per_set_12hr ) #p=0.00418

ratios_per_set_18hr<- ratios_per_set %>%
  filter(Time == "18" )
wilcox.test(C1_C2_d13c_diff ~ Condition, data = ratios_per_set_18hr ) #p=0.00158

ratios_per_set_24hr<- ratios_per_set %>%
  filter(Time == "24" )
wilcox.test(C1_C2_d13c_diff ~ Condition, data = ratios_per_set_24hr ) #p=0.0371

ratios_per_set_30hr<- ratios_per_set %>%
  filter(Time == "30" )
wilcox.test(C1_C2_d13c_diff ~ Condition, data = ratios_per_set_30hr )

ratios_per_set_36hr<- ratios_per_set %>%
  filter(Time == "36" )
wilcox.test(C1_C2_d13c_diff ~ Condition, data = ratios_per_set_36hr )

ratios_per_set_42hr<- ratios_per_set %>%
  filter(Time == "42" )
wilcox.test(C1_C2_d13c_diff ~ Condition, data = ratios_per_set_42hr )

ratios_per_set_48hr<- ratios_per_set %>%
  filter(Time == "48" )
wilcox.test(C1_C2_d13c_diff ~ Condition, data = ratios_per_set_48hr )
