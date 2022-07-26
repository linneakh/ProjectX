
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

###plot just 13C-CO2 from pyruvate over time
data_sub_13C_flux %>%
  mutate(Condition=factor(Condition, levels = c("pre_drought", "drought"))) %>%
  mutate(Site = case_when(
    startsWith(Label, "S1") ~ "Site1",
    startsWith(Label, "S2") ~ "Site2",
    startsWith(Label, "S3") ~ "Site3"
  )) %>%
  filter(fCO2_13_py < 1) %>%
  ggplot(aes(x = Trel, y = fCO2_13_py, color = Condition)) +
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(aes(group = interaction(Condition)), span = 0.5) +
  facet_wrap(~ Pyruv) +
  #ylim(0, 0.5) +
  scale_color_manual(values = colors,
                     breaks = c("pre_drought", "drought"),
                     labels = c("Pre Drought", "Drought")) + 
  labs(x = "Time post pyruvate (hr)", y= bquote("Soil 13C-" ~ CO[2]~"efflux (" *mu*mol ~m^-2*s^-1*")")) +
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        legend.title = element_blank())
ggsave("Figures/CO2-VOCs/C1_C2_13C_CO2_f_py.png", width = 4, height = 3, dpi = 1000)


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
  #filter(Set != "Set4") %>%
  ggplot(aes(x = Time, y = C1_C2_d13c_diff, fill = Condition)) +
  geom_boxplot() +
  labs(x = "Time post pyruvate (hr)", y = bquote(delta^13 *CO[2] * (C1 - C2))) +
  facet_wrap(. ~ Site) + 
  scale_fill_manual(values = colors,
                    breaks = c("pre_drought", "drought"),
                    labels = c("Pre Drought", "Drought")) + 
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        legend.title = element_blank())
ggsave("Figures/CO2-VOCs/C1-C2-difference-all-sites_w_set4.png", width = 4, height = 4, dpi = 1000)

###C1/C2 ratios, separate by site and divide C1:C2 ratios of nearby locations (paired locations)
# modify data_sub table to add time intervals and calculate C1 and C2 ratios per site 
ratios_per_set <- data_sub_13C_flux %>%
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
  select(-c(Trel, D13C, Flux, AtF, AtFe)) %>%
  group_by(Condition, Site, Set, Pyruv, Time) %>%
  summarise_all(funs(mean)) %>%
  pivot_wider(names_from = Pyruv, values_from = c(fCO2_13, fCO2_13_py)) %>%
  mutate(C1_C2_diff_fCO2_13_py = fCO2_13_py_C1 - fCO2_13_py_C2) %>%
  mutate(C1_C2_ratio_sum_fCO2_13_py = fCO2_13_py_C1/(fCO2_13_py_C1 + fCO2_13_py_C2))

ratios_per_set$Condition <- factor(ratios_per_set$Condition, c("pre_drought", "drought"))

# graph of all differences (flux13)
ratios_per_set %>%
  #filter(Set != "Set4") %>%
  ggplot(aes(x = Time, y = C1_C2_diff_fCO2_13, fill = Condition)) +
  geom_boxplot() +
  labs(x = "Time post pyruvate (hr)", y= bquote("Soil 13C-" ~ CO[2]~"efflux (C1 - C2)")) +
  #facet_wrap(. ~ Site) + 
  scale_fill_manual(values = colors,
                    breaks = c("pre_drought", "drought"),
                    labels = c("Pre Drought", "Drought")) + 
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        legend.title = element_blank())
ggsave("Figures/CO2-VOCs/C1-C2-difference-13C-CO2.png", width = 5, height = 3, dpi = 1000)

# graph of all ratios (flux13)
ratios_per_set %>%
  filter(Set != "Set4") %>%
  ggplot(aes(x = Time, y = C1_C2_ratio_fCO2_13, fill = Condition)) +
  geom_boxplot() +
  labs(x = "Time post pyruvate (hr)", y= bquote("Soil 13C-" ~ CO[2]~"efflux (C1/C2)")) +
  #facet_wrap( Site~ Set) + 
  scale_fill_manual(values = colors,
                    breaks = c("pre_drought", "drought"),
                    labels = c("Pre Drought", "Drought")) + 
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        legend.title = element_blank())
ggsave("Figures/CO2-VOCs/C1-C2-ratio-13C-CO2.png", width = 5, height = 3, dpi = 1000)

# graph of all differences (flux13Conly)
ratios_per_set %>%
  #filter(Set != "Set4") %>%
  ggplot(aes(x = Time, y = C1_C2_diff_fCO2_13_py, fill = Condition)) +
  geom_boxplot() +
  labs(x = "Time post pyruvate (h)", y= expression(""^13*"C-CO"[2-C1]*" - "^13*"C-CO"[2-C2]*" ("*mu*"mol m"^-2*"s"^-1*")")) +
  #facet_wrap(. ~ Site) + 
  scale_fill_manual(values = colors,
                    breaks = c("pre_drought", "drought"),
                    labels = c("Pre Drought", "Drought")) + 
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        legend.title = element_blank())
ggsave("Figures/CO2-VOCs/C1-C2-difference-13C-CO2-py-w-set4.png", width = 4, height = 4, dpi = 1000)

# graph of all ratios (flux13Conly)
ratios_per_set %>%
  filter(Set != "Set4") %>%
  filter(Site != "Site2" |
           Set != "Set1", .preserve = TRUE) %>%
  ggplot(aes(x = Time, y = C1_C2_ratio_fCO2_13_py, fill = Condition)) +
  geom_boxplot() +
  labs(x = "Time post pyruvate (hr)", y= bquote("Soil 13C-" ~ CO[2]~"efflux (C1/C2)")) +
  #facet_wrap(~ Site) + 
  scale_fill_manual(values = colors,
                    breaks = c("pre_drought", "drought"),
                    labels = c("Pre Drought", "Drought")) + 
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        legend.title = element_blank())
ggsave("Figures/CO2-VOCs/C1-C2-ratio-13C-CO2-py.png", width = 5, height = 3, dpi = 1000)

# graph of all ratios (flux 13C from pyruvate, C1/(C1+C2))
ratios_per_set %>%
  #filter(Set != "Set4") %>%
  ggplot(aes(x = Time, y = C1_C2_ratio_sum_fCO2_13_py, fill = Condition)) +
  geom_boxplot() +
  labs(x = "Time post pyruvate (h)", y= expression(""^13*"C-CO"[2-C1]*" / ("^13*"C-CO"[2-C1]*" + "^13*"C-CO"[2-C2]*") ")) +
  #facet_wrap(~Site) + 
  scale_fill_manual(values = colors,
                    breaks = c("pre_drought", "drought"),
                    labels = c("Pre Drought", "Drought")) + 
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        legend.title = element_blank())
ggsave("Figures/CO2-VOCs/C1-C2-ratio-sum-13C-CO2-py-w-set4.png", width = 5, height = 4, dpi = 1000)

####statistical analysis
#Flux - Condition, in C1 and C2 separated
data_sub_c1 <- data_sub %>%
 filter(Pyruv == "C1")

data_sub_c2 <- data_sub %>%
  filter(Pyruv == "C2")

data_sub_pre_drought <- data_sub %>%
  filter(Condition == "pre_drought")

data_sub_drought <- data_sub %>%
  filter(Condition == "drought")


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

# d13C in predrought between C1 and c2
t.test(D13C ~ Pyruv, data = data_sub_pre_drought)
#data:  D13C by Pyruv
#t = 16.112, df = 746.32, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  2764.003 3531.028
#sample estimates:
#  mean in group C1 mean in group C2 
#4740.906         1593.390 

wilcox.test(D13C ~ Pyruv, data = data_sub_pre_drought)
#data:  D13C by Pyruv
#W = 245832, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0

# d13C in drought between C1 and c2
t.test(D13C ~ Pyruv, data = data_sub_drought)
#data:  D13C by Pyruv
#t = 11.077, df = 811.59, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  1017.691 1456.069
#sample estimates:
#  mean in group C1 mean in group C2 
#2418.629         1181.749 

wilcox.test(D13C ~ Pyruv, data = data_sub_drought)
#data:  D13C by Pyruv
#W = 190837, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0


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

# ratios for differences in C1 - C2 deltaC-CO2r - new calculations
#t-test
ratios_per_set_3hr<- ratios_per_set %>%
  filter(Time == "3" )
t.test(C1_C2_diff_fCO2_13 ~ Condition, data = ratios_per_set_3hr )

ratios_per_set_6hr<- ratios_per_set %>%
  filter(Time == "6" )
t.test(C1_C2_diff_fCO2_13 ~ Condition, data = ratios_per_set_6hr ) #p = 0.02934

ratios_per_set_12hr<- ratios_per_set %>%
  filter(Time == "12" )
t.test(C1_C2_diff_fCO2_13 ~ Condition, data = ratios_per_set_12hr ) #p=0.00742

ratios_per_set_18hr<- ratios_per_set %>%
  filter(Time == "18" )
t.test(C1_C2_diff_fCO2_13 ~ Condition, data = ratios_per_set_18hr ) #p=0.003106

ratios_per_set_24hr<- ratios_per_set %>%
  filter(Time == "24" )
t.test(C1_C2_diff_fCO2_13 ~ Condition, data = ratios_per_set_24hr )

ratios_per_set_30hr<- ratios_per_set %>%
  filter(Time == "30" )
t.test(C1_C2_diff_fCO2_13 ~ Condition, data = ratios_per_set_30hr )

ratios_per_set_36hr<- ratios_per_set %>%
  filter(Time == "36" )
t.test(C1_C2_diff_fCO2_13 ~ Condition, data = ratios_per_set_36hr )

ratios_per_set_42hr<- ratios_per_set %>%
  filter(Time == "42" )
t.test(C1_C2_diff_fCO2_13 ~ Condition, data = ratios_per_set_42hr )

ratios_per_set_48hr<- ratios_per_set %>%
  filter(Time == "48" )
t.test(C1_C2_diff_fCO2_13 ~ Condition, data = ratios_per_set_48hr )

#wilcoxan test
ratios_per_set_3hr<- ratios_per_set %>%
  filter(Time == "3" )
wilcox.test(C1_C2_diff_fCO2_13 ~ Condition, data = ratios_per_set_3hr )

ratios_per_set_6hr<- ratios_per_set %>%
  filter(Time == "6" )
wilcox.test(C1_C2_diff_fCO2_13 ~ Condition, data = ratios_per_set_6hr ) 

ratios_per_set_12hr<- ratios_per_set %>%
  filter(Time == "12" )
wilcox.test(C1_C2_diff_fCO2_13 ~ Condition, data = ratios_per_set_12hr ) #p=0.003702

ratios_per_set_18hr<- ratios_per_set %>%
  filter(Time == "18" )
wilcox.test(C1_C2_diff_fCO2_13 ~ Condition, data = ratios_per_set_18hr ) #p=0.007898

ratios_per_set_24hr<- ratios_per_set %>%
  filter(Time == "24" )
wilcox.test(C1_C2_diff_fCO2_13 ~ Condition, data = ratios_per_set_24hr ) 

ratios_per_set_30hr<- ratios_per_set %>%
  filter(Time == "30" )
wilcox.test(C1_C2_diff_fCO2_13 ~ Condition, data = ratios_per_set_30hr )

ratios_per_set_36hr<- ratios_per_set %>%
  filter(Time == "36" )
wilcox.test(C1_C2_diff_fCO2_13 ~ Condition, data = ratios_per_set_36hr )

ratios_per_set_42hr<- ratios_per_set %>%
  filter(Time == "42" )
wilcox.test(C1_C2_diff_fCO2_13 ~ Condition, data = ratios_per_set_42hr )

ratios_per_set_48hr<- ratios_per_set %>%
  filter(Time == "48" )
wilcox.test(C1_C2_diff_fCO2_13 ~ Condition, data = ratios_per_set_48hr )

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

# C1/C2 ratios of 13C-CO2-from-pyruvate- statistical tests
#t-test
ratios_per_set_3hr<- ratios_per_set %>%
  filter(Time == "3" )
t.test(C1_C2_ratio_fCO2_13_py ~ Condition, data = ratios_per_set_3hr )

ratios_per_set_6hr<- ratios_per_set %>%
  filter(Time == "6" )
t.test(C1_C2_ratio_fCO2_13_py ~ Condition, data = ratios_per_set_6hr ) #p = 0.004931

ratios_per_set_12hr<- ratios_per_set %>%
  filter(Time == "12" )
t.test(C1_C2_ratio_fCO2_13_py ~ Condition, data = ratios_per_set_12hr ) #p=0.001838

ratios_per_set_18hr<- ratios_per_set %>%
  filter(Time == "18" )
t.test(C1_C2_ratio_fCO2_13_py ~ Condition, data = ratios_per_set_18hr ) #p=0.0004503

ratios_per_set_24hr<- ratios_per_set %>%
  filter(Time == "24" )
t.test(C1_C2_ratio_fCO2_13_py ~ Condition, data = ratios_per_set_24hr ) #0.03411

ratios_per_set_30hr<- ratios_per_set %>%
  filter(Time == "30" )
t.test(C1_C2_ratio_fCO2_13_py ~ Condition, data = ratios_per_set_30hr )

ratios_per_set_36hr<- ratios_per_set %>%
  filter(Time == "36" )
t.test(C1_C2_ratio_fCO2_13_py ~ Condition, data = ratios_per_set_36hr )

ratios_per_set_42hr<- ratios_per_set %>%
  filter(Time == "42" )
t.test(C1_C2_ratio_fCO2_13_py ~ Condition, data = ratios_per_set_42hr )

ratios_per_set_48hr<- ratios_per_set %>%
  filter(Time == "48" )
t.test(C1_C2_ratio_fCO2_13_py ~ Condition, data = ratios_per_set_48hr )

#wilcoxan test
ratios_per_set_3hr<- ratios_per_set %>%
  filter(Time == "3" )
wilcox.test(C1_C2_ratio_fCO2_13_py ~ Condition, data = ratios_per_set_3hr ) # 0.03595

ratios_per_set_6hr<- ratios_per_set %>%
  filter(Time == "6" )
wilcox.test(C1_C2_ratio_fCO2_13_py ~ Condition, data = ratios_per_set_6hr ) # 0.007898

ratios_per_set_12hr<- ratios_per_set %>%
  filter(Time == "12" )
wilcox.test(C1_C2_ratio_fCO2_13_py ~ Condition, data = ratios_per_set_12hr ) #p=0.001563

ratios_per_set_18hr<- ratios_per_set %>%
  filter(Time == "18" )
wilcox.test(C1_C2_ratio_fCO2_13_py ~ Condition, data = ratios_per_set_18hr ) #p=0.002468

ratios_per_set_24hr<- ratios_per_set %>%
  filter(Time == "24" )
wilcox.test(C1_C2_ratio_fCO2_13_py ~ Condition, data = ratios_per_set_24hr ) #0.02065

ratios_per_set_30hr<- ratios_per_set %>%
  filter(Time == "30" )
wilcox.test(C1_C2_ratio_fCO2_13_py ~ Condition, data = ratios_per_set_30hr )

ratios_per_set_36hr<- ratios_per_set %>%
  filter(Time == "36" )
wilcox.test(C1_C2_ratio_fCO2_13_py ~ Condition, data = ratios_per_set_36hr )

ratios_per_set_42hr<- ratios_per_set %>%
  filter(Time == "42" )
wilcox.test(C1_C2_ratio_fCO2_13_py ~ Condition, data = ratios_per_set_42hr )

ratios_per_set_48hr<- ratios_per_set %>%
  filter(Time == "48" )
wilcox.test(C1_C2_ratio_fCO2_13_py ~ Condition, data = ratios_per_set_48hr )

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

