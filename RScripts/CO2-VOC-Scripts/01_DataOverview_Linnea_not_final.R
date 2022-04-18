
# Visualize Pyruvate Co2 Fluxes

# Johannes Ingrisch, modified by Linnea Honeker
# 2021-08-26

# load packages
library(tidyverse)
library(plotly)
library(ggthemes)
library(patchwork)

theme_custom <- theme_few()

theme_set(theme_custom)

colors = c("#2ECC71" ,"#F39C12")


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
  filter(Trel > -10, Trel < 50, Flux > 0, D13C > -100)

# Isotope Sig
filename=paste("Figures/CO2-VOCs/C1_C2_delta_co2-update.png", sep = "")
png(filename ,width=5, height=3, unit='in', res = 1000)

data_sub %>%
  gather(Variable, Value, -Label, -Condition, -Pyruv, -Trel) %>% 
  mutate(Condition=factor(Condition, levels = c("pre_drought", "drought"))) %>%
  filter(Variable == "D13C") %>% 
  mutate(label_state = ifelse(Trel >0, "post", "pre")) %>% 
  mutate(Site = case_when(
    startsWith(Label, "S1") ~ "Site1",
    startsWith(Label, "S2") ~ "Site2",
    startsWith(Label, "S3") ~ "Site3"
  )) %>%
  ggplot(aes(x = Trel, y = Value, color = Pyruv)) +
  facet_grid(Site~Condition) + 
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(aes(group = interaction(label_state, Pyruv)), span = 0.5) +
  #scale_color_manual(values = colors,
  #                   breaks = c("pre_drought", "drought"),
  #                   labels = c("Pre Drought", "Drought")) + 
  labs(x = "Hours after labeling", y = bquote(delta^13 *CO[2])) +
  theme(text = element_text(size = 12,
                           family = "Arial",
                           color = "black"),
    legend.position = "bottom",
    legend.title = element_blank())
dev.off()


filename=paste("Figures/CO2-VOCs/C1_comb_C2_flux-update.png", sep = "")
png(filename ,width=9, height=5, unit='in', res = 1000)

data_sub %>% 
  gather(Variable, Value, -Label, -Condition, -Pyruv, -Trel) %>% 
  mutate(Condition=factor(Condition, levels = c("pre_drought", "drought"))) %>% 
  filter(Variable == "Flux") %>% 
  filter(Label != "S2.P1") %>%
  mutate(label_state = ifelse(Trel >0, "post", "pre")) %>% 
  mutate(Site = case_when(
    startsWith(Label, "S1") ~ "Site1",
    startsWith(Label, "S2") ~ "Site2",
    startsWith(Label, "S3") ~ "Site3"
  )) %>%
  ggplot(aes(x = Trel, y = Value, color = Condition)) +
  facet_grid(Site ~ Pyruv) + 
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

###C1/C2 ratios
data_sub %>%
  mutate(label_state = ifelse(Trel >0, "post", "pre")) %>% 
  filter(label_state == "post") %>%
  select(-(c(Trel, label_state))) %>%
  filter(Label != "S2.P1") %>%
  mutate(Site = case_when(
    startsWith(Label, "S1") ~ "Site1",
    startsWith(Label, "S2") ~ "Site2",
    startsWith(Label, "S3") ~ "Site3"
  )) %>%
  group_by(Site, Condition, Pyruv) %>%
  summarise_all(funs(mean)) %>%
  select(-Label) %>%
  pivot_wider(names_from = Pyruv, values_from = c(Flux, D13C)) %>%
  mutate(C1_C2_ratio = Flux_C1/Flux_C2) %>%
  mutate(C1_C2_d13c_ratio = D13C_C1/D13C_C2) %>%
  ggplot(aes(x=Condition, y=C1_C2_d13c_ratio, fill = Condition)) +
  geom_boxplot() +
  #facet_grid(~Site) +
  scale_x_discrete(labels = c("Pre Drought", "Drought")) +
  scale_fill_manual(values = colors,
                    breaks = c("pre_drought", "drought"),
                    labels = c("Pre Drought", "Drought")) + 
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        legend.title = element_blank()) 
 
data_sub %>%
  mutate(label_state = ifelse(Trel >0, "post", "pre")) %>% 
  filter(label_state == "post") %>%
  select(-label_state) %>%
  filter(Label != "S2.P1") %>%
  mutate(Site = case_when(
    startsWith(Label, "S1") ~ "Site1",
    startsWith(Label, "S2") ~ "Site2",
    startsWith(Label, "S3") ~ "Site3"
  )) %>%
  mutate(Time = case_when(
    (Trel < 5) ~ "5hr",
    (Trel > 5 & Trel < 10) ~ "10hr",
    (Trel > 10 & Trel < 20) ~ "20hr",
    (Trel > 20 & Trel < 40) ~ "40hr",
    (Trel > 40 & Trel < 48) ~ "48hr"
  )) %>%
  mutate(Time=factor(Time, levels = c("5hr", "10hr", "20hr", "40hr", "48hr"))) %>% 
  group_by(Site, Condition, Pyruv, Time) %>%
  summarise_all(funs(mean)) %>%
  select(-c(Label, Trel)) %>%
  pivot_wider(names_from = Pyruv, values_from = c(Flux, D13C)) %>%
  mutate(C1_C2_ratio = Flux_C1/Flux_C2) %>%
  mutate(C1_C2_d13c_ratio = D13C_C1/D13C_C2) %>%
  ggplot(aes(x = Time, y = C1_C2_d13c_ratio, fill = Condition)) +
  geom_boxplot() +
  #facet_wrap(. ~ Site) + 
  scale_fill_manual(values = colors,
                    breaks = c("pre_drought", "drought"),
                    labels = c("Pre Drought", "Drought")) + 
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        legend.title = element_blank()) 

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