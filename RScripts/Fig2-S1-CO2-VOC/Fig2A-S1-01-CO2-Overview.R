
# Visualize Pyruvate Co2 Fluxes

# Johannes Ingrisch, modified by Linnea Honeker
# 2021-08-26

# load packages
library(tidyverse)
library(plotly)
library(ggthemes)
library(patchwork)

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

# factor Condition
data$Condition <- factor(data$Condition, c(
  "pre_drought", "drought"))

# Isotope Sig
filename=paste("Figures/Fig2-S1-CO2-VOCs/Fig2A-C1_C2_delta_co2-update.png", sep = "")
png(filename ,width=4, height=4, unit='in', res = 1000)

data_sub %>%
  gather(Variable, Value, -Label, -Condition, -Pyruv, -Trel) %>% 
  mutate(Condition=factor(Condition, levels = c("pre_drought", "drought"))) %>%
  filter(Variable == "D13C") %>% 
  mutate(label_state = ifelse(Trel >0, "post", "pre")) %>% 
  ggplot(aes(x = Trel, y = Value, color = Condition)) +
  facet_grid(. ~ Pyruv) + 
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(aes(group = interaction(label_state, Condition)), span = 0.5) +
  scale_color_manual(values = colors,
                     breaks = c("pre_drought", "drought"),
                     labels = c("Pre Drought", "Drought")) + 
  labs(x = "Time post pyruvate (h)", y = expression(delta^13 *"C"[CO[2]]*" (\u2030, VPDB)")) +
  theme(text = element_text(size = 12,
                           family = "Arial",
                           color = "black"),
    legend.position = "bottom",
    legend.title = element_blank())
dev.off()


filename=paste("Figures/Fig2-S1-CO2-VOCs/FigS1-C1_comb_C2_flux-update.png", sep = "")
png(filename ,width=9, height=5, unit='in', res = 1000)

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
  labs(x = "Time post pyruvate (h)", y = bquote("Soil" ~ CO[2]~"efflux (" *mu*mol ~m^-2*s^-1*")")) +
  theme(text = element_text(size = 17,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        legend.title = element_blank()) +
  ylim(0,6)
dev.off()


# Plotwise
data_sub %>%
  gather(Variable, Value, -Label, -Condition, -Pyruv, -Trel) %>% 
  mutate(Condition=factor(Condition, levels = c("pre_drought", "drought"))) %>%
  filter(Variable == "D13C" & Pyruv == "C1") %>% 
  mutate(label_state = ifelse(Trel >0, "post", "pre")) %>% 
  ggplot(aes(x = Trel, y = Value, color = Condition)) +
  facet_wrap(. ~ Label) + 
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  geom_point(alpha = 0.3) +
  geom_smooth(aes(group = interaction(label_state, Condition)), span = 0.5) +
  scale_color_manual(values = colors,
                     breaks = c("pre_drought", "drought"),
                     labels = c("Pre Drought", "Drought")) + 
  labs(x = "Hours after labeling", y = bquote(delta^13 *CO[2]), subtitle = "C1-Pyruvate") +
  theme(legend.position = "bottom")
ggsave("./Figures/Fig2-S1-CO2-VOCs/01_Pyruvate_SoilResp_C1_Plotwise_update.png")

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
ggsave("./Figures/Fig2-S1-CO2-VOCs/01_Pyruvate_SoilResp_C2_Plotwise_update.png")

####statistical analysis, comparing differences in treatments based on mean flux data
#Flux - Condition, in C1 and C2 separated
data_sub_c1 <- data_sub %>%
 filter(Pyruv == "C1") %>%
  filter(Trel > 0)

data_sub_c2 <- data_sub %>%
  filter(Pyruv == "C2")%>%
  filter(Trel > 0)

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