
# Visualize Pyruvate voc Fluxes

# Linnea Honeker
# 2021-08-26, modified 4/1/22

# load packages
library(tidyverse)
library(plotly)
library(ggthemes)
library(patchwork)

theme_custom <- theme_few()

theme_set(theme_custom)

colors = c("#2ECC71" ,"#F39C12")

data.13C <- read.csv("./Output/CO2-VOCs/diacetyl_flux-13C-all.csv") 
data.13C.12C <- read.csv("./Output/CO2-VOCs/diacetyl_flux-13C-12C-all.csv")

# re order factors for condition
data.13C.12C$Condition <- factor(data.13C.12C$Condition, c(
  "Pre-Drought", "Drought"))

# filter out times outside of -12 hour to 48 hour and water samples
data_sub <- data.13C.12C %>%
  filter(Time > -0.3 & Time <2 ) %>%
  filter(Pyr == "C1" |
           Pyr == "C2")

# Isotope Sig diacetyl
png(filename ,width=4, height=4, unit='in', res = 1000)

data_sub %>%
  mutate(Time_hours = Time*24) %>%
  ggplot(aes(x = Time_hours, y = Flux, color = Condition)) +
  facet_grid(. ~ Pyr) + 
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  #ylim(-0.03, 0.12) +
  geom_point(alpha = 0.3, size = 1) +
  geom_smooth(aes(group = Condition), span = 0.5) +
  scale_color_manual(values = colors,
                     breaks = c("Pre-Drought", "Drought"),
                     labels = c("Pre Drought", "Drought")) + 
  labs(x = "Time post pyruvate (h)", y = expression(""^13*"C/("^12*"C + "^13*"C) flux (m"^-2*" h"^-1*")")) +
  theme(text = element_text(size = 12,
                            color = "black"),
        legend.position = "bottom",
        legend.title = element_blank())
filename=paste("./Figures/Fig2-S2-CO2-VOCs/corrected_voc_data/Fig2F-diacetyl-corrected-13C-12C.pdf", sep = "")
ggsave(filename ,width=4, height=4, units='in', dpi = 300)



# Plotwise by chamber
data_sub %>%
  mutate(Time_hours = Time*24) %>%
  filter(Pyr == "C1") %>%
  ggplot(aes(x = Time_hours, y = Flux, color = Condition)) +
  facet_wrap(. ~ Label, ncol = 5) + 
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  geom_point(alpha = 0.3) +
  geom_smooth(aes(group =  Condition), span = 0.5) +
  scale_color_manual(values = colors,
                     breaks = c("Pre-Drought", "Drought"),
                     labels = c("Pre Drought", "Drought")) + 
  labs(x = "Hours after labeling", y = "13C/(12C + 13C) flux", subtitle = "C1 pyurvate") +
  theme(legend.position = "bottom")
ggsave("./Figures/Fig2-S2-CO2-VOCs/corrected_voc_data/01_pyruvate-diacetyl_C1_Plotwise.png")

data_sub %>%
  mutate(Time_hours = Time*24) %>%
  filter(Pyr == "C2") %>%
  ggplot(aes(x = Time_hours, y = Flux, color = Condition)) +
  facet_wrap(. ~ Label, ncol = 5) + 
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  geom_point(alpha = 0.3) +
  geom_smooth(aes(group = Condition), span = 0.5) +
  scale_color_manual(values = colors,
                     breaks = c("Pre-Drought", "Drought"),
                     labels = c("Pre Drought", "Drought")) + 
  labs(x = "Hours after labeling", y = "13C/(12C + 13C) flux", subtitle = "C2 pyurvate") +
  theme(legend.position = "bottom")
ggsave("./Figures/Fig2-S2-CO2-VOCs/corrected_voc_data/01_Pyruvate_diacetyl_C2_Plotwise.png")


####statistical analysis, comparing differences in treatments based on mean flux data
#Flux - Condition, in C1 and C2 separated
data_c1 <- data_sub %>%
 filter(Pyr == "C1") %>%
  filter(Time > 0)

data_c2 <- data_sub %>%
  filter(Pyr == "C2") %>%
  filter(Time > 0) 

wilcox.test(Flux ~ Condition, data = data_c1, paired = FALSE) #there are an uneven number of samples so paired = TRUE returns an error
#data:  Flux by Condition
#W = 75848, p-value 0.00659
#alternative hypothesis: true location shift is not equal to 0


wilcox.test(Flux ~ Condition, data = data_c2, paired = FALSE)
#data:  Flux by Condition
#W = 40978, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0

#Flux - Condition, in C1 and C2 combined

wilcox.test(Flux ~ Condition, data = data_sub, paired = FALSE)
#data:  Flux by Condition
#W = 180858, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0

t.test(Flux ~ Condition, data = data_sub)
#data:  Flux by Condition
#t =-23.028, df = 899.09, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -78.68654 -66.32762
#sample estimates:
#  mean in group pre-drought     mean in group drought 
# -1.830744                 70.676335 

