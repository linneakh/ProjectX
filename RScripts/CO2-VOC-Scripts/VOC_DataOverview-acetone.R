

# Visualize Pyruvate voc Fluxes

# Johannes Ingrisch, modified by Linnea Honeker
# 2021-08-26, modified 4/1/22

# load packages
library(tidyverse)
library(plotly)
library(ggthemes)
library(patchwork)

theme_custom <- theme_few()

theme_set(theme_custom)

colors = c("#2ECC71" ,"#F39C12")

data <- read.csv("./Data/voc_raw/acetone-for-ggplot.csv") 




# re order factors for condition

data$Condition <- factor(data$Condition, c(
  "pre-drought", "drought"))

# filter out times outside of -12 hour to 48 hour

data_sub <- data %>%
  filter(Time > -0.3 & Time <2 )

# Isotope Sig acetate
filename=paste("Figs/acetone.png", sep = "")
png(filename ,width=3, height=4, unit='in', res = 1000)

data_sub %>%
  filter(Type == "sample") %>% 
  mutate(Time_hours = Time*24) %>%
  ggplot(aes(x = Time_hours, y = Flux, color = Condition)) +
  facet_grid(. ~ Pyruv) + 
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  geom_point(alpha = 0.3, size = 1.5) +
  geom_smooth(aes(group = Condition), span = 0.5) +
  scale_color_manual(values = colors,
                     breaks = c("pre-drought", "drought"),
                     labels = c("Pre Drought", "Drought")) + 
  labs(x = "Hours after labeling", y = "13C/(12C + 13C) flux") +
  ylim(-150, 300) +
  theme(text = element_text(size = 12,
                           family = "Arial",
                           color = "black"),
    legend.position = "bottom",
    legend.title = element_blank())
dev.off()




# Plotwise
data_sub %>%
  filter(Type == "sample") %>% 
  mutate(Time_hours = Time*24) %>%
  filter(Pyruv == "C1") %>%
  ggplot(aes(x = Time_hours, y = Flux, color = Condition)) +
  facet_wrap(. ~ Sample, ncol = 5) + 
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  geom_point(alpha = 0.3) +
  geom_smooth(aes(group =  Condition), span = 0.5) +
  scale_color_manual(values = colors,
                     breaks = c("pre-drought", "drought"),
                     labels = c("Pre Drought", "Drought")) + 
  labs(x = "Hours after labeling", y = "13C/(12C + 13C) flux", subtitle = "C1 pyurvate") +
  theme(legend.position = "bottom")
ggsave("./Figs/01_Pyruvate_acetone_C1_Plotwise.png")

data_sub %>%
  filter(Type == "sample") %>% 
  mutate(Time_hours = Time*24) %>%
  filter(Pyruv == "C2") %>%
  ggplot(aes(x = Time_hours, y = Flux, color = Condition)) +
  facet_wrap(. ~ Sample, ncol = 5) + 
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  geom_point(alpha = 0.3) +
  geom_smooth(aes(group = Condition), span = 0.5) +
  scale_color_manual(values = colors,
                     breaks = c("pre-drought", "drought"),
                     labels = c("Pre Drought", "Drought")) + 
  labs(x = "Hours after labeling", y = "13C/(12C + 13C) flux", subtitle = "C2 pyurvate") +
  theme(legend.position = "bottom")
ggsave("./Figs/01_Pyruvate_acetone_C2_Plotwise.png")

# Isotope Sig - boxplots
filename=paste("Figs/aacetone-boxplot.png", sep = "")
png(filename ,width=4, height=4, unit='in', res = 1000)

data_sub %>%
  filter(Type == "sample") %>% 
  mutate(label_state = ifelse(Time >0 & Time < 2, "post", "pre")) %>% 
  mutate(Time_hours = Time*24) %>%
  filter(label_state == "post") %>%
  ggplot(aes(x = Condition, y = Flux, fill = Condition)) +
  facet_grid(. ~ Pyruv) + 
   geom_boxplot() +
  scale_x_discrete(labels = c("Pre Drought", "Drought")) +
  scale_fill_manual(values = colors,
                     breaks = c("pre-drought", "drought"),
                     labels = c("Pre Drought", "Drought")) + 
  labs(y = "13C/(12C+13C) flux" ) +
  theme(text = element_text(size = 10,
                            family = "Arial",
                            color = "black"),
        legend.position = "bottom",
        axis.title.x = element_blank(),
        legend.title = element_blank())
dev.off()

# Isotope Sig - boxplots - C1
filename=paste("Figs/aacetone-boxplot-C1.png", sep = "")
png(filename ,width=2, height=4, unit='in', res = 1000)

data_sub %>%
  filter(Type == "sample") %>% 
  mutate(label_state = ifelse(Time >0 & Time < 2, "post", "pre")) %>% 
  mutate(Time_hours = Time*24) %>%
  filter(label_state == "post") %>%
  filter(Pyruv == "C1") %>%
  ggplot(aes(x = Condition, y = Flux, fill = Condition)) +
  facet_grid(. ~ Pyruv) + 
  geom_boxplot() +
  scale_x_discrete(labels = c("Pre Drought", "Drought")) +
  scale_fill_manual(values = colors,
                    breaks = c("pre-drought", "drought"),
                    labels = c("Pre Drought", "Drought")) + 
  ylim(-150, 300) +
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank())
dev.off()

# Isotope Sig - boxplots - C2
filename=paste("Figs/aacetone-boxplot-C2.png", sep = "")
png(filename ,width=2, height=4, unit='in', res = 1000)

data_sub %>%
  filter(Type == "sample") %>% 
  mutate(label_state = ifelse(Time >0 & Time < 2, "post", "pre")) %>% 
  mutate(Time_hours = Time*24) %>%
  filter(label_state == "post") %>%
  filter(Pyruv == "C2") %>%
  ggplot(aes(x = Condition, y = Flux, fill = Condition)) +
  facet_grid(. ~ Pyruv) + 
  geom_boxplot() +
  scale_x_discrete(labels = c("Pre Drought", "Drought")) +
  scale_fill_manual(values = colors,
                    breaks = c("pre-drought", "drought"),
                    labels = c("Pre Drought", "Drought")) + 
  ylim(-150, 300) +
  theme(text = element_text(size = 12,
                            family = "Arial",
                            color = "black"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank())
dev.off()

####statistical analysis
#Flux - Condition, in C1 and C2 separated
data_c1 <- data_sub %>%
 filter(Pyruv == "C1")

data_c2 <- data_sub %>%
  filter(Pyruv == "C2")


wilcox.test(Flux ~ Condition, data = data_c1)
#data:  Flux by Condition
#W = 98902, p-value = 0.592
#alternative hypothesis: true location shift is not equal to 0


wilcox.test(Flux ~ Condition, data = data_c2)
#data:  Flux by Condition
#W = 192527, p-value < 2.2e-16
#alternative hypothesis: true location shift is not equal to 0


