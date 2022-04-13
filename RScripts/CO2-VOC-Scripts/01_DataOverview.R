

# Visualize Pyruvate Co2 Fluxes

# Johannes Ingrisch
# 2021-08-26

# load packages
library(tidyverse)
library(plotly)
library(ggthemes)
library(patchwork)

theme_custom <- theme_few()

theme_set(theme_custom)



data <- read.csv("./Data/01_raw/S3_PyruvateLabeling_qc.csv") %>% 
  select(-X)

data <- data %>% 
  mutate(Timestamp = lubridate::ymd_hms(Timestamp_MST, tz = "America/Phoenix"), 
         Timestamp_Label_MST = lubridate::ymd_hms(Timestamp_Label_MST, tz = "America/Phoenix")) %>% 
  select(-Timestamp_MST)


# Quickview


data <- data %>% 
  select(Label, Campaign, Pyruv, Trel, Flux, D13C, Qc_iso, Qc_flux)


data_sub <- data %>% 
  filter(Qc_flux == 1, Qc_iso == 1)  %>% 
  select(-Qc_iso, -Qc_flux) %>% 
  filter(Trel > -10, Trel < 50, Flux > 0, D13C > -100)

# Isotope Sig

p1 <- data_sub %>% 
  gather(Variable, Value, -Label, -Campaign, -Pyruv, -Trel) %>% 
  filter(Variable == "D13C") %>% 
  mutate(label_state = ifelse(Trel >0, "post", "pre")) %>% 
  ggplot(aes(x = Trel, y = Value, color = Campaign)) +
  facet_grid(. ~ Pyruv) + 
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  geom_point(alpha = 0.3) +
  geom_smooth(aes(group = interaction(label_state, Campaign)), span = 0.5) +
  scale_color_manual(values = c("steelblue", "orange")) + 
  labs(x = "Hours after labeling", y = bquote(delta^13 *CO[2])) +
  theme(legend.position = "bottom")
p1

p2 <- data_sub %>% 
  gather(Variable, Value, -Label, -Campaign, -Pyruv, -Trel) %>% 
  filter(Variable == "Flux") %>% 
  mutate(label_state = ifelse(Trel >0, "post", "pre")) %>% 
  ggplot(aes(x = Trel, y = Value, color = Campaign)) +
  facet_grid(. ~ Pyruv) + 
  geom_point(alpha = 0.3) +
  geom_smooth(span = 0.5) +
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  scale_color_manual(values = c("steelblue", "orange")) + 
  labs(x = "Hours after labeling", y = bquote("Soil" ~ CO[2]~"efflux (" *mu*mol ~m^-2*s^-1*")")) +
  theme(legend.position = "bottom") +
  ylim(0,6)
p2


p1 / p2 / guide_area() + plot_layout(guides = "collect", heights = c(1,1,0.3))
ggsave("./Figs/01_Pyruvate_SoilResp_Overview.png",)  


# Plotwise
data_sub %>% 
  gather(Variable, Value, -Label, -Campaign, -Pyruv, -Trel) %>% 
  filter(Variable == "D13C", Pyruv == "C1") %>% 
  mutate(label_state = ifelse(Trel >0, "post", "pre")) %>% 
  ggplot(aes(x = Trel, y = Value, color = Campaign)) +
  facet_wrap(. ~ Label) + 
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  geom_point(alpha = 0.3) +
  geom_smooth(aes(group = interaction(label_state, Campaign)), span = 0.5) +
  scale_color_manual(values = c("steelblue", "orange")) + 
  labs(x = "Hours after labeling", y = bquote(delta^13 *CO[2]), subtitle = "C1-Pyruvate") +
  theme(legend.position = "bottom")
ggsave("./Figs/01_Pyruvate_SoilResp_C1_Plotwise.png")

data_sub %>% 
  gather(Variable, Value, -Label, -Campaign, -Pyruv, -Trel) %>% 
  filter(Variable == "D13C", Pyruv == "C2") %>% 
  mutate(label_state = ifelse(Trel >0, "post", "pre")) %>% 
  ggplot(aes(x = Trel, y = Value, color = Campaign)) +
  facet_wrap(. ~ Label) + 
  scale_x_continuous(breaks = c(-12, 0, 12, 24, 36, 48)) +
  geom_point(alpha = 0.3) +
  geom_smooth(aes(group = interaction(label_state, Campaign)), span = 0.5) +
  scale_color_manual(values = c("steelblue", "orange")) + 
  labs(x = "Hours after labeling", y = bquote(delta^13 *CO[2]), subtitle = "C2-Pyruvate") +
  theme(legend.position = "bottom")
ggsave("./Figs/01_Pyruvate_SoilResp_C2_Plotwise.png")
