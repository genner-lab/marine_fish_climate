#================================
# Analysis B
#================================

# Load packages
library(dplyr)
library(ggplot2)
library(lme4)
library(effects)
library(ggpubr)

#-------------------------------------
# 1. Read in data 
#-------------------------------------

# Load up file again
All_regions <- read.csv('final_abundance_changes.csv', sep = ',') 

# Create region subsets
Oceania <- All_regions %>% filter(region == 'Australia' | region == 'New Zealand' | region == 'Indonesia')
Eastern_US <- All_regions %>% filter(region == 'Eastern_US' | region == 'Belize' | region == 'Panama' | region == 'Canada' | region == 'Caribbean')
Europe <- All_regions %>% filter(region == 'Europe' | region == 'Spain')
Western_US <- All_regions %>% filter(region == 'Western_US')
Other <- All_regions %>% filter(region == 'Antartica' | region == 'Pacific' | region == 'American Samoa' | region == 'South_America') 

#-------------------------------------
# 2 Time series of 10+ year duration
#-------------------------------------

# Filter for 10+ year duration
All_ten <- All_regions %>% filter(duration > 9)

# Create region subsets
Oceania_ten <- Oceania %>% filter(duration > 9)
Eastern_US_ten <- Eastern_US %>% filter(duration > 9)
Europe_ten <- Europe %>% filter(duration > 9)
Western_US_ten <- Western_US %>% filter(duration > 9)

#-------------------------------------
# 3 Time series of 20+ year duration
#-------------------------------------

# Filter for 20+ year duration
All_twenty <- All_regions %>% filter(duration > 19)

# Create region subsets
Oceania_twenty <- Oceania %>% filter(duration > 19)
Eastern_US_twenty <- Eastern_US %>% filter(duration > 19)
Europe_twenty <- Europe %>% filter(duration > 19)
Western_US_twenty <- Western_US %>% filter(duration > 19)

#-------------------------------------
# 4 Time series of 30+ year duration
#-------------------------------------

# Filter for 30+ year duration
All_thirty <- All_regions %>% filter(duration > 29)

# Create region subsets
Oceania_thirty <- Oceania %>% filter(duration > 29)
Eastern_US_thirty <- Eastern_US %>% filter(duration > 29)
Europe_thirty <- Europe %>% filter(duration > 29)
Western_US_thirty <- Western_US %>% filter(duration > 29)

#-------------------------------------
# 5 Density plot_All_ten
#-------------------------------------

# ---- Parameters ----
n_bins <- 300   # thin bars
bar_width <- 0.005
y_scale <- 100   # multiply diff to show percentage

# ---- Compute difference per bin ----
df_bin <- All_ten %>%
  mutate(bin = cut(position_in_range, breaks = n_bins)) %>%
  group_by(bin) %>%
  summarise(
    position = mean(position_in_range),
    diff = mean(2*abundance_temp_binary - 1),
    .groups = "drop"
  ) %>%
  mutate(fill_color = ifelse(diff >= 0, "Positive", "Negative"))

# ---- Top plot: 

top_plot <- ggplot() +
  geom_col(data = df_bin, aes(x = position, y = diff * y_scale, fill = fill_color),
           width = bar_width) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  scale_fill_manual(values = c("Positive" = "red3", "Negative" = "#3B7EA1")) +
  scale_y_continuous(
    limits = c(-40, 40),
    breaks = seq(-40, 40, by = 10)
  ) +
  labs(x = NULL, y = "% positive − % negative") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.4),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    panel.grid = element_blank()
  )

# ---- Bottom plot: 
binwidth_bottom <- (1/300)
bottom_plot <- ggplot(All_ten, aes(x = position_in_range)) +
  geom_histogram(binwidth = binwidth_bottom, fill = "black", colour = "black", linewidth = 0.1) +
  labs(x = "Position in range", y = "Count") +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.4),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    panel.grid = element_blank()
  )

# ---- Combine plots All_regions--- ----
All_ten_plot <- ggpubr::ggarrange(
  top_plot,
  bottom_plot,
  ncol = 1,
  heights = c(3, 1),   # top plot taller
  align = "v"
)
All_ten_plot

#-------------------------------------
# 6 Density plot NE Pacific
#-------------------------------------

# ---- Parameters ----
n_bins <- 300   # thin bars
bar_width <- 0.005
y_scale <- 100   # multiply diff to show percentage

# ---- Compute difference per bin ----
df_bin <- Western_US_ten %>%
  mutate(bin = cut(position_in_range, breaks = n_bins)) %>%
  group_by(bin) %>%
  summarise(
    position = mean(position_in_range),
    diff = mean(2*abundance_temp_binary - 1),
    .groups = "drop"
  ) %>%
  mutate(fill_color = ifelse(diff >= 0, "Positive", "Negative"))

# ---- Top plot

top_plot <- ggplot() +
  geom_col(data = df_bin, aes(x = position, y = diff * y_scale, fill = fill_color),
           width = bar_width) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  scale_fill_manual(values = c("Positive" = "red3", "Negative" = "#3B7EA1")) +
  scale_y_continuous(
    limits = c(-80, 80),
    breaks = seq(-80, 80, by = 20)
  ) +
  labs(x = NULL, y = "% positive − % negative") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.4),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    panel.grid = element_blank()
  )

# ---- Bottom plot: 

binwidth_bottom <- (1/300)
bottom_plot <- ggplot(Western_US_ten, aes(x = position_in_range)) +
  geom_histogram(binwidth = binwidth_bottom, fill = "black", colour = "black", linewidth = 0.1) +
  labs(x = "Position in range", y = "Count") +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.4),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    panel.grid = element_blank()
  )

# ---- Combine plots All_regions--- ----
NEPacificplot <- ggpubr::ggarrange(
  top_plot,
  bottom_plot,
  ncol = 1,
  heights = c(3, 1),   # top plot taller
  align = "v"
)
NEPacificplot

#-------------------------------------
# 7 Density plot NW Atlantic
#-------------------------------------

# ---- Parameters ----
n_bins <- 300   # thin bars
bar_width <- 0.005
y_scale <- 100   # multiply diff to show percentage

# ---- Compute difference per bin ----
df_bin <- Eastern_US_ten %>%
  mutate(bin = cut(position_in_range, breaks = n_bins)) %>%
  group_by(bin) %>%
  summarise(
    position = mean(position_in_range),
    diff = mean(2*abundance_temp_binary - 1),
    .groups = "drop"
  ) %>%
  mutate(fill_color = ifelse(diff >= 0, "Positive", "Negative"))

# ---- Top plot

top_plot <- ggplot() +
  geom_col(data = df_bin, aes(x = position, y = diff * y_scale, fill = fill_color),
           width = bar_width) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  scale_fill_manual(values = c("Positive" = "red3", "Negative" = "#3B7EA1")) +
  scale_y_continuous(
    limits = c(-60, 60),
    breaks = seq(-60, 60, by = 30)
  ) +
  labs(x = NULL, y = "% positive − % negative") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.4),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    panel.grid = element_blank()
  )

# ---- Bottom plot: 

binwidth_bottom <- (1/300)
bottom_plot <- ggplot(Eastern_US_ten, aes(x = position_in_range)) +
  geom_histogram(binwidth = binwidth_bottom, fill = "black", colour = "black", linewidth = 0.1) +
  labs(x = "Position in range", y = "Count") +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.4),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    panel.grid = element_blank()
  )

# ---- Combine plots All_regions--- ----
NWAtlanticplot <- ggpubr::ggarrange(
  top_plot,
  bottom_plot,
  ncol = 1,
  heights = c(3, 1),   # top plot taller
  align = "v"
)
NWAtlanticplot


#-------------------------------------
# 8 Density plot Europe
#-------------------------------------

# ---- Parameters ----
n_bins <- 300   # thin bars
bar_width <- 0.005
y_scale <- 100   # multiply diff to show percentage

# ---- Compute difference per bin ----
df_bin <- Europe_ten %>%
  mutate(bin = cut(position_in_range, breaks = n_bins)) %>%
  group_by(bin) %>%
  summarise(
    position = mean(position_in_range),
    diff = mean(2*abundance_temp_binary - 1),
    .groups = "drop"
  ) %>%
  mutate(fill_color = ifelse(diff >= 0, "Positive", "Negative"))

# ---- Top plot

top_plot <- ggplot() +
  geom_col(data = df_bin, aes(x = position, y = diff * y_scale, fill = fill_color),
           width = bar_width) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  scale_fill_manual(values = c("Positive" = "red3", "Negative" = "#3B7EA1")) +
  scale_y_continuous(
    limits = c(-80, 80),
    breaks = seq(-80, 80, by = 20)
  ) +
  labs(x = NULL, y = "% positive − % negative") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.4),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    panel.grid = element_blank()
  )

# ---- Bottom plot: 

binwidth_bottom <- (1/300)
bottom_plot <- ggplot(Europe_ten, aes(x = position_in_range)) +
  geom_histogram(binwidth = binwidth_bottom, fill = "black", colour = "black", linewidth = 0.1) +
  labs(x = "Position in range", y = "Count") +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.4),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    panel.grid = element_blank()
  )

# ---- Combine plots All_regions--- ----
NEAtlanticplot <- ggpubr::ggarrange(
  top_plot,
  bottom_plot,
  ncol = 1,
  heights = c(3, 1),   # top plot taller
  align = "v"
)
NEAtlanticplot

#-------------------------------------
# 9 Density plot Oceania
#-------------------------------------

# ---- Parameters ----
n_bins <- 300   # thin bars
bar_width <- 0.005
y_scale <- 100   # multiply diff to show percentage

# ---- Compute difference per bin ----
df_bin <- Oceania_ten %>%
  mutate(bin = cut(position_in_range, breaks = n_bins)) %>%
  group_by(bin) %>%
  summarise(
    position = mean(position_in_range),
    diff = mean(2*abundance_temp_binary - 1),
    .groups = "drop"
  ) %>%
  mutate(fill_color = ifelse(diff >= 0, "Positive", "Negative"))

# ---- Top plot

top_plot <- ggplot() +
  geom_col(data = df_bin, aes(x = position, y = diff * y_scale, fill = fill_color),
           width = bar_width) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  scale_fill_manual(values = c("Positive" = "red3", "Negative" = "#3B7EA1")) +
  scale_y_continuous(
    limits = c(-100, 100),
    breaks = seq(-100, 100, by = 50)
  ) +
  labs(x = NULL, y = "% positive − % negative") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.4),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    panel.grid = element_blank()
  )

# ---- Bottom plot: 

binwidth_bottom <- (1/300)
bottom_plot <- ggplot(Oceania_ten, aes(x = position_in_range)) +
  geom_histogram(binwidth = binwidth_bottom, fill = "black", colour = "black", linewidth = 0.1) +
  labs(x = "Position in range", y = "Count") +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.4),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    panel.grid = element_blank()
  )

# ---- Combine plots All_regions--- ----
Oceanaplot <- ggpubr::ggarrange(
  top_plot,
  bottom_plot,
  ncol = 1,
  heights = c(3, 1),   # top plot taller
  align = "v"
)
Oceanaplot

#--------------------------------
# 10 Models plot of effects of time series duration
#--------------------------------

model_all_ten <- glmer(abundance_temp_binary ~ position_in_range + (1|accepted_name), family = binomial("logit"), data = All_ten)
effect_model_all_ten <- Effect("position_in_range", model_all_ten)
effect_model_all_ten <- as.data.frame(effect_model_all_ten)
summary(model_all_ten)
        
model_all_twenty <- glmer(abundance_temp_binary ~ position_in_range + (1|accepted_name), family = binomial("logit"), data = All_twenty)
effect_model_all_twenty <- Effect("position_in_range",model_all_twenty)
effect_model_all_twenty <- as.data.frame(effect_model_all_twenty)
summary(model_all_twenty)

model_all_thirty <- glmer(abundance_temp_binary ~ position_in_range + (1|accepted_name), family = binomial("logit"), data = All_thirty)
effect_model_all_thirty <- Effect("position_in_range",model_all_thirty)
effect_model_all_thirty <- as.data.frame(effect_model_all_thirty)
summary(model_all_thirty)

effect_theme <- theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 10, vjust=0.1),
  axis.title.y = element_text(size = 10),
  axis.text = element_text(size = 8, colour='black'),
  axis.text.x = element_text(vjust = 0),
  plot.title = element_text(lineheight=1, size = 11, hjust=0.5, vjust=0.000001),
  panel.background = element_rect(fill = "white",colour = "white"),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', linewidth=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', linewidth=0.5, linetype='solid'),
  legend.justification = c("right", "bottom"),legend.position = c(.90,.10),legend.text = element_text(size = 10))

cols <- c("10+ Years" = "#A40519","20+ Years" = "black","30+ Years" = "#4D67F8")

plotmodel_all_times <- ggplot(
  data = effect_model_all_ten,
  aes(x = position_in_range, y = fit)) + 
  scale_y_continuous(
    name = "Abundance Change",
    breaks = c(0.35, 0.5, 0.65),
    limits = c(0.35, 0.65)) +
  scale_x_continuous(
    name = "Position in Range",
    breaks = c(0.0, 0.25, 0.5, 0.75, 1.0),
    limits = c(0, 1)) +
  geom_line(
    data = effect_model_all_ten,
    aes(colour = "10+ Years"),
    linewidth = 1) +
  geom_ribbon(
    data = effect_model_all_ten,
    aes(ymin = lower, ymax = upper, fill = "10+ Years"),
    alpha = 0.3,
    colour = NA) +
  geom_line(
    data = effect_model_all_twenty,
    aes(colour = "20+ Years"),
    linewidth = 1) +
  geom_ribbon(
    data = effect_model_all_twenty,
    aes(ymin = lower, ymax = upper, fill = "20+ Years"),
    alpha = 0.3,
    colour = NA) +
  geom_line(
    data = effect_model_all_thirty,
    aes(colour = "30+ Years"),
    linewidth = 1) +
  geom_ribbon(
    data = effect_model_all_thirty,
    aes(ymin = lower, ymax = upper, fill = "30+ Years"),
    alpha = 0.3,
    colour = NA) +
  scale_color_manual(
    values = cols,
    name = "Length of timeseries") +
  scale_fill_manual(
    values = cols,
    guide = "none") +
  theme(axis.title.x = element_text(hjust = -0.1)) +
  effect_theme

plotmodel_all_times

#--------------------------------
# 10 Regional Models of full time series duration
#--------------------------------

model_WesternUS_ten <- glmer(abundance_temp_binary ~ position_in_range + (1|accepted_name), family = binomial("logit"), data = Western_US_ten)
summary(model_WesternUS_ten)
model_WesternUS_twenty <- glmer(abundance_temp_binary ~ position_in_range + (1|accepted_name), family = binomial("logit"), data = Western_US_twenty)
summary(model_WesternUS_twenty)
model_WesternUS_thirty <- glmer(abundance_temp_binary ~ position_in_range + (1|accepted_name), family = binomial("logit"), data = Western_US_thirty)
summary(model_WesternUS_thirty)

model_EasternUS_ten <- glmer(abundance_temp_binary ~ position_in_range + (1|accepted_name), family = binomial("logit"), data = Eastern_US_ten)
summary(model_EasternUS_ten)
model_EasternUS_twenty <- glmer(abundance_temp_binary ~ position_in_range + (1|accepted_name), family = binomial("logit"), data = Eastern_US_twenty)
summary(model_EasternUS_twenty)
model_EasternUS_thirty <- glmer(abundance_temp_binary ~ position_in_range + (1|accepted_name), family = binomial("logit"), data = Eastern_US_thirty)
summary(model_EasternUS_thirty)

model_Europe_ten <- glmer(abundance_temp_binary ~ position_in_range + (1|accepted_name), family = binomial("logit"), data = Europe_ten)
summary(model_Europe_ten)
model_Europe_twenty <- glmer(abundance_temp_binary ~ position_in_range + (1|accepted_name), family = binomial("logit"), data = Europe_twenty)
summary(model_Europe_twenty)
model_Europe_thirty <- glmer(abundance_temp_binary ~ position_in_range + (1|accepted_name), family = binomial("logit"), data = Europe_thirty)
summary(model_Europe_thirty)

model_Oceania_ten <- glmer(abundance_temp_binary ~ position_in_range + (1|accepted_name), family = binomial("logit"), data = Oceania_ten)
summary(model_Oceania_ten)
model_Oceania_twenty <- glmer(abundance_temp_binary ~ position_in_range + (1|accepted_name), family = binomial("logit"), data = Oceania_twenty)
summary(model_Oceania_twenty)
model_Oceania_thirty <- glmer(abundance_temp_binary ~ position_in_range + (1|accepted_name), family = binomial("logit"), data = Oceania_thirty)
summary(model_Oceania_thirty)

#--------------------------------
# 11 Core models, using correlation instead of binary
#--------------------------------

model_All_ten_corr <- lmer(abundance_temp_correlation ~ position_in_range + (1|accepted_name), data = All_ten)
summary(model_All_ten_corr)

model_WesternUS_ten_corr <- lmer(abundance_temp_correlation ~ position_in_range + (1|accepted_name), data = Western_US_ten)
summary(model_WesternUS_ten_corr)

model_EasternUS_ten_corr <- lmer(abundance_temp_correlation ~ position_in_range + (1|accepted_name), data = Eastern_US_ten)
summary(model_EasternUS_ten_corr)

model_Europe_ten_corr <- lmer(abundance_temp_correlation ~ position_in_range + (1|accepted_name), data = Europe_ten)
summary(model_Europe_ten_corr)

model_Oceania_ten_corr <- lmer(abundance_temp_correlation ~ position_in_range + (1|accepted_name), data = Oceania_ten)
summary(model_Oceania_ten_corr)

# end of script