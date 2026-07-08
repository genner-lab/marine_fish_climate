#================================
# Analysis C
#================================

# Load packages
library(ggplot2)
library(lme4)
library(scales)
library(sjPlot)
library(ggstatsplot)  

#-------------------------------------
# 1. Read in data 
#-------------------------------------

# Load up file again
All_regions <- read.csv('final_abundance_changes.csv', sep = ',') 

# Make habitat and larval development factors
All_regions$Habitat <- as.factor(All_regions$Habitat)
All_regions$PlaceofDevelopment <- as.factor(All_regions$PlaceofDevelopment)

#-------------------------------------
# 2 Define theme for plots
#-------------------------------------

# Define theme for plots
effect_theme_fig3 <- theme(
  text = element_text(size = 14),
  axis.title.x = element_text(size = 23, vjust=0.1),
  axis.title.y = element_text(size = 22),
  axis.text = element_text(size = 25, colour='black'),
  axis.text.x = element_text(vjust = 0),
  plot.title = element_text(lineheight=1, size = 19, hjust=0.5, vjust=0.000001),
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', linewidth=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', linewidth=0.5, linetype='solid'),
  legend.justification = c(1, 0), legend.position = c(0.9,0.1),legend.text = element_text(size = 15))

#-------------------------------------
# 3 Run models
#-------------------------------------

# Habitat model                          
habitat <- glmer(abundance_temp_binary ~ position_in_range * Habitat + (1 | accepted_name_numeric), family = binomial("logit"), data = All_regions)
summary(habitat)

# Larval habitat model
larval_habitat <- glmer(abundance_temp_binary ~ position_in_range * PlaceofDevelopment + (1 | accepted_name_numeric), family = binomial("logit"), data = All_regions)
summary(larval_habitat)

# Depth model
depth <- glmer(abundance_temp_binary ~ position_in_range * logDepth + (1 | accepted_name_numeric), family = binomial("logit"), data = All_regions)
summary(depth)

# Trophic level model
trophic <- glmer(abundance_temp_binary ~ position_in_range * Troph + (1|accepted_name_numeric), family = binomial("logit"), data = All_regions)
summary(trophic)

# Body_length model
body_length <- glmer(abundance_temp_binary ~ position_in_range * logLength + (1 | accepted_name_numeric), family = binomial("logit"), data = All_regions)
summary(body_length)

# Fecundity model
fecundity <- glmer(abundance_temp_binary ~ position_in_range * logFecundity + (1 | accepted_name_numeric), family = binomial("logit"), data = All_regions)
summary(fecundity)

# Growth_rate model
growth_rate <- glmer(abundance_temp_binary ~ position_in_range * logK + (1 | accepted_name_numeric), family = binomial("logit"), data = All_regions)
summary(growth_rate)


#-------------------------------------
# 4 Fig 3 - Ecological trait effect size plots
#-------------------------------------

# Create plot for different habitats
habitat_plot <- plot_model(habitat, type = "pred", terms = c("position_in_range [all]", "Habitat"), title = '') + ylab("") + xlab("") + scale_fill_manual(name = "Habitat type", values = c("#0072B2", "darkgrey", "chocolate", "purple")) +
  theme_classic() + scale_color_manual(name = "Habitat type", values = c("#0072B2", "black")) +
  scale_y_continuous(labels = number_format(accuracy = 0.01)) + effect_theme_fig3 +
  theme(legend.justification = c("right", "bottom"), legend.position = c(0.99, .03), legend.title = element_text(size = 23), legend.text = element_text(size = 23))
habitat_plot

# Create plot for different larval phase habitats
larval_plot <- plot_model(larval_habitat, type = "int", terms = c("position_in_range [all]", "PlaceofDevelopment"), title = '') + ylab("") + xlab("") + scale_fill_manual(name = "Habitat type", values = c("#0072B2", "darkgrey", "chocolate", "purple")) +
  theme_classic() + scale_color_manual(name = "Larval habitat type", values = c("#0072B2", "darkgrey")) +
  scale_y_continuous(labels = number_format(accuracy = 0.01)) + effect_theme_fig3 +
  theme(legend.justification = c("right", "bottom"), legend.position = c(.96, .08), legend.title = element_text(size = 23), legend.text = element_text(size = 24)) 
larval_plot

# Create plot for different depths
depth_plot <- plot_model(depth, mdrt.values = "quart", type = "pred", terms = c("position_in_range [all]", "logDepth"), title = '') +
  ylab("") + xlab("") + scale_fill_manual(values = c("#0072B2", "darkgrey", "chocolate")) + 
  scale_color_manual(values = c("#0072B2", "darkgrey", "chocolate"), labels = c("Low", "Medium", "High")) +
  labs(color = "Depth", fill  = "Depth") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +  effect_theme_fig3 +
  theme(legend.justification = c("right", "bottom"), legend.position = c(.9, .1), legend.title = element_text(size = 23), legend.text = element_text(size = 24))
depth_plot

# Create plot for different trophic levels
trophic_plot <- plot_model(trophic, mdrt.values = "quart", type = "int", terms = c("position_in_range [all]", "trophic"), title = '') +
  ylab("") + xlab("") + scale_fill_manual(values = c("#0072B2", "darkgrey", "chocolate")) + 
  scale_color_manual(values = c("#0072B2", "darkgrey", "chocolate"), labels = c("Low", "Medium", "High")) +
  labs(color = "Trophic level", fill  = "Trophic level") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +  effect_theme_fig3 +
  theme(legend.justification = c("right", "bottom"), legend.position = c(.9, .1), legend.title = element_text(size = 23), legend.text = element_text(size = 24))
trophic_plot

# Create plot for different body lengths
body_length_plot <- plot_model(body_length, mdrt.values = "quart", type = "int", terms = c("position_in_range [all]", "Body length"), title = '') +
  ylab("") + xlab("") + scale_fill_manual(values = c("#0072B2", "darkgrey", "chocolate")) + 
  scale_color_manual(values = c("#0072B2", "darkgrey", "chocolate"), labels = c("Short", "Medium", "Long")) +
  labs(color = "Body length", fill  = "Body length") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + effect_theme_fig3 + 
  theme(legend.justification = c("right", "bottom"), legend.position = c(.865, .1), legend.title = element_text(size = 24), legend.text = element_text(size = 25))
body_length_plot

# Create plot for different fecundities
fecundity_plot <- plot_model(fecundity, mdrt.values = "quart", type = "int", terms = c("position_in_range [all]", "Fecundity"), title = '') +
  ylab("") + xlab("") + scale_fill_manual(values = c("#0072B2", "darkgrey", "chocolate")) + 
  scale_color_manual(values = c("#0072B2", "darkgrey", "chocolate"), labels = c("Low", "Medium", "High")) +
  labs(color = "Fecundity", fill  = "Fecundity") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + effect_theme_fig3 +
  theme(legend.justification = c("right", "bottom"), legend.position = c(.84, .1), legend.title = element_text(size = 24), legend.text = element_text(size = 25))
fecundity_plot

# Create plot for different growth rates
growth_rate_plot <- plot_model(growth_rate, mdrt.values = "quart", type = "int", terms = c("position_in_range [all]", "growth rate"), title = '') +
  ylab("") + xlab("") + scale_fill_manual(values = c("#0072B2", "darkgrey", "chocolate")) + 
  scale_color_manual(values = c("#0072B2", "darkgrey", "chocolate"), labels = c("Slow", "Medium", "Fast")) +
  labs(color = "Growth Rates", fill  = "Growth Rates") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) + effect_theme_fig3 +
  theme(legend.justification = c("right", "bottom"), legend.position = c(.9, .1), legend.title = element_text(size = 24), legend.text = element_text(size = 25))
growth_rate_plot

# correlogram

life_history <- read.csv('life_history_31March2026.csv', sep = ',')

unique_accepted_names <- as.data.frame(unique(All_regions$accepted_name))
names(unique_accepted_names)[1] <- "accepted_name"

life_history_filtered <- life_history %>%
  filter(accepted_name %in% unique_accepted_names$accepted_name)

colnames(life_history_filtered)[2] <- c("Maximum length")
colnames(life_history_filtered)[4] <- c("Maximum depth")
colnames(life_history_filtered)[8] <- c("Trophic level")
colnames(life_history_filtered)[9] <- c("Growth rate")

ggstatsplot::ggcorrmat(
  data = life_history_filtered[,c("Maximum length", "Maximum depth", "Fecundity", "Trophic level", "Growth rate")],
  type = "nonparametric",
  colors = c("darkred", "white", "steelblue") # change default colors
)

# end of script