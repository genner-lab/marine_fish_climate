#================================
# Analysis D
#================================

# Load packages

library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)
library(purrr)
library(maps)
library(devtools)
install_github("ltorgo/DMwR2",ref="master")
library(DMwR2)

#-------------------------------------
# 1. Read in data 
#-------------------------------------

# Load up file again
All_regions <- read.csv('final_abundance_changes.csv', sep = ',') 

#-------------------------------------
# 2 Defining theme for plots
#-------------------------------------

# Define theme for plots
effect_theme_fig4 <- theme(
  text = element_text(size = 14),
  axis.title.x = element_text(size = 41, vjust=0.3),
  axis.title.y = element_text(size = 43),
  axis.text = element_text(size = 38, colour='black'),
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
#  3 Community responses
#-------------------------------------

# # Get the community responses

Community <- All_regions %>%
  group_by(lat_cell, lon_cell, survey, region) %>%
  summarise(
    corr = cor(abundance_temp_correlation, position_in_range,
               use = "complete.obs", method = "spearman"),
    n_sp = sum(!is.na(abundance_temp_correlation) & !is.na(position_in_range)),
    .groups = "drop",
    start_year = mean (start_year),
    end_year = mean (end_year),
    duration = end_year - start_year -1
  )

Community <- Community %>% filter(duration > 9)

#-------------------------------------
#  4 Make world map
#-------------------------------------

# Get the world outline 
world <- map_data("world")

# Average corr by unique grid cell

# Plot the community responses onto the world map
hex_map_mean <- ggplot(Community, aes(x = lon_cell, y = lat_cell)) +
  geom_polygon(data = world,
               aes(x = long, y = lat, group = group),
               fill   = "grey",
               colour = "grey70",
               linewidth = 0.2,
               alpha  = 0.3) +
  stat_summary_hex(
    aes(z = corr),
    fun = mean,
    bins = 100,
    color = "grey20",
    linewidth = 0.25
  ) +
  scale_fill_gradientn(
    colours = c("#4D67F8", "#8094CC", "#D9E1F8", "#D52B2B", "#A40519"),
    limits = c(-1, 1),                  
    oob    = scales::squish,            
    breaks = c(-1, 0, 1),               
    labels = c("-1", "0", "1"),         
    name   = "Response\nmagnitude",
    guide  = guide_colorbar(
      barheight   = unit(7,  "mm"),
      barwidth    = unit(100, "mm"),
      direction   = "horizontal",
      title.theme = element_text(margin = margin(r = 13, b = 23))
    )) +
  theme_void() +
  theme(
    legend.position      = c(0.37, 0.2),
    legend.justification = c(0.0, 0.0),
    legend.text          = element_text(size = 15),
    legend.title         = element_text(size = 16),
    plot.title           = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.subtitle        = element_text(hjust = 0.5, size = 14))

# Display hex map
print(hex_map_mean)

# Export hex map
ggsave(hex_map_mean, file="community.pdf", width = 40, height = 20, units = "cm", dpi = 600, device = cairo_pdf)

#-------------------------------------
# 5 Test if the response is related to the absolute amount of temperature variation
#-------------------------------------

# Make list of start and end years, for each survey, and each cell

Survey_years <- All_regions %>%
  group_by(region, survey, lat_cell, lon_cell) %>%
  summarise(
    start_year = min(start_year, na.rm = TRUE),
    end_year = max(end_year, na.rm = TRUE),
    .groups = "drop")

Survey_years_complete <- Survey_years %>%
  mutate(year = purrr::map2(start_year, end_year, seq)) %>%
  unnest(year) %>%
  select(region, survey, lon_cell, lat_cell, year)

# Load in climate data
climate <- read.csv('HADISST.csv', head = T, sep = ',') %>%
  rename(lat_cell = lat, lon_cell = lon)

# Join climate data to abundance data
Survey_years_complete <- Survey_years_complete %>%
  left_join(climate, by = c('lat_cell', 'lon_cell', 'year'))

# select the unique grid cells and the required columns for imputation
columns_for_imputation <- Survey_years_complete %>%
  select(lon_cell, lat_cell, year, max_temperature, mean_temperature) %>%
  distinct()

# Perform imputation to fill in missing temperature values
imputed_data <- knnImputation(columns_for_imputation, k = 1) 

# Replace the old temperature values with the imputed ones
Survey_years_complete <- Survey_years_complete %>%
  select(-max_temperature, -mean_temperature) %>%
  left_join(imputed_data, by = c('lat_cell', 'lon_cell', 'year'))

write.csv(Survey_years_complete, file="Survey_years_complete.csv")

###Summarise the temperature data, and add to the community data

Survey_years_complete_summary <- Survey_years_complete %>%
  group_by(survey, lon_cell, lat_cell) %>%
  summarise(
    slope_max = coef(lm(max_temperature ~ year, na.action = na.omit))[2],
    slope_mean = coef(lm(mean_temperature ~ year, na.action = na.omit))[2],
    .groups = "drop"
  )

#Finally, add to community

Community <- Community %>%
  left_join(Survey_years_complete_summary, by = c('lat_cell', 'lon_cell', 'survey'))

#Run models to see if survey duration, temperature change or diversity predicts community response

community_response_model <- lm(corr ~ duration + slope_max + slope_mean + n_sp, data=Community)
summary(community_response_model)

community_response_model2 <- lm(corr ~ duration, data=Community)
summary(community_response_model2)

##Plot the against time...

Duration_Community <- ggplot(Community, aes(x = duration, y = corr)) +
  geom_point(aes(color = corr)) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_gradientn(
    colours = c("#4D67F8", "#8094CC", "#D9E1F8", "#D52B2B", "#A40519"),
    values  = scales::rescale(c(-1, 0, 1)),
    limits  = c(-1, 1),
    oob     = scales::squish
  ) +
  coord_cartesian(ylim = c(-1, 1), xlim = c(10, 60)) +
  scale_x_continuous(breaks = seq(10, 60, by = 10)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),      # tick labels
    axis.title = element_text(size = 16)      # axis titles
  )

Duration_Community
Duration_Community
ggsave(Duration_Community, file="Duration_Community.pdf", width = 12, height = 8, units = "cm", dpi = 600, device = cairo_pdf)


#-------------------------------------
# 7.3.2 Plotting the community response for each region
#-------------------------------------

# Create region subsets
#Oceania <- All_regions %>% filter(region == 'Australia' | region == 'New Zealand' | region == 'Indonesia')
#Europe <- All_regions %>% filter(region == 'Europe' | region == 'Spain')
#Eastern_US <- All_regions %>% filter(region == 'Eastern_US' | region == 'Belize' | region == 'Panama' | region == 'Canada' | region == 'Caribbean')
#Other <- All_regions %>% filter(region == 'Antartica' | region == 'Pacific' | region == 'American Samoa' | region == 'South_America') 
#Western_US <- All_regions %>% filter(region == 'Western_US')


# All regions
All_community <- ggplot(Community, aes(x = corr)) +
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.1, boundary = 0, closed   = "left", color    = "black") +
  geom_vline(xintercept = 0, color = "#D52B2B", linewidth  = 2, linetype   = "dashed") +
  scale_fill_gradientn(colours = c("#4D67F8", "#8094CC", "#D9E1F8", "#D52B2B", "#A40519"), values  = scales::rescale(c(-1, 0, 1)),
    limits  = c(-1, 1), oob = scales::squish, guide = "none") +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1)) +
  scale_y_continuous(breaks = c(0,200,400,600), limits = c(0, 600)) +
  labs(x = NULL, y = NULL) +
  #labs(x = "Response magnitude", y = "Number of communities") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 60), axis.title.y = element_text(size = 60), axis.text = element_text(size = 45))
All_community

# Europe

Europe_comm <- Community %>% filter(region == 'Europe' | region == 'Spain')

Europe_community <- ggplot(Europe_comm, aes(x = corr)) +
  geom_histogram(aes(fill = ..x..), binwidth = 0.1, boundary = 0, closed   = "left", color    = "black") +
  geom_vline(xintercept = 0, color = "#D52B2B", linewidth  = 2, linetype   = "dashed") +
  scale_fill_gradientn(colours = c("#4D67F8", "#8094CC", "#D9E1F8", "#D52B2B", "#A40519"), values  = scales::rescale(c(-1, 0, 1)),
                       limits  = c(-1, 1), oob = scales::squish, guide = "none") +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1)) +
  scale_y_continuous(breaks = c(0,100, 200), limits = c(0, 250)) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 45), axis.title.y = element_text(size = 45), axis.text = element_text(size = 45))
Europe_community

# Oceania

Oceania_com <- Community %>% filter(region == 'Australia' | region == 'New Zealand' | region == 'Indonesia')

Oceania_community <- ggplot(Oceania_com, aes(x = corr)) +
  geom_histogram(aes(fill = ..x..), binwidth = 0.1, boundary = 0, closed   = "left", color    = "black") +
  geom_vline(xintercept = 0, color = "#D52B2B", linewidth  = 2, linetype   = "dashed") +
  scale_fill_gradientn(colours = c("#4D67F8", "#8094CC", "#D9E1F8", "#D52B2B", "#A40519"), values  = scales::rescale(c(-1, 0, 1)),
                       limits  = c(-1, 1), oob = scales::squish, guide = "none") +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1)) +
  scale_y_continuous(breaks = c(0,25,50), limits = c(0, 60)) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 45), axis.title.y = element_text(size = 45), axis.text = element_text(size = 45))
Oceania_community

# Northwest Atlantic

Eastern_US_com <- Community %>% filter(region == 'Eastern_US' | region == 'Belize' | region == 'Panama' | region == 'Canada' | region == 'Caribbean')

Eastern_US_community <- ggplot(Eastern_US_com, aes(x = corr)) +
  geom_histogram(aes(fill = ..x..), binwidth = 0.1, boundary = 0, closed   = "left", color    = "black") +
  geom_vline(xintercept = 0, color = "#D52B2B", linewidth  = 2, linetype   = "dashed") +
  scale_fill_gradientn(colours = c("#4D67F8", "#8094CC", "#D9E1F8", "#D52B2B", "#A40519"), values  = scales::rescale(c(-1, 0, 1)),
                       limits  = c(-1, 1), oob = scales::squish, guide = "none") +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1)) +
  scale_y_continuous(breaks = c(0,100, 200), limits = c(0, 250)) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 45), axis.title.y = element_text(size = 45), axis.text = element_text(size = 45))
Eastern_US_community

Western_US_com <- Community %>% filter(region == 'Western_US')

# Northwest Pacific
Western_US_community <- ggplot(Western_US_com, aes(x = corr)) +
  geom_histogram(aes(fill = ..x..), binwidth = 0.1, boundary = 0, closed   = "left", color    = "black") +
  geom_vline(xintercept = 0, color = "#D52B2B", linewidth  = 2, linetype   = "dashed") +
  scale_fill_gradientn(colours = c("#4D67F8", "#8094CC", "#D9E1F8", "#D52B2B", "#A40519"), values  = scales::rescale(c(-1, 0, 1)),
                       limits  = c(-1, 1), oob = scales::squish, guide = "none") +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1, 1)) +
  scale_y_continuous(breaks = c(0,75, 150), limits = c(0, 150)) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 45), axis.title.y = element_text(size = 45), axis.text = element_text(size = 45))
Western_US_community

# Export plots
ggsave(All_community, file="All_community.pdf", width = 35, height = 30, units = "cm", dpi = 600, device = cairo_pdf)
ggsave(Europe_community, file="Europe_community.pdf", width = 35, height = 30, units = "cm", dpi = 600, device = cairo_pdf)
ggsave(Oceania_community, file="Oceania_community.pdf", width = 35, height = 30, units = "cm", dpi = 600, device = cairo_pdf)
ggsave(Eastern_US_community, file="Eastern_US_community.pdf", width = 35, height = 30, units = "cm", dpi = 600, device = cairo_pdf)
ggsave(Western_US_community, file="Western_US_community.pdf", width = 35, height = 30, units = "cm", dpi = 600, device = cairo_pdf)

#some statistics

binom.test(sum(Community$corr > 0, na.rm = TRUE),
           nrow(Community),
           p = 0.5)

binom.test(sum(Western_US_com$corr > 0, na.rm = TRUE),
           nrow(Western_US_com),
           p = 0.5)

binom.test(sum(Eastern_US_com$corr > 0, na.rm = TRUE),
           nrow(Eastern_US_com),
           p = 0.5)

binom.test(sum(Europe_comm$corr > 0, na.rm = TRUE),
           nrow(Europe_comm),
           p = 0.5)

binom.test(sum(Oceania_com$corr > 0, na.rm = TRUE),
           nrow(Oceania_com),
           p = 0.5)

#Example plot, for vcell 61.0, -177.5 EBS

ExamplePlotData <- All_regions %>%
  filter(survey == "EBS", lat_cell == "58", lon_cell == "-165")

ggplot(ExamplePlotData, aes(x = position_in_range, y = abundance_temp_correlation)) +
  geom_point(size = 3) +
  theme_classic() + geom_smooth(method = "lm",  color = "black", fill = "grey80") +
  coord_cartesian(ylim = c(-0.4, 0.4), xlim = c(0.3, 1)) +
  scale_x_continuous(breaks = seq(0.3, 1, by = 0.1)) +
  labs(
    x = "position_in_range",
    y = "abundance_temp_correlation") +
      theme(
        axis.text = element_text(size = 14),      # tick labels
        axis.title = element_text(size = 16)      # axis titles
  )

#End of Code