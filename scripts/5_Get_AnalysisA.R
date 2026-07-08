#================================
# Analysis A
#================================

# Load packages
library(dplyr)
library(ggplot2)
library(data.table)
library(lme4)
library(effects)
library(lmerTest)
library(ggpubr)
library(grid)         
library(maps)   

#-------------------------------------
# 1. Read in data 
#-------------------------------------

# Load up file again
All_regions <- read.csv('final_abundance_changes.csv', sep = ',') 

#-------------------------------------
# 2 Time series duration on world map
#-------------------------------------

# Get the map outline
world <- map_data("world")

# Create hex map
hex_map <- ggplot(All_regions, aes(x = lon_cell, y = lat_cell)) +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), 
               fill="grey", colour="grey70", linewidth=0.2, alpha=0.3) +
  geom_hex(bins=100, color = "grey20",
           linewidth = 0.25) +
  scale_fill_gradientn(
    colors = c("#D9E1F8", "#8094CC", "#DA3A3A"),
    trans = "log10",
    name = "Number of time series",
    guide = guide_colorbar(
      barheight = unit(5, "mm"),   
      barwidth = unit(100, "mm"),  
      direction = "horizontal",    
      title.theme = element_text(margin = margin(r = 17, b = 19))
    )
  ) +
  theme_void() +
  theme(
    legend.position = c(0.39, 0.2),  # 0.39, 0.2 Adjust these values to move the legend
    legend.justification = c(0.0, 0.0),  # Aligns the legend within the plot
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14)
  )

# Display hex map
hex_map

# Export hex map
ggsave(hex_map, file="hex_map.pdf",width = 40, height = 20, unit = "cm", dpi = 600, device = cairo_pdf)

#-------------------------------------
# 3 Taxonomic order frequency plot
#-------------------------------------

# Filter to include only orders that occur in more than 2000 time series
taxa <- All_regions %>%
  group_by(Order) %>%
  summarise(n = n()) %>%
  filter(n > 1800)

# Create order plot
orderbarplot <- ggplot(data = taxa, aes(x = reorder(Order, -n), y=n))+ 
  geom_bar(stat = "identity", width=0.5,fill='red4') +
  coord_flip() + theme_classic() +
  xlab("Taxonomic order")+ylab("Number of time series")+
  scale_y_continuous(limits = c(0,42000), expand = c(0,0)) +
  theme(axis.title = element_text(size = 20, color = 'black'),
        axis.text = element_text(size = 15, color = '#2C2323'),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = -10)))

# Display order plot
orderbarplot

# Export order plot
ggsave(orderbarplot, file="orderbarplot.pdf",width = 40, height = 20, unit = "cm", dpi = 600, device = cairo_pdf)

#-------------------------------------
# 4 - Temperature change over time
#-------------------------------------

# Load in temperature and year values for grid cells 
temp <- fread('abundances.csv') %>%
  rename(lat_cell = latitude_gridcell, lon_cell = longitude_gridcell) %>%
  select(lat_cell, lon_cell, mean_temperature, year) %>%
  mutate(grid_cell = paste(lat_cell, lon_cell, sep = '_')) %>%
  select(-lat_cell, -lon_cell) 

# Get the post-filtered grid cells 
correct_grid_cells <- fread('final_abundance_changes.csv') %>%
  select(lat_cell, lon_cell) %>%
  mutate(grid_cell = paste(lat_cell, lon_cell, sep = '_')) %>%
  select(-lat_cell, -lon_cell) %>%
  unique()

# Join the temperature data to the grid cell data
temp_data <- correct_grid_cells %>%
  left_join(temp, by = 'grid_cell') %>%
  unique()

# Obtain the estimate of temperature change over time
temp_change_model <- lmer(mean_temperature ~ year + (1|grid_cell),data=temp_data)
summary(temp_change_model)

# Convert estimate to an effect for plotting 
effect_temp_change <- as.data.frame(effect("year", temp_change_model))

# Define theme for plot
temp_plot_theme <- theme(
  axis.title.x = element_text(size = 40, vjust=1.5, hjust=0.5, colour='black'),
  axis.title.y = element_text(size = 35, colour='black'),
  axis.text = element_text(size = 32, colour='black'),
  axis.text.x = element_text(size = 32, vjust = 0.5, colour='black'),
  plot.title = element_text(lineheight=.8, face="bold", size = 10, hjust=0.5, vjust=0.2),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', linewidth=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', linewidth=0.5, linetype='solid'),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 13),
  legend.position.inside = c(0.9,0.1),
  legend.justification = c("right","bottom"))

# Create plot
temp_plot <- ggplot(data=effect_temp_change, aes(x=year,y=fit)) + 
  geom_line(colour="#C7340C") +
  geom_ribbon(aes(ymin=lower, ymax=upper), linetype=0, alpha=0.3, colour="#EBC220", fill="#EBC220") +
  scale_y_continuous(name="Mean SST (°C)", limits = c(9, 12)) +
  scale_x_continuous(name="Year") +
  theme_bw() + temp_plot_theme

# Display plot
temp_plot

# Export plot
ggsave(temp_plot, file="temp_plot.pdf",width = 30, height = 30, unit = "cm", dpi = 600, device = cairo_pdf)

# end of script