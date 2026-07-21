#================================
# Analysis C
#================================

# Load packages
library(glmmTMB)      #for glmmTMB
library(DHARMa)       #for residual diagnostics
library(MuMIn)        #for model diagnostics
library(effects)      #for partial effects plots
library(ggeffects)    #for partial effects plots
library(emmeans)      #for estimating marginal means
library(modelr)       #for auxillary modelling functions
library(patchwork)    #for multiple plots
library(ggpubr)       #for multiple plots
library(tidyverse)    #for data wrangling

#-------------------------------------
# 1. Read in data 
#-------------------------------------

# Load up file
All_regions <- read.csv('final_abundance_changes.csv', sep = ',')

All_regions <- All_regions |>
  mutate(region = as.factor(region),
         Order = as.factor(Order),
         PlaceofDevelopment = factor(PlaceofDevelopment, levels = c("Pelagic", "Demersal")),
         Habitat = factor(Habitat, levels = c("Pelagic", "Demersal")),
         accepted_name = as.factor(accepted_name)) |>
  droplevels()

#-------------------------------------
# 2 Run model and get responses
#-------------------------------------

mod.a.re <- glmmTMB(
  abundance_temp_binary ~
    position_in_range *
    #Includes the interaction with each trait but not between them:
    (logDepth + Troph + logLength + logFecundity + logK + #numerical
       Habitat + PlaceofDevelopment) + #categorical
    (1 | survey) +
    (1 | accepted_name_numeric),
  family = binomial(link = "logit"),
  data = All_regions,
  REML = TRUE
)

summary(mod.a.re)


#-------------------------------------
# 3 Define theme for plots
#-------------------------------------

# Define theme for plots
theme_fig3 <- theme_light() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = .5, linetype = 'solid'),
        axis.ticks = element_line(color = "black", linewidth=0.25, linetype='solid'),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", colour = "black", hjust = 0, vjust = 0),
        strip.placement = "inside",
        panel.spacing = unit(2, "mm"),
        legend.justification = c(1, 0), 
        legend.position = c(0.975,0.025), 
        legend.background = element_rect(fill = "transparent"),
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 7),
        legend.title = element_blank(), 
        axis.text = element_text(size = 7)
  )
  
#-------------------------------------
# 4 Plot figure Ecological trait effect size plots
#-------------------------------------

# Get trait response and contrast high and low
get_trait_response <- function(mod, dat, trait,
                               mode = c("range", "quantile")) {
  
  mode <- match.arg(mode)
  
  x <- dat[[trait]]
  grid <- list(position_in_range = c(0, 1))
  
  if (is.numeric(x)) {
    
    if (mode == "range") {
      trait_values <- seq_range(x, n = 2)
      trait_range <- c("low", "high")
    } else {
      trait_values <- quantile(x, c(.25, .75), na.rm = TRUE)
      trait_range <- c("q25", "q75")
    }
    
    grid[[trait]] <- trait_values
    
  } else if (is.factor(x)) {
    
    grid[[trait]] <- levels(x)
    trait_range <- levels(x)
    
  } else {
    stop(trait, " must be numeric or factor.")
  }
  
  emm <- mod |>
    emmeans(reformulate(c("position_in_range", trait)), at = grid) |>
    regrid()
  
  prob <- emm |>
    as_tibble() |>
    rename(TraitValue = all_of(trait)) |>
    mutate(TraitRange = rep(trait_range, each = 2))
  
  cont <- contrast(emm, interaction = "revpairwise") |>
    summary(infer = TRUE) |>
    rename(Trait_revpairwise = all_of(paste0(trait, "_revpairwise")))
  
  list(
    emm = emm,
    prob = prob,
    contrast = cont
  )
}

# Plot individual trait response
plot_trait_response <- function(mod, dat, trait, trait.name, trait_labels,
                                mode = c("range", "quantile")) {
  
  mode <- match.arg(mode)
  xgrid <- seq(0, 1, length.out = 100)
  grid <- list(position_in_range = xgrid)
  x <- dat[[trait]]
  
  if (is.numeric(x)) {
    
    if (mode == "range") {
      grid[[trait]] <- seq_range(x, n = 2)
    } else {
      grid[[trait]] <- quantile(x, c(.25, .75), na.rm = TRUE)
    }
    
    if (is.null(trait_labels[[trait]])) {
      stop("No trait labels supplied for numeric trait: ", trait)
    } else {
      traitvalue <- trait_labels[[trait]]
    }
    
  } else if (is.factor(x)) {
    
    grid[[trait]] <- levels(x)
    
    if (is.null(trait_labels[[trait]])) {
      traitvalue <- levels(x)
    } else {
      traitvalue <- trait_labels[[trait]]
    }
    
  } else {
    
    stop(trait, " must be numeric or factor.")
    
  }
  
  pred <- mod |>
    emmeans(reformulate(c("position_in_range", trait)), at = grid) |>
    regrid() |>
    as_tibble() |>
    mutate(
      traitvalue = rep(traitvalue, each = length(xgrid)),
      trait = trait,
      trait.name = trait.name
    ) |>
    dplyr::select(-all_of(trait))
  
  ggplot(pred, aes(x = position_in_range, y = prob,
                   colour = traitvalue, fill = traitvalue)) +
    geom_ribbon( aes(ymin = asymp.LCL, ymax = asymp.UCL),
                 alpha = 0.2, colour = NA) +
    geom_line(linewidth = .5) +
    facet_wrap(~trait.name) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1), breaks = c(0, .5, 1), labels = c("0", "0.5", "1")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1), breaks = c(0, .5, 1), labels = c("0", "0.5", "1")) +
    scale_colour_manual( values = c("#A40519", "#4D67F8")) +
    scale_fill_manual(values = c("#A40519", "#4D67F8")) +
    theme_fig3 +
    labs(x = NULL, y = NULL)
}

# Prepare traits and plotting criteria
traits <- c("logDepth", "Troph", "logLength", "logFecundity", "logK", "Habitat", "PlaceofDevelopment")
trait.name <- c("Depth", "Trophic level", "Body length", "Fecundity", "Growth rate", "Adult habitat", "Larval habitat")
trait.name.p1 <- c("e. Depth", "d. Trophic level", "a. Body length", "f. Fecundity", "c. Growth rate", "b. Adult habitat", "g. Larval habitat")
trait_labels <- list(
  logDepth     = c("Shallow", "Deep"),
  Troph        = c("Low", "High"), 
  logLength    = c("Small", "Large"),
  logFecundity = c("Low", "High"),
  logK         = c("Slow", "Fast"),
  Habitat      = c("Pelagic", "Demersal"),
  PlaceofDevelopment = c("Pelagic", "Demersal")
)


# Plot individual traits
plots_q <- tibble(
  trait = traits,
  trait.name = trait.name.p1
) |>
  mutate(
    plot = pmap(
      list(trait, trait.name),
      ~ plot_trait_response(
        mod = mod.a.re,
        dat = All_regions,
        trait = ..1,
        trait.name = ..2,
        trait_labels = trait_labels,
        mode = "quantile"
      )
    )
  )


# Plot effect sizes
results <- tibble(
  trait = traits,
  trait.name = trait.name,
  results = map(traits, ~get_trait_response(mod.a.re, All_regions, .x, mode = "quantile")) ) |>
  mutate(
    emm = map(results, "emm"),
    prob = map(results, "prob"),
    cont = map(results, "contrast") ) |>
  dplyr::select(-results)

plot.effects <- results |> unnest(cont) |> 
  mutate(title = "h. Poleward advantage") |> 
  mutate(trait.name = fct_reorder(trait.name, -estimate)) |> 
  ggplot(aes(y = trait.name, x = -estimate)) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "black", linewidth = 0.25) +
  geom_pointrange(aes(xmin = -asymp.LCL, xmax = -asymp.UCL), size = .2) +
  #geom_text(aes(label = round(-estimate,2), vjust = -1)) +
  #scale_colour_gradient(low = "#A40519", high = "#4D67F8") +
  facet_wrap(~title) +
  theme_fig3 +
  theme(legend.position = "none") +
  #theme(panel.grid = element_blank()) +
  labs(y = "", x = expression(Trait~effect~(ΔΔP)))

# Check results
results[1:5,] |> unnest(prob)
results[6:7,] |> unnest(prob)
results[5:7,] |> unnest(cont)

# Figure 3
left <- wrap_plots( plots_q$plot[[3]] + theme(axis.text.x = element_blank()),
                    plots_q$plot[[6]] + theme(axis.text = element_blank()), 
                    plots_q$plot[[5]] + theme(axis.text = element_blank()),
                    plot_spacer(), 
                    plots_q$plot[[2]], 
                    plots_q$plot[[1]] + theme(axis.text.y = element_blank()),
                    plots_q$plot[[4]] + theme(axis.text.y = element_blank()),
                    plots_q$plot[[7]] + theme(axis.text.y = element_blank()), 
                    ncol = 4, nrow = 2) &
  #theme(aspect.ratio = 1) &
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(.1, .25, .1, .25),"lines") )

left <- wrap_elements(left) + 
  labs(caption = expression(paste(bold("Position in range "), italic(" (0 = equatorward limit, 1 = poleward limit)"))),
       tag = expression(bold("Probability of positive response to warming"))) +
  theme(plot.caption = element_text(size = 8, color = "black", hjust = 0.5, vjust = 1),    
        plot.tag = element_text(size = 8, color = "black", angle = 90, hjust = 0.5),
        plot.tag.position = "left")

right <- plot.effects + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        plot.margin = unit(c(0,0,0,0),"lines"))

right <- wrap_elements(right) + 
  labs(caption = expression(paste(bold("                    Trait effect"), italic(" (ΔΔP)")))) +
  theme(plot.caption = element_text(size = 8, color = "black", hjust = 0.5, vjust = 1)) 

fig3 <- wrap_plots(left, right, ncol = 2) +
  plot_layout(widths = c(2.6, 1))
fig3
ggsave(filename = "figures/Fig3_v3.png", plot = fig3, width = 13*1.6, height = 9, units = "cm", dpi = 300)


#Tables
## Table 2
results |> unnest(cont) |> as_tibble() |> 
  arrange(estimate) |> 
  mutate(deltaP = -round(estimate*100,1), 
         CI = paste(-round(asymp.UCL*100,1), -round(asymp.LCL*100,1), sep = " _ "),
         z.ratio = -round(z.ratio, 2), 
         p.value = round(p.value, 3)) |> 
  dplyr::select(-c(trait, emm, prob, position_in_range_revpairwise, 
            Trait_revpairwise, estimate, SE, df, asymp.LCL, asymp.UCL)) |> 
  write_csv(file = "figures/table2_new.csv")

## Table S4
results[1:5,] |> unnest(prob) |> 
  mutate(TraitValue = as.factor(TraitValue)) |> 
  bind_rows(results[6:7,] |> unnest(prob)) |> as_tibble() |> 
  dplyr::select(-c(emm, cont, SE, df)) |> 
  mutate(prob = round(prob*100, 1), 
         asymp.LCL = round(asymp.LCL*100, 1), 
         asymp.UCL = round(asymp.UCL*100, 1)) |> 
  relocate(TraitRange, .before = TraitValue) |> 
  #filter(position_in_range == 1) #|> 
  write_csv("figures/tableS4_new.csv")


# end of script