### VISUALISE AND PLOT THE OUTCOME AND COVARIATES DATA 

# Load packages 
library(readr)
library(dplyr)
library(ggplot2)
library(sf)
library(ggnewscale)
library(tidyr)
library(purrr)

# Load fitting data 
dat <- read_csv("01_data_processing/fitting_data.csv")

## VISUALISE THE OUTCOME 

# View the outcome data in space and time 
ggplot(dat, aes(x = Year, y = Number_confirmed_events, fill = AEZ_Name)) +
  geom_col() +
  facet_wrap(~ AEZ_Name) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Year",
    y = "Number of Confirmed Events"
  ) + 
  theme(text = element_text(face = "bold"))

## VISUALISE THE LIVESTOCK DATA 

# Load the data 
map <- read_sf(
  "01_data_processing/Kenya_Agro-Ecological_Zones_Data/Kenya_Agro-Ecological_Zones_Data.shp", 
  stringsAsFactors = F)
livestock <- read_csv("01_data_processing/livestock_density_aez.csv")

# Streamline the AEZ catagories 
map <- map %>%
  mutate(AEZ_Name = case_when(
    AEZ_Name == "Inner Lowland" ~ "Lowland",
    
    AEZ_Name %in% c("Lower Highland", "Upper Highland", "Tropical Alpine") ~ "Highland",
    
    AEZ_Name %in% c("Lower Midland", "Upper Midland") ~ "Midland",
    
    AEZ_Name %in% c("Waterbody - Inner Lowland",
                    "Waterbody - Lower Highland",
                    "Waterbody - Lower Midland",
                    "Waterbody - Upper Midland",
                    "Waterbody - Lake Victoria",
                    "Waterbody - Indian Ocean") ~ "Waterbody",
    
    AEZ_Name == "Coastal Lowland" ~ "Coastal",
    
    AEZ_Name == "Nairobi City" ~ "Nairobi City",
    
    TRUE ~ AEZ_Name   # safety fallback
  ))

# Merge the data 
map <- left_join(map, livestock)

# View livestock density data in space 
ggplot() +
  
  # AEZ
  geom_sf(
    data = map %>% mutate(panel = "AEZ"),
    aes(fill = AEZ_Name),
    colour = NA
  ) +
  scale_fill_brewer(
    palette = "Pastel2",
    name = "AEZ",
    guide = guide_legend(order = 1)
  ) +
  
  new_scale_fill() +
  
  # Sheep (order 4)
  geom_sf(
    data = map %>% mutate(panel = "Sheep density"),
    aes(fill = sheep_mean),
    colour = NA
  ) +
  scale_fill_distiller(
    palette = "Oranges",
    direction = 1,
    name = "Sheep density",
    guide = guide_colourbar(order = 4)
  ) +
  
  new_scale_fill() +
  
  # Goats (order 3)
  geom_sf(
    data = map %>% mutate(panel = "Goat density"),
    aes(fill = goat_mean),
    colour = NA
  ) +
  scale_fill_distiller(
    palette = "Blues",
    direction = 1,
    name = "Goat density",
    guide = guide_colourbar(order = 3)
  ) +
  
  new_scale_fill() +
  
  # Cattle (order 2)
  geom_sf(
    data = map %>% mutate(panel = "Cattle density"),
    aes(fill = cattle_mean),
    colour = NA
  ) +
  scale_fill_distiller(
    palette = "Purples",
    direction = 1,
    name = "Cattle density",
    guide = guide_colourbar(order = 2)
  ) +
  
  facet_wrap(~panel, ncol = 2) +
  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    text = element_text(face = "bold"),
    strip.text = element_text(size = 12),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.position = "right"
  )

## VISUALISE THE COVARIATE DATA 

# View the covariate data in space and time 
dat %>%
  pivot_longer(
    cols = c(tmp, pre, pdsi, pet, sm),
    names_to = "covariate",
    values_to = "value"
  ) %>%
  mutate(
    covariate = factor(
      covariate,
      levels = c("tmp", "pre", "pdsi", "pet", "sm"),
      labels = c("Temperature", "Precipitation", "Drought (PDSI)",
                 "Potential evapotranspiration (PET)", "Soil moisture")
    )
  ) %>%
  ggplot(aes(x = Year, y = value, colour = AEZ_Name, group = AEZ_Name)) +
  geom_line(linewidth = 0.8, alpha = 0.9) +
  facet_wrap(~ covariate, scales = "free_y", ncol = 2) +
  scale_colour_brewer(palette = "Set2", name = "AEZ") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    text = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    legend.position = "right"
  ) +
  labs(x = "Year", y = NULL)

# Visualise the outcome v covariates per AEZ 
plots_by_aez <- dat %>%
  mutate(AEZ_Name = as.character(AEZ_Name)) %>%
  split(.$AEZ_Name) %>%
  imap(~ .x %>%
         pivot_longer(
           cols = c(tmp, pre, pdsi, pet, sm),
           names_to = "covariate",
           values_to = "value"
         ) %>%
         mutate(
           covariate = factor(
             covariate,
             levels = c("tmp", "pre", "pdsi", "pet", "sm"),
             labels = c("Temperature", "Precipitation", "Drought (PDSI)",
                        "Potential evapotranspiration (PET)", "Soil moisture")
           )
         ) %>%
         ggplot(aes(x = Number_confirmed_events, y = value, colour = covariate)) +
         geom_point(alpha = 0.7, size = 1.6) +
         scale_colour_brewer(palette = "Pastel2") +
         facet_wrap(~ covariate, scales = "free_y", ncol = 2) +
         labs(
           title = paste("AEZ:", .y),
           x = "Number confirmed events",
           y = NULL
         ) +
         theme_minimal() +
         theme(
           panel.grid = element_blank(),
           text = element_text(face = "bold"),
           strip.text = element_text(face = "bold"),
           plot.title = element_text(face = "bold"),
           legend.position = "none"
         )
  )

# Save the plots 
dir.create("02_data_exploration/covariate_v_outcome_by_AEZ", showWarnings = FALSE)
iwalk(plots_by_aez, ~ ggsave(
  filename = file.path("covariate_scatter_by_AEZ",
                       paste0(gsub("[^A-Za-z0-9]+", "_", .y), ".png")),
  plot = .x,
  width = 12, height = 10, dpi = 300
))





