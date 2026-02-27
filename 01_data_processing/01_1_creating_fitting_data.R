### MERGE THE RVF & ENVIRONMENTAL DATA 

# Load packages 
library(readr)
library(dplyr)

# Load dataset
env <- read_csv("01_data_processing/env_clean_by_zone.csv")
rvf <- read_csv("01_data_processing/rvf_clean_by_zone.csv")
livestock <- read_csv("01_data_processing/livestock_density_aez.csv")

## CLEAN THE RVF DATA 

# Add up cases by year and AEZ
rvf <- rvf %>%
  group_by(Year, AEZ_Name) %>%
  summarise(Number_confirmed_events = sum(Number_confirmed_events, na.rm = TRUE),
            .groups = "drop")

# Assess the data for distribution
# Subset to >1980
rvf <- rvf %>% filter(Year > 1979)

## CLEAN THE COVARIATE DATA 

# Average the covariates by year and AEZ
env <- env %>%
  select(-...1) %>%
  group_by(year, AEZ_Name) %>%
  summarise(
    across(c(tmp, pre, pdsi, pet, sm), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )
env <- env %>% filter(year > 1979) %>% rename(Year = "year")

## MERGE TO A COMPLETE DATASET 

# Merge 
dat <- left_join(env, rvf)
dat <- left_join(dat, livestock)

# Convert NAs to zeros for cases
dat$Number_confirmed_events[is.na(dat$Number_confirmed_events)] <- 0

# Save
write_csv(dat, "01_data_processing/fitting_data.csv")


