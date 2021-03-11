# This Script -------------------------------------------------------------
# Author: Laura Zwyssig
# Goal: Create uninsured and cardiovascular health shock indicators
# Last edited: 11-03-2021 Pia Arce

# Preliminaries -----------------------------------------------------------
xlibrary <- c("tidyverse", "readr")

# load libraries and automatically install all packages if missing
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = xlibrary)

options(xtable.floating = FALSE)
options(xtable.timestamp = "")

# Clear workspace
rm(list=ls())

# Load data ---------------------------------------------------------------
# Set path to main folder
path <- "C:/Users/parce/Desktop/GeiGhei/GxInsurance/"
hrs <- readRDS(file.path(path,"1_data/interim/hrs2.rds"))

# NOT done in this version: dropping observations which were created
# during reshape, but actually don't represent actual observed rounds
# Still having them makes working with lags easier in the following
# These observations will fall out anyway later when we impose the
# restriction that the outcome cannot be missing
# Change data format --------------------------------------------------
hrs <- zap_formats(hrs)

# Create uninsured dummy --------------------------------------------------
hrs <- hrs %>%
  mutate(uninsured = ifelse(
    covr == 0 &
      covs == 0 &
      hiothp == 0 &
      higov == 0,
    1, 0
  ))

# Create CV dummy ---------------------------------------------------------
# Only leads to a health shock if other illness in group didn't occur in a previous wave
hrs <- hrs %>% 
  group_by(hhidpn) %>% 
  mutate(lag_hearte = dplyr::lag(hearte, n = 1, default = NA),
         lag_stroke = dplyr::lag(stroke, n = 1, default = NA),
         lag_hearte = ifelse(
           hearte == 0, 0, lag_hearte),
         lag_stroke = ifelse(
           stroke == 0, 0, lag_stroke)
         ) %>% 
  mutate(cv = ifelse(
    lag_hearte == 0 &
    lag_stroke == 0 &
      (
        hearte == 1 |
          stroke == 1
      ),
    1, 0)) %>% 
  ungroup()



# Create mortality dummy --------------------------------------------------
hrs <- hrs %>% 
  mutate(alive = ifelse(iwstat == 5 | iwstat == 6 & iwstat!=0, 0, 1))


# Create a mean over cognitive test scores for each individual
cog_ind <- hrs[c("hhidpn","cogtot")]
mean_cog <- aggregate(cog_ind[2],by = list(cog_ind$hhidpn), mean, na.rm = TRUE)
colnames(mean_cog) <- c("hhidpn", "cog_mean")
hrs <- left_join(hrs,mean_cog, by = c("hhidpn"))


# Create a mean over risk aversion scores for each individual
risk_ind <- hrs[c("hhidpn","risk")]
mean_risk <- aggregate(risk_ind[2],by = list(risk_ind$hhidpn),mean, na.rm = TRUE)
colnames(mean_risk) <- c("hhidpn", "risk_mean")
hrs <- left_join(hrs,mean_risk, by = c("hhidpn"))


# Remove the duplicated individuals
hrs$d <- hrs$hhidpn + hrs$year + hrs$round #548145
dupe <- hrs[duplicated(hrs$d),] # 9492 dupe
hrs <- hrs[!duplicated(hrs$d), ] #538653 

# Save data ---------------------------------------------------------------
 saveRDS(hrs,file.path(path, "1_data/processed/hrs_final2.rds"))

