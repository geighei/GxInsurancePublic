# This Script -------------------------------------------------------------
# Author: Laura Zwyssig
# Goal: Reshape RAND data file to long format for relevant variables
# Last edited: 20.06.2018

# Preliminaries -----------------------------------------------------------
library(plyr)
library(tidyverse)
library(readr)
library(haven)
library(stringr)
library(reshape2)

# clear workspace
rm(list=ls())

# Load data ---------------------------------------------------------------
#path <- "~/Dropbox/GxInsurance"
path <- "C:/Users/parce/Desktop/GeiGhei/GxInsurance"

df <- readRDS(file.path(path,"1_data/interim/randhrs1992_2016v2.rds"))

# Prepare reshape ---------------------------------------------------------
# get base names for all time-varying variables
vars_varying_base <- c("iwendm", "iwendy", "agem_e", "agey_e", "iwstat", 
                       "wthh", "wtresp",
                       "drink",
                       "smokev", "smoken", "drink", "drinkr", "drinkd", "drinkn",
                       "shlt", "vigact", "vgactx", "mdactx", "ltactx",
                       "cancr", "lung", "heart", "hibp", "diab", "strok", "psych", "arthr",
                       "higov", "govmr", "govmd", "govva",
                       "covr", "covs",
                       "henum",
                       "hiothp",
                       "oopmd",
                       "iearn",
                       "lbrf",
                       "mstat", "mpart", "risk", "cogtot"
)

# Select all time-varying variables
vars_varying <- c()
for (v in vars_varying_base) {
  vars_varying <- c(vars_varying, grep(paste0("r[0-9]+", v), names(df), value=TRUE))
}

vars_varying <- c(vars_varying, grep(paste0("h[0-9]+", "itot"), names(df), value=TRUE))
vars_varying <- c(vars_varying, grep(paste0("h[0-9]+", "atota"), names(df), value=TRUE))
vars_varying <- c(vars_varying, grep(paste0("h[0-9]+", "atotb"), names(df), value=TRUE))
vars_varying <- c(vars_varying, grep(paste0("h[0-9]+", "atotw"), names(df), value=TRUE))
vars_varying <- c(vars_varying, grep(paste0("inw", "[0-9]+"), names(df), value=TRUE))

vars_varying_new <- vars_varying

rm(vars_varying)

# Select all constant individual-specific variables
vars_fixed <- c("hhidpn",
                "hacohort", "racohbyr",
                "rawtsamp", "raestrat", "raehsamp",
                "rabmonth", "rabyear", "rabdate", 
                "ragender", "raracem", "rahispan",
                "rameduc", "rafeduc", "raeduc",
                "raedyrs", "raedegrm",
                "radyear", "radmonth", "raddate"
)

# Create subset of data with only relevant variables
df <- df[c(vars_fixed, vars_varying_new)]

# Begin reshaping
df <- melt(df, id.vars = vars_fixed)
df$round <- as.numeric(str_extract(df$variable, "[0-9]+"))
df$variable <- sub('[0-9]+', '', df$variable)
df$variable <- gsub('^r', '',df$variable)

# Cast variables to colum names
df <- dcast(df, hhidpn + round + hacohort + racohbyr +
              rawtsamp + raestrat + raehsamp +
              rabmonth + rabyear + rabdate +
              ragender + raracem + rahispan +
              rameduc + rafeduc + raeduc +
              raedyrs + raedegrm + 
              radyear + radmonth + raddate
            ~ variable, value.var='value')

# create year variable (used in merge with FAT files)
# for round-year connections, see randhrs_P.pdf documentation file
df$year <- NA
df$year[df$round==1] <- 1992
df$year[df$round==2] <- 1994
df$year[df$round==3] <- 1996
df$year[df$round==4] <- 1998
df$year[df$round==5] <- 2000
df$year[df$round==6] <- 2002
df$year[df$round==7] <- 2004
df$year[df$round==8] <- 2006
df$year[df$round==9] <- 2008
df$year[df$round==10] <- 2010
df$year[df$round==11] <- 2012
df$year[df$round==12] <- 2014
df$year[df$round==13] <- 2016

# arrange data
df <- arrange(df, hhidpn, round)

# save dataset
saveRDS(df, file=file.path(path,"1_data/interim/rand_reshaped.rds"))
