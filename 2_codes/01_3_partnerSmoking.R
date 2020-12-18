# This Script -------------------------------------------------------------
# Author: Pia Arce
# Goal: Add smoking partner variable from RAND data 
# Last edited: 

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

# Partner ever smoke and smoking now: Smoken Smokev
vars_varying_base <- c("smokev", "smoken")

# Select all time-varying variables
vars_varying <- c()
for (v in vars_varying_base) {
  vars_varying <- c(vars_varying, grep(paste0("s[0-9]+", v), names(df), value=TRUE))
}
ids <- c("hhidpn")
# Create subset of data with only relevant variables
df <- df[c(ids,vars_varying)]
df <- melt(df, id.vars = ids)
df$round <- as.numeric(str_extract(df$variable, "[0-9]+"))
df$variable <- sub('[0-9]+', '', df$variable)
df$variable <- gsub('^s', '',df$variable)

# Cast variables to colum names
df <- dcast(df, hhidpn + round 
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
saveRDS(df, file=file.path(path,"1_data/interim/rand_partnerData.rds"))

