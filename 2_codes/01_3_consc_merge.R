# This Script -------------------------------------------------------------
# Author: Regina Seibel
# Goal: Merge conscientiousness questions from LB to hrs data file
# Last edited: 22.02.2019

# Preliminaries -----------------------------------------------------------
library(tidyverse)
library(haven)
library(readxl)
library(stringr)
library(labelled)


# clear workspace
rm(list=ls())

# load data
#path <- "../../../../Dropbox/GxInsurance/"
#path <- "T:/econ/biroli/geighei/code/GxInsurance"
path <- "C:/Users/parce/Desktop/GeiGhei/GxInsurance/"

#hrs <- readRDS(file.path(path,"1_data/processed/hrs_final.rds"))
hrs <- readRDS(file.path(path,"1_data/processed/hrs_final2.rds"))

consc <- readRDS(file.path(path,"1_data/interim/consc.rds"))

# create the conscientiousness indicator
consc <- filter(consc,!is.na(as.integer(get("033D"))))
consc <- filter(consc,!is.na(as.integer(get("033Z"))))


rev <- c("033D", "033H", "033M", "033Z")
new <- consc[rev]
new <- apply(new, 2,function(x) 5 - x) 
consc_only <- cbind(new, consc["033T"])
consc_only$mean <- apply(consc_only,1,mean)
consc_only <- cbind(consc[c("hhidpn","year")],consc_only)

# calculate the average over years (if mulitple observation per individual)
consc_only <- aggregate(consc_only[3:8], by = list(consc_only$hhidpn), mean)

colnames(consc_only) <- c("hhidpn", "33d", "33h", "33m", "33z", "33t", "consc_mean")

# merge and keep all hrs observations
hrs_consc <- left_join(hrs,consc_only, by = c("hhidpn"))

# save dataset
#saveRDS(hrs_consc, file = file.path(path,"1_data/processed/hrs_consc.rds"))
saveRDS(hrs_consc, file = file.path(path,"1_data/processed/hrs_consc2.rds"))
