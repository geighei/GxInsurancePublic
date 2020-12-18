# This Script -------------------------------------------------------------
# Author: Laura Zwyssig
# Goal: Convert RAND Stata file to an .rds file
# Last edited: 20.06.2018

# Preliminaries -----------------------------------------------------------
library(tidyverse)
library(haven)
library(labelled)
library(readstata13)
# clear workspace
rm(list=ls())

## RAND Data
# Remove Stata specific inofrmation ---------------------------------------
# Data source: https://ssl.isr.umich.edu/hrs/files2.php?versid=34 (RAND Version P data file)

# Load data
#path <- "~/Dropbox/GxInsurance"
#path <- "T:/econ/biroli/geighei/data/HRS/RANDfiles/randhrs1992_2016v2_STATA/"
path <- "C:/Users/parce/Desktop/GeiGhei/GxInsurance/1_data"

#path <- "/home/debian/biroli/geighei/code/GxInsurance/1_data"

#rand <- readstata13::read.dta13(file.path(path,"raw/randhrs1992_2016v2.dta"),
#                                generate.factors=T, nonint.factors = FALSE)
rand <- read_dta(file.path(path,"raw/randhrs1992_2016v2.dta"))

# Remove Stata labels
rand <- remove_labels(rand)

# Remove special Stata format
rand <- zap_formats(rand)

# Save as R data file (.rds)
saveRDS(rand, file.path(path,"interim/randhrs1992_2016v2.rds"))

# HRS - Leave-behind questionnaires
# Load data - 2006
LB06 <- read_stata(file.path(path,"raw/H06LB_R.dta"))
# generate person identifier and keep the variables on conscientiousness only
LB06 <- LB06 %>% mutate(hhidpn = 1000* as.numeric(HHID) + as.numeric(PN))
LB06$year <- 2006

names(LB06) <- gsub("KLB","",names(LB06))
consc <- c( "hhidpn", "year", "033A","033B","033C","033D","033E","033F",
            "033G","033H","033I", "033J","033K","033L","033M","033N","033O",
            "033P","033Q","033R", "033S","033T","033U","033V","033W","033X",
            "033Y","033Z")
LB06 <- LB06[,consc]

# make it compatible with years 2010/12 (new questions were added)
LB06$`033Z_2` <- NA
LB06$`033Z_3` <- NA
LB06$`033Z_4` <- NA
LB06$`033Z_5` <- NA
LB06$`033Z_6` <- NA


# Load data - 2008
LB08 <- read_stata(file.path(path,"raw/H08LB_R.dta"))
# generate person identifier and keep the variables on conscientiousness only
LB08 <- LB08 %>% mutate(hhidpn = 1000* as.numeric(HHID) + as.numeric(PN))
LB08$year <- 2008

names(LB08) <- gsub("LLB","",names(LB08))
consc <- c( "hhidpn", "year", "033A","033B","033C","033D","033E","033F",
            "033G","033H","033I", "033J","033K","033L","033M","033N","033O",
            "033P","033Q","033R", "033S","033T","033U","033V","033W","033X",
            "033Y","033Z")
LB08 <- LB08[,consc]

# make it compatible with years 2010/12 (new questions were added)
LB08$`033Z_2` <- NA
LB08$`033Z_3` <- NA
LB08$`033Z_4` <- NA
LB08$`033Z_5` <- NA
LB08$`033Z_6` <- NA


# Load data - 2010
LB10 <- read_stata(file.path(path,"raw/H10LB_R.dta"))
# generate person identifier and keep the variables on conscientiousness only
LB10 <- LB10 %>% mutate(hhidpn = 1000* as.numeric(HHID) + as.numeric(PN))
LB10$year <- 2010

names(LB10) <- gsub("MLB","",names(LB10))
consc <- c("hhidpn", "year","033A","033B","033C","033D","033E","033F","033G", 
           "033H","033I", "033J","033K","033L","033M","033N","033O","033P",
           "033Q","033R", "033S","033T","033U","033V","033W","033X","033Y",
           "033Z", "033Z_2", "033Z_3","033Z_4","033Z_5","033Z_6")
LB10 <- LB10[,consc]

# Load data - 2012
LB12 <- read_stata(file.path(path,"raw/H12LB_R.dta"))
# generate person identifier and keep the variables on conscientiousness only
LB12 <- LB12 %>% mutate(hhidpn = 1000* as.numeric(HHID) + as.numeric(PN))
LB12$year <- 2012

names(LB12) <- gsub("NLB","",names(LB12))
consc <- c("hhidpn", "year","033A","033B","033C","033D","033E","033F","033G", 
           "033H","033I", "033J","033K","033L","033M","033N","033O","033P",
           "033Q","033R", "033S","033T","033U","033V","033W","033X","033Y",
           "033Z", "033Z_2", "033Z_3","033Z_4","033Z_5","033Z_6")
LB12 <- LB12[,consc]

# Load data - 2014
LB14 <- read_stata(file.path(path,"raw/H14LB_R.dta"))
# generate person identifier and keep the variables on conscientiousness only
LB14 <- LB14 %>% mutate(hhidpn = 1000* as.numeric(HHID) + as.numeric(PN))
LB14$year <- 2014
consc14 <- c("hhidpn", "year", "OLB031A", "OLB031B", "OLB031C", "OLB031D", "OLB031E", "OLB031F", 
             "OLB031G","OLB031H","OLB031I", "OLB031J", "OLB031K", "OLB031L", "OLB031M", "OLB031N",
             "OLB031O", "OLB031P", "OLB031Q", "OLB031R", "OLB031S", "OLB031T", "OLB031U", "OLB031V",
             "OLB031W", "OLB031X", "OLB031Y", "OLB031Z_1", "OLB031Z_2", "OLB031Z_3", "OLB031Z_4", 
             "OLB031Z_5","OLB031Z_6")
LB14 <- LB14[,consc14]
# Rename to match the other variables from prevoius rounds
# Now with smart quotes 
consc <- c('hhidpn', 'year','033A','033B','033C','033D','033E','033F','033G', 
           '033H','033I', '033J','033K','033L','033M','033N','033O','033P',
           '033Q','033R', '033S','033T','033U','033V','033W','033X','033Y',
           '033Z', '033Z_2', '033Z_3','033Z_4','033Z_5','033Z_6')
for(i in 1:length(consc)){
  colnames(LB14)[i] <- consc[i]
}

# Load data - 2016
LB16 <- read_stata(file.path(path,"raw/H16LB_R.dta"))
# generate person identifier and keep the variables on conscientiousness only
LB16 <- LB16 %>% mutate(hhidpn = 1000* as.numeric(HHID) + as.numeric(PN))
LB16$year <- 2016
consc16 <- c("hhidpn", "year", "PLB031A", "PLB031B", "PLB031C", "PLB031D", "PLB031E", "PLB031F", 
             "PLB031G", "PLB031H", "PLB031I", "PLB031J", "PLB031K", "PLB031L", "PLB031M", "PLB031N",
             "PLB031O", "PLB031P", "PLB031Q", "PLB031R", "PLB031S", "PLB031T", "PLB031U", "PLB031V", 
             "PLB031W", "PLB031X", "PLB031Y", "PLB031Z_1", "PLB031Z_2", "PLB031Z_3", "PLB031Z_4", 
             "PLB031Z_5", "PLB031Z_6")
LB16 <- LB16[,consc16]
# Rename to match the other variables from prevoius rounds
for(i in 1:length(consc)){
  colnames(LB16)[i] <- consc[i]
}


# bind everything together
LB <- rbind(LB06, LB08, LB10, LB12, LB14, LB16)

# Remove Stata labels
LB <- remove_labels(LB)

# Remove special Stata format
LB <- zap_formats(LB)

# Save as R data file (.rds)
saveRDS(LB, file = file.path(path,"interim/consc.rds"))

# HRS smoke variables cleaned and processed by Michelle in the group folder "gheighei/data/HRS/temp1/SMOKE_HRS_ALL.dta"
# load data
smoke <- read_stata(file.path(path,"raw/SMOKE_HRS_ALL.dta"))

# keep only variables relevant for merging and analysis
smoke <- smoke[c("HHID","PN","YEAR","wave","cigsDay","cigsPastDay","cigsAll_Current","cigsAll_Past")]

# Remove Stata labels
smoke <- remove_labels(smoke)

# Remove special Stata format
smoke <- zap_formats(smoke)

# Save as R data file (.rds)
saveRDS(smoke, file = file.path(path,"interim/SMOKE_HRS_ALL.rds"))
