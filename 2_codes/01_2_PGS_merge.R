# This Script -------------------------------------------------------------
# Author: Pia Arce
# Goal: Merge HRS polygenic scores to reshaped RAND data file
# Last edited: 28.08.2020

# Preliminaries -----------------------------------------------------------

xlibrary <- c("tidyverse", "haven", "readxl", "stringr", "labelled" )
# load libraries and automatically install all packages if missing
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = xlibrary)
# clear workspace
rm(list=ls())

# Load and prepare data ---------------------------------------------------
# Load processed (reduced & reshaped) rand version P data file

# Set here the path for "interim" folder
path <- "C:/Users/parce/Desktop/GeiGhei/GxInsurance/"
df_rand <- readRDS(file.path(path,"1_data/interim/rand_reshaped.rds"))

# Load publicly available polygenic scores----------------
# Set here the path for publicly available polygenic scores
pgs <- read_stata(file.path(path,"1_data/raw/PGENSCOREv2/PGENSCOREE_R.dta"))

# Remove Stata labels
pgs <- remove_labels(pgs)
# Remove stata specific format
pgs <- zap_formats(pgs)
# Create hhidpn
pgs <- pgs %>% mutate(hhidpn = 1000* as.numeric(HHID) + as.numeric(PN))

# Load the new PGS constructed by Amr here https://github.com/geighei/HRS_cleaning/tree/master/PGS_construct -----------
pgs_self <- read.csv(file.path(path,"1_data/raw/PGS_HRS/all_pgs_hrs_combined_threshold.csv"))
# Create hhidpn
pgs_self <- pgs_self %>% mutate(hhidpn = 1000* as.numeric(HHID) + as.numeric(PN))

# Load PGS for risk and and educational attainment as well and merge with new PGS
pgs_more <- read.delim(file.path(path,"1_data/raw/PGS_HRS/hrs_prsice_pgs.txt"), sep = " ")

pgs_more <- pgs_more[c("FID","IID","risk_gwas_pval_1","risk_pc1_pval_1","ea_pval_1","cog_pval_1","noncog_pval_1")]
names(pgs_more)[2] <- "ID"

pgs_more <- full_join(pgs_self,pgs_more,by = c("FID","ID"))

# merge all the PGS
pgs_all <- left_join(pgs_more, pgs, by=c("hhidpn"))
tempcheck <- anti_join(pgs_more, pgs, by=c("hhidpn"))
dim(tempcheck)
rm(tempcheck)

# Merge all PGS data to reshaped RAND data
df <- left_join(df_rand, pgs_all, by=c("hhidpn"))

# Give intuitive definition to RAND variables -----------------------------
df$ragender <- as.numeric(df$ragender)
df <- df %>% mutate(female = ragender - 1)

df <- df %>% mutate(age = agey_e)

# Merge Smoking data, which was missing from RAND data, but Michelle has constructed in the group folder
# "gheighei/data/HRS/temp1/SMOKE_HRS_ALL.dta"
df_smoke <- readRDS(file.path(path,"1_data/interim/SMOKE_HRS_ALL.rds"))

df_smoke <- df_smoke %>% mutate(hhidpn = 1000* as.numeric(HHID) + as.numeric(PN))

names(df_smoke)[3] <- "year"
names(df_smoke)[5] <- "cigpdnow"
names(df_smoke)[6] <- "cigpdpast"
names(df_smoke)[7] <- "infl_cigpdnow"
names(df_smoke)[8] <- "infl_cigpdpast"


df_smoke <- df_smoke[c("hhidpn","year","cigpdnow","cigpdpast","infl_cigpdnow","infl_cigpdpast")]

df_final <- left_join(df, df_smoke, by = c("hhidpn","year"))
tempchck <- anti_join(df, df_smoke, by = c("hhidpn","year"))
dim(tempchck)
rm(tempchck)
saveRDS(df_final, file = file.path(path,"1_data/interim/hrs.rds"))

# Load new PGS constructed by Jeremy (2020)----------------------------------------
pgspath <- "T:/econ/biroli/geighei/data/HRS/PGS/PRScs"
pgs_maxCPD <- read.delim(file.path(pgspath,"/maxCPD/maxCPD_PGS.profile"), header = TRUE, sep = "")
pgs_maxCPD <- pgs_maxCPD[c("FID","IID","SCORE")]
names(pgs_maxCPD)[names(pgs_maxCPD) == "SCORE"] <- "maxCPD_1"

pgs_smokeInit <- read.delim(file.path(pgspath,"/smokeInit/original_qc/smokeInit_PGS.sscore"), header = TRUE, sep = "")
names(pgs_smokeInit)[names(pgs_smokeInit) == "X.FID"] <- "FID"
names(pgs_smokeInit)[names(pgs_smokeInit) == "SCORE1_AVG"] <- "SCORE"
pgs_smokeInit <- pgs_smokeInit[c("FID","IID","SCORE")]
names(pgs_smokeInit)[names(pgs_smokeInit) == "SCORE"] <- "smokeInit_1"

pgs_smokeInit2 <- read.delim(file.path(pgspath,"/smokeInit/smokeInit_PGS.profile"), header = TRUE, sep = "")
pgs_smokeInit2 <- pgs_smokeInit2[c("FID","IID","SCORE")]
names(pgs_smokeInit2)[names(pgs_smokeInit2) == "SCORE"] <- "smokeInitSmall_1"

pgs_BMI <- read.delim(file.path(pgspath,"/bmi/bmi_PGS.profile"), header = TRUE, sep = "")
pgs_BMI <- pgs_BMI[c("FID","IID","SCORE")]
names(pgs_BMI)[names(pgs_BMI) == "SCORE"] <- "Bmi_1"

###
pgs_newSmoke <- full_join(pgs_maxCPD ,pgs_smokeInit,by = c("FID","IID"))
temp <- anti_join(pgs_maxCPD,pgs_smokeInit,by = c("FID","IID"))
dim(temp)
rm(temp)

pgs_newSmoke <- full_join(pgs_newSmoke ,pgs_smokeInit2,by = c("FID","IID"))
temp <- anti_join(pgs_newSmoke, pgs_smokeInit2,by = c("FID","IID"))
dim(temp)
rm(temp)

pgs_newSmoke <- full_join(pgs_BMI,pgs_newSmoke, by = c("FID","IID") )
temp <- anti_join(pgs_BMI,pgs_newSmoke, by = c("FID","IID") )
dim(temp)
rm(temp)

names(pgs_newSmoke)[names(pgs_newSmoke) == "IID"] <- "iid"
pgs_newSmoke$iid <- as.numeric(pgs_newSmoke$iid)

# Extract middle-step ID
pathhrsid <- "T:/econ/biroli/geighei/data/HRS/genomeraw/archive"
hrs_id <- read_stata(file.path(pathhrsid,"/HRS_id.dta"))
hrs_id <- hrs_id[c("HHID","PN","LOCAL_ID")]
names(hrs_id)[names(hrs_id) == "LOCAL_ID"] <- "iid"
hrs_id$iid <- as.numeric(hrs_id$iid) 
hrs_id <- as.data.frame(hrs_id)

# Merge by iid
#pgs_allNew <- merge(pgs_newSmoke, hrs_id) # 15566
pgs_allNew <- left_join(pgs_newSmoke, hrs_id, by = c("iid")) # 15566
tempcheck <- anti_join(pgs_newSmoke, hrs_id, by=c("iid"))
dim(tempcheck)
rm(tempcheck)
pgs_allNew <- pgs_allNew[complete.cases(pgs_allNew$smokeInit_1), ] # 15566

pgs_allNew <- pgs_allNew %>% mutate(hhidpn = 1000* as.numeric(HHID) + as.numeric(PN))
#pgs_allNew$hhidpn2 <- as.numeric(paste0(pgs_allNew$HHID, 0, pgs_allNew$PN))

#df <- df_final
df_final2 <- left_join( df_final,pgs_allNew, by=c("hhidpn")) #449940
tempcheck <- anti_join(df_final, pgs_allNew, by=c("hhidpn")) 
#[1] 263148    274
dim(tempcheck)
rm(tempcheck)

# Standardize new PGSs
df_final2 <- df_final2 %>% mutate_at(c("smokeInit_1","smokeInitSmall_1", "maxCPD_1","Bmi_1"), ~(scale(.) %>% as.numeric))

# New PGS are flipped 
# NOTE: Change this if it is fixed 
df_final2$maxCPD_1 <- -1 * df_final2$maxCPD_1
df_final2$Bmi_1 <- -1 * df_final2$Bmi_1
# Save dataset ------------------------------------------------------------
saveRDS(df_final2, file = file.path(path,"1_data/interim/hrs2.rds"))

