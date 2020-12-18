# This Script -------------------------------------------------------------
# Author: Laura Zwyssig, Pietro Biroli, Regina Seibel
# Goal: Contains functions used in 02_1_analysis.Rmd
# Last edited: March 2020

# Preliminaries -----------------------------------------------------------
xlibrary <- c("readr", "haven", "stringr", "stargazer", "plm", "sandwich", "lmtest", "car", "ggplot2", "scales", "multcomp",
        "broom", "rdrobust", "xtable", "foreign", "lfe", "starpolishr", "rlang", "readstata13","tidyverse", "dplyr")

lapply(xlibrary, require, character.only = TRUE)

# Function: bold.sometext is for making some cells bold in the tables later on
bold.sometext <-
  function(x) gsub('BOLD(.*)',paste('\\\\textbf{\\1','}'),x)


# Function: Creates a dataframe for analysis ------------------------------
# Two special inputs: outcome for different y-variables, indep_var for different x-variables
# on which health shock effect might differ
create_df <- function(df_raw, min_sample_age, max_sample_age, ptu, indep_var , outcome = "smoken",  high_cutoff = 0.33333, only_baseline = TRUE, drop_6566_cv = TRUE, drop_6566_nsr = FALSE, drop_non_gen = TRUE){
  ## Inputs:
  # df_raw:             dataset (string)
  # min_sample_age      minumum age to be part of sample (integer, e.g. 60)
  # max_sample_age      maximum age to be part of sample (integer, e.g. 70)
  # ptu                 percentage of waves without insurance to be considered unisurance (integere, e.g. 100)
  # outcome             main outcome variable (e.g. "smoken")
  # indep_var           main PGS (e.g. = "si_1")
  # high_cutoff         percentile cutoff above wich you're considered high PGS (default = 0.33)
  # only_baseline       = TRUE if dropping never smokers at baseline
  # drop_6566_cv        = TRUE if dropping individuals with reported health shocks at age 65/66
  # drop_6566_nsr       = TRUE if dropping non-smoking related health shocks at age 65/66
  # drop_non_gen        = TRUE if dropping missing values in outcome, indep var, and cardiovascular disease (cv)

  # Individuals need to be observed for two different waves minimally
  df <- df_raw %>%
    filter(age >= min_sample_age,
           age <= max_sample_age) %>%
    group_by(hhidpn) %>%
    filter(n() >= 2) %>%
    ungroup()

  # Use only individuals with non-missing data in any variable needed for the analysis
  # Exception: unins, post65 (variables created further below)
  if(drop_non_gen){
    df <- df %>% filter(!is.na(get(indep_var)),
                        !is.na(cv),
                        !is.na(get(outcome))) # need get() so it realizes that this is a variable
  }

  # Exclude the ones who reported having never smoked at baseline (first observation per hhidpn)
  if (outcome == "smoken") {
  if(only_baseline) {
    smk_baseline <- df %>%
    group_by(hhidpn) %>%
    arrange(year) %>%
    filter(row_number() == 1) %>%
    filter(smokev == 1) %>%
    ungroup() %>%
    select(hhidpn)
  stopifnot(anyDuplicated(smk_baseline$hhidpn) == 0)

  df <- df %>% filter(hhidpn %in% smk_baseline$hhidpn)
  }
  }

  # Create persistent uninsured dummy, Medicare eligibility dummy, and pre-65 dummy
  uninsured_hhidpn <- df %>%
    group_by(hhidpn) %>%
    filter(age < 65) %>%
    mutate(unins = ifelse(sum(uninsured) >= ptu/100*n(), 1, 0)) %>%
    ungroup() %>%
    filter(unins == 1) %>%
    select(hhidpn) %>%
    distinct() %>%
    pull(hhidpn)

  df <- df %>%
    mutate(unins = ifelse(hhidpn %in% uninsured_hhidpn, 1, 0),
           post65 = ifelse(age >= 65, 1, 0),
           pre65 = ifelse(age < 65, 1, 0))

  # Keep individuals with non-missing pre-65 insurance status
  df <- df %>% filter(!is.na(unins),
                      !is.na(post65))

  # Exclude individuals with reported health shocks at age 65/66
  if (drop_6566_cv) {
    df <- df %>%
      group_by(hhidpn) %>%
      mutate(cv_6566 = ifelse(
        (
          age == 65 |
            age == 66
        ) & cv == 1, 1, 0)) %>%
      mutate(cv_6566 = max(cv_6566, na.rm = TRUE)) %>%
      ungroup() %>%
      filter(cv_6566 != 1) %>%
      select(-cv_6566)
  }

  # Exclude non-smoking related health shocks at age 65/66
  if (drop_6566_nsr) {
    df <- df %>%
      group_by(hhidpn) %>%
      mutate(nsr_6566 = ifelse(
        (
          age == 65 |
            age == 66
        ) & nsr == 1, 1, 0)) %>%
      mutate(nsr_6566 = max(nsr_6566, na.rm = TRUE)) %>%
      ungroup() %>%
      filter(nsr_6566 != 1) %>%
      select(-nsr_6566)
  }

  # Create low and high PGS dummies and lagged CV dummy
    df <- df %>%
    mutate(high_pgs = ifelse(get(indep_var) >  quantile(get(indep_var), high_cutoff , na.rm = TRUE), 1, 0),
           low_pgs  = ifelse(get(indep_var) <= quantile(get(indep_var), high_cutoff , na.rm = TRUE), 1, 0)) %>%
    group_by(hhidpn) %>%
    mutate(lag1_cv = dplyr::lag(cv, n = 1, default = NA),
           lag2_cv = dplyr::lag(cv, n = 2, default = NA),
           cv_2w = ifelse(
             lag1_cv == 1 |
               cv==1, 1, 0
           ),
           cv_3w = ifelse(
             lag1_cv == 1 |
               lag2_cv == 1 |
               cv == 1, 1, 0
           )) %>%
    ungroup()

  # Create age polynomials (of demeaned age to avoid multicollinearity problems)
  df <- df %>% mutate(agem = agem_e - mean(agem_e, na.rm = TRUE),
                                      agem2 = agem^2,
                                      agem3 = (agem^3/100))
  return(df)
}
# Function: Creates a dataframe for analysis ------------------------------
# Two special inputs: outcome for different y-variables, indep_var for different x-variables
# on which health shock effect might differ
create_df_median <- function(df_raw, min_sample_age, max_sample_age, ptu, indep_var , outcome = "smoken",   only_baseline = TRUE, drop_6566_cv = TRUE, drop_6566_nsr = FALSE, drop_non_gen = TRUE){
  ## Inputs:
  # df_raw:             dataset (string)
  # min_sample_age      minumum age to be part of sample (integer, e.g. 60)
  # max_sample_age      maximum age to be part of sample (integer, e.g. 70)
  # ptu                 percentage of waves without insurance to be considered unisurance (integere, e.g. 100)
  # outcome             main outcome variable (e.g. "smoken")
  # indep_var           main PGS (e.g. = "si_1")
  # high_cutoff         percentile cutoff above wich you're considered high PGS (default = 0.33)
  # only_baseline       = TRUE if dropping never smokers at baseline
  # drop_6566_cv        = TRUE if dropping individuals with reported health shocks at age 65/66
  # drop_6566_nsr       = TRUE if dropping non-smoking related health shocks at age 65/66
  # drop_non_gen        = TRUE if dropping missing values in outcome, indep var, and cardiovascular disease (cv)
  
  
  # Individuals need to be observed for two different waves minimally
  df <- df_raw %>%
    filter(age >= min_sample_age,
           age <= max_sample_age) %>%
    group_by(hhidpn) %>%
    filter(n() >= 2) %>%
    ungroup()
  
  # Use only individuals with non-missing data in any variable needed for the analysis
  # Exception: unins, post65 (variables created further below)
  if(drop_non_gen){
    df <- df %>% filter(!is.na(get(indep_var)),
                        !is.na(cv),
                        !is.na(get(outcome))) # need get() so it realizes that this is a variable
  }
  
  # Exclude the ones who reported having never smoked at baseline (first observation per hhidpn)
  if (outcome == "smoken") {
    if(only_baseline) {
      smk_baseline <- df %>%
        group_by(hhidpn) %>%
        arrange(year) %>%
        filter(row_number() == 1) %>%
        filter(smokev == 1) %>%
        ungroup() %>%
        select(hhidpn)
      stopifnot(anyDuplicated(smk_baseline$hhidpn) == 0)
      
      df <- df %>% filter(hhidpn %in% smk_baseline$hhidpn)
    }
  }
  
  # Create persistent uninsured dummy, Medicare eligibility dummy, and pre-65 dummy
  uninsured_hhidpn <- df %>%
    group_by(hhidpn) %>%
    filter(age < 65) %>%
    mutate(unins = ifelse(sum(uninsured) >= ptu/100*n(), 1, 0)) %>%
    ungroup() %>%
    filter(unins == 1) %>%
    select(hhidpn) %>%
    distinct() %>%
    pull(hhidpn)
  
  df <- df %>%
    mutate(unins = ifelse(hhidpn %in% uninsured_hhidpn, 1, 0),
           post65 = ifelse(age >= 65, 1, 0),
           pre65 = ifelse(age < 65, 1, 0))
  
  # Keep individuals with non-missing pre-65 insurance status
  df <- df %>% filter(!is.na(unins),
                      !is.na(post65))
  
  # Exclude individuals with reported health shocks at age 65/66
  if (drop_6566_cv) {
    df <- df %>%
      group_by(hhidpn) %>%
      mutate(cv_6566 = ifelse(
        (
          age == 65 |
            age == 66
        ) & cv == 1, 1, 0)) %>%
      mutate(cv_6566 = max(cv_6566, na.rm = TRUE)) %>%
      ungroup() %>%
      filter(cv_6566 != 1) %>%
      select(-cv_6566)
  }
  
  # Exclude non-smoking related health shocks at age 65/66
  if (drop_6566_nsr) {
    df <- df %>%
      group_by(hhidpn) %>%
      mutate(nsr_6566 = ifelse(
        (
          age == 65 |
            age == 66
        ) & nsr == 1, 1, 0)) %>%
      mutate(nsr_6566 = max(nsr_6566, na.rm = TRUE)) %>%
      ungroup() %>%
      filter(nsr_6566 != 1) %>%
      select(-nsr_6566)
  }
  
  # Create low and high PGS dummies and lagged CV dummy
  df <- df %>%
    mutate(high_pgs = ifelse(get(indep_var) >  median(get(indep_var) , na.rm = TRUE), 1, 0),
           low_pgs  = ifelse(get(indep_var) <= median(get(indep_var) , na.rm = TRUE), 1, 0)) %>%
    group_by(hhidpn) %>%
    mutate(lag1_cv = dplyr::lag(cv, n = 1, default = NA),
           lag2_cv = dplyr::lag(cv, n = 2, default = NA),
           cv_2w = ifelse(
             lag1_cv == 1 |
               cv==1, 1, 0
           ),
           cv_3w = ifelse(
             lag1_cv == 1 |
               lag2_cv == 1 |
               cv == 1, 1, 0
           )) %>%
    ungroup()
  
  # Create age polynomials (of demeaned age to avoid multicollinearity problems)
  df <- df %>% mutate(agem = agem_e - mean(agem_e, na.rm = TRUE),
                      agem2 = agem^2,
                      agem3 = agem^3)
  return(df)
}
# Function: Counts individuals --------------------------------------------
count_individuals <- function(df, by_round=FALSE, simple=FALSE){
  print_count <- function(df_hhidpn) {
    print(df_hhidpn %>% distinct() %>% tally())
  }
  if (simple) {
    df %>% select(hhidpn) %>% distinct() %>% tally()
  } else if (!simple) {
    if (!by_round) {
      print("All")
      print_count(df %>% select(hhidpn))
      print("With genetic data")
      print_count(df %>% filter(!is.na(si_1)) %>% select(hhidpn))
    } else {
      print("All")
      print_count(df %>% group_by(round) %>% select(round, hhidpn))
      print("With genetic data")
      print_count(df %>% filter(!is.na(si_1)) %>%
                    group_by(round) %>% select(round, hhidpn))
    }
  }
}

# Function: Creates a samplesize table ------------------------------------
table_counts <- function(df){
  table_counts <- tibble(number_individuals = rep(c("all", "with genetic", "with low PGS", "with high PGS"), 3))

  vars_shocks <- c("cv", "cv", "cv")
  for (i in c(1, 5, 9)){
    table_counts[["final_sample"]][i] <- count_individuals(df, simple = TRUE)$n
    table_counts[["final_sample"]][i + 1] <- count_individuals(df %>% filter(!is.na(si_1)), simple = TRUE)$n
    table_counts[["final_sample"]][i + 2] <- count_individuals(df %>% filter(high_pgs==0), simple = TRUE)$n
    table_counts[["final_sample"]][i + 3] <- count_individuals(df %>% filter(high_pgs==1), simple = TRUE)$n

    j <- ifelse(i==1, 1,
                ifelse(i==5, 2, 3))

    df$var_of_interest <- df[[vars_shocks[j]]]
    # Number of individuals with health shock pre age 65 and persistently uninsured before
    table_counts[["shock_pre_65_uninsured"]][i] <- count_individuals(df %>%
                                                                       group_by(hhidpn) %>%
                                                                       filter(age < 65, unins == 1) %>%
                                                                       summarise(max_voi = max(var_of_interest)) %>%
                                                                       filter(max_voi == 1) %>%
                                                                       ungroup(), simple = TRUE)$n

    # Number of individuals with health shock post age 65 and persistently uninsured before + genetic data
    table_counts[["shock_pre_65_uninsured"]][i + 1] <- count_individuals(df %>%
                                                                           filter(!is.na(si_1)) %>%
                                                                           group_by(hhidpn) %>%
                                                                           filter(age < 65, unins == 1) %>%
                                                                           summarise(max_voi = max(var_of_interest)) %>%
                                                                           filter(max_voi == 1) %>%
                                                                           ungroup(), simple = TRUE)$n

    # Number of individuals with health shock post age 65 and persistently uninsured before + low PGS
    table_counts[["shock_pre_65_uninsured"]][i + 2] <- count_individuals(df %>%
                                                                           filter(high_pgs==0) %>%
                                                                           group_by(hhidpn) %>%
                                                                           filter(age < 65, unins == 1) %>%
                                                                           summarise(max_voi = max(var_of_interest)) %>%
                                                                           filter(max_voi == 1) %>%
                                                                           ungroup(), simple = TRUE)$n

    # Number of individuals with health shock post age 65 and persistently uninsured before + high PGS
    table_counts[["shock_pre_65_uninsured"]][i + 3] <- count_individuals(df %>%
                                                                           filter(high_pgs==1) %>%
                                                                           group_by(hhidpn) %>%
                                                                           filter(age < 65, unins == 1) %>%
                                                                           summarise(max_voi = max(var_of_interest)) %>%
                                                                           filter(max_voi == 1) %>%
                                                                           ungroup(), simple = TRUE)$n

    # Number of individuals with health shock post age 65 and persistently uninsured before
    table_counts[["shock_post_65_uninsured"]][i] <- count_individuals(df %>%
                                                                        filter(unins == 1) %>%
                                                                        group_by(hhidpn) %>%
                                                                        filter(age >= 65) %>%
                                                                        summarise(max_voi = max(var_of_interest)) %>%
                                                                        filter(max_voi == 1) %>%
                                                                        ungroup(), simple = TRUE)$n

    # Number of individuals with health shock post age 65 and persistently uninsured before + genetic data
    table_counts[["shock_post_65_uninsured"]][i + 1] <- count_individuals(df %>%
                                                                            filter(unins == 1, !is.na(si_1)) %>%
                                                                            group_by(hhidpn) %>%
                                                                            filter(age >= 65) %>%
                                                                            summarise(max_voi = max(var_of_interest)) %>%
                                                                            filter(max_voi == 1) %>%
                                                                            ungroup(), simple = TRUE)$n

    # Number of individuals with health shock post age 65 and persistently uninsured before + low PGS
    table_counts[["shock_post_65_uninsured"]][i + 2] <- count_individuals(df %>%
                                                                            filter(unins == 1, high_pgs==0) %>%
                                                                            group_by(hhidpn) %>%
                                                                            filter(age >= 65) %>%
                                                                            summarise(max_voi = max(var_of_interest)) %>%
                                                                            filter(max_voi == 1) %>%
                                                                            ungroup(), simple = TRUE)$n

    # Number of individuals with health shock post age 65 and persistently uninsured before + high PGS
    table_counts[["shock_post_65_uninsured"]][i + 3] <- count_individuals(df %>%
                                                                            filter(unins == 1, high_pgs==1) %>%
                                                                            group_by(hhidpn) %>%
                                                                            filter(age >= 65) %>%
                                                                            summarise(max_voi = max(var_of_interest)) %>%
                                                                            filter(max_voi == 1) %>%
                                                                            ungroup(), simple = TRUE)$n

  }

  stargazer(table_counts, summary=FALSE, type = "text", rownames = FALSE)
}

# Function: Make sure that there are no missing values for a var ----------
check_no_na <- function(var_vector) {
  stopifnot(sum(is.na(var_vector)) == 0)
}

# Function: Average cessation rate among baseline smokers ------------------
mean_cess_baseline_smk <- function(df_original, df_temp, return_df = FALSE, return_se = FALSE, pgs = "si_1") {

    # Get hhidpn of all baseline smokers in
  # given dataset
  if (pgs == "si_1") {
    outcome = "smoken"
  } else if (pgs == "dpw_1") {
    outcome = "drink"
  }
    baseline <- df_temp %>%
      filter(get(outcome) == 1) %>%
      pull(hhidpn)

  # Create cessation rate
  # (based on full dataset passed to outer function)
  df_cess <- df_original %>%
    filter(hhidpn %in% baseline) %>%
    group_by(hhidpn) %>%
    arrange(year) %>%
    mutate(lag = dplyr::lag(get(outcome), n = 1, default = NA)) %>%
    mutate(cess = ifelse(
      lag == 1 &
        get(outcome) == 0,
      1, 0)) %>%
    ungroup()

  # Get average cessation rate
  check_no_na(df_cess$cess)
  avg_cess_rate <- round(
    df_cess %>% filter(!is.na(cess)) %>%
      summarize(mean_cess = mean(cess)) %>%
      mutate(mean_cess = mean_cess * 100),
    digits = 2)


  # Get standard deviation of cessation rate
  sd_cess_rate <- round(
    df_cess %>% filter(!is.na(cess)) %>%
      summarize(sd_cess = sd(cess)),
    digits = 2)


  # Get number of observations (for calculating standard error later on)
  no_cess_rate <- df_cess %>% filter(!is.na(cess)) %>%
    summarize(no_cess = n())


  if (return_df) {
    return(list("rate" = avg_cess_rate, "df" = df_cess))
  } else if(return_se) {
    return(list("rate" = avg_cess_rate, "sd" = sd_cess_rate, "no" = no_cess_rate))
  } else {
    return(avg_cess_rate)
  }

}


# Function: Rounded mean and rounded standard deviation pasted ------------
paste_mean_std <- function(df_sum, col) {
  check_no_na(df_sum[[col]])
  return(paste0(round(mean(df_sum[[col]]), digits = 2),
                " (", round(sd(df_sum[[col]]), digits = 2), ")"))
}

# Function: Rounded mean in percent ---------------------------------------
mean_perc <- function(df_sum, col) {
  check_no_na(df_sum[[col]])
  return(round(mean(df_sum[[col]]) * 100, digits = 2))
}

# Function: T test --------------------------------------------------------
perform_ttest <- function(left_df, right_df, col) {
  return(round(t.test(left_df[[col]], right_df[[col]])$p.value, digits = 4))
}


# Function: Creates a descriptive table split by high-PGS dummy -----------
descriptive_table_pgs <- function(df, name){
  # Create skeleton of table ------------------------------------------------
  # Will later on be filled with values
  table <- tibble(variable = c("",
                               "Age (baseline)",
                               "Smoking PGS",
                               "Years of education",
                               "Income (nominal \\$ 1000)",
                               "No. waves present",
                               "",
                               "Female",
                               "Smoking (baseline)",
                               "Persistently uninsured",
                               "CV health shock",
                               "Avg. cessation rate (baseline smokers)",
                               "No. of individuals"),
                  All = "",
                  lowPGS = "",
                  highPGS = ""
  )

  # Create baseline datasets ------------------------------------------------
  # "All"
  baseline_all <- df %>% group_by(hhidpn) %>%
    arrange(year) %>%
    mutate(no_waves = n(),
           max_cv = max(cv)) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    select(hhidpn,
           high_pgs,
           age,
           si_1,
           raedyrs,
           iearn,
           no_waves,
           female,
           smoken,
           uninsured,
           max_cv,
           unins) %>%
    mutate(iearn = iearn/1000)

  # "low PGS"
  baseline_low <- baseline_all %>%
    filter(high_pgs == 0)

  # "high PGS"
  baseline_high <- baseline_all %>%
    filter(high_pgs == 1)

  # Fill columns ------------------------------------------------------------
  for(i in 2:4){
    # In order of columns
    groups <- list(NA, baseline_all, baseline_low, baseline_high)
    # Select data for this column
    baseline_df <- groups[[i]]


    # Mean and standard deviations --------------------------------------------
    # Age (baseline)
    sum_age <- paste_mean_std(baseline_df, "age")

    # PGS
    sum_pgs <- paste_mean_std(baseline_df, "si_1")

    # Years of education
    # This variable contains missing values!
    sum_edyrs <- paste_mean_std(baseline_df %>%
      filter(!is.na(raedyrs)), "raedyrs") #dropping some NAs in the education variable

    # Individual income (nominal dollars)
    sum_earnings <- paste_mean_std(baseline_df, "iearn")

    # No. waves present
    sum_no_waves <- paste_mean_std(baseline_df, "no_waves")


    # Mean in percent ---------------------------------------------------------
    # Female
    sum_female <- mean_perc(baseline_df, "female")

    # Smoking (baseline)
    sum_smoking <- mean_perc(baseline_df, "smoken")

    # Persistently uninsured
    sum_unins <- mean_perc(baseline_df, "unins")

    # CV health shock
    sum_cv <- mean_perc(baseline_df, "max_cv")

    # No. of individuals
    n_individuals <- count_individuals(baseline_df, simple = TRUE)

    # Add calculated statistics to table
    # CAUTION: If a variable is added, moved, or deleted in the following line,
    # the row numbers in the "average cessation rate" just below might also needs to be updated
    column <- c("Mean (SD)",
                sum_age,
                sum_pgs,
                sum_edyrs,
                sum_earnings,
                sum_no_waves,
                "\\%",
                sum_female,
                sum_smoking,
                sum_unins,
                sum_cv,
                "",
                n_individuals)
    table[[i]] <- column
  }
  
  # Average cessation rate among baseline smokers -----------------------------------------
  # Calculated and added outside of above loop
  
  # COLUMN "All"
  table[12, 2] <- mean_cess_baseline_smk(df, baseline_all)
  
  
  # COLUMN "low PGS"
  low_cess <- mean_cess_baseline_smk(df, baseline_low, return_df = TRUE)
  table[12, 3] <- low_cess$rate
  
  
  # COLUMN "high PGS"
  high_cess <- mean_cess_baseline_smk(df, baseline_high, return_df = TRUE)
  table[12, 4] <- high_cess$rate
  
  # Change column names for table
  names(table)[1] <- "BOLD"
  names(table)[2] <- "BOLD All"
  names(table)[3] <- "BOLD Low PGS"
  names(table)[4] <- "BOLD High PGS"
  
  # Adding the P values ---------------------------------------------
  
  # T tests
  test_vars <- c("age",
                 "si_1",
                 "raedyrs",
                 "iearn",
                 "no_waves",
                 "female",
                 "smoken",
                 "unins",
                 "max_cv")
  p_2 <- c()
  for (v in test_vars) {
    p_2 <- c(p_2, format(round(perform_ttest(baseline_low, baseline_high, v), digits = 2), nsmall =2))
  }
  p_2 <- c("",p_2[1:5],"",p_2[6:9], round(t.test(low_cess$df$cess, high_cess$df$cess)$p.value, digits = 2),"")
  
  # Add to table
  table <-cbind(table,p_2)
  names(table)[5] <- "BOLD P value"
  
  return(table)
}

descriptive_table_pgs_short <- function(df, name){
  # Create skeleton of table ------------------------------------------------
  # Will later on be filled with values
  table <- tibble(variable = c("",
                               "Age (baseline)",
                               "Smoking PGS",
                               "No. waves present",
                               "",
                               "Female",
                               "Smoking (baseline)",
                               "Persistently uninsured",
                               "CV health shock",
                               "No. of individuals",
                               "Person-year observations"),
                  All = "",
                  lowPGS = "",
                  highPGS = ""
  )
  
  # Create baseline datasets ------------------------------------------------
  # "All"
  baseline_all <- df %>% group_by(hhidpn) %>%
    arrange(year) %>%
    mutate(no_waves = n(),
           max_cv = max(cv)) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    select(hhidpn,
           high_pgs,
           age,
           si_1, 
           no_waves,
           female,
           smoken,
           uninsured,
           max_cv,
           unins)
  
  # "low PGS"
  baseline_low <- baseline_all %>%
    filter(high_pgs == 0)
  
  # "high PGS"
  baseline_high <- baseline_all %>%
    filter(high_pgs == 1)
  
  #All sample for n of person year individuals
  dfAll <- df
  dfLow <- df %>% filter(high_pgs == 0)
  dfHigh <- df %>% filter(high_pgs == 1)
  
  # Fill columns ------------------------------------------------------------
  for(i in 2:4){
    # In order of columns
    groups <- list(NA, baseline_all, baseline_low, baseline_high)
    groupsdf <- list(NA, dfAll, dfLow, dfHigh)
    # Select data for this column
    baseline_df <- groups[[i]]
    all_df <- groupsdf[[i]]
    
    
    # Mean and standard deviations --------------------------------------------
    # Age (baseline)
    sum_age <- paste_mean_std(baseline_df, "age")
    
    # PGS
    sum_pgs <- paste_mean_std(baseline_df, "si_1") #si_1
    
    # No. waves present
    sum_no_waves <- paste_mean_std(baseline_df, "no_waves")
    
    
    # Mean in percent ---------------------------------------------------------
    # Female
    sum_female <- mean_perc(baseline_df, "female")
    
    # Smoking (baseline)
    sum_smoking <- mean_perc(baseline_df, "smoken")
    
    # Persistently uninsured
    sum_unins <- mean_perc(baseline_df, "unins")
    
    # CV health shock
    sum_cv <- mean_perc(baseline_df, "max_cv")
    
    # No. of individuals
    n_individuals <- count_individuals(baseline_df, simple = TRUE)
    
    # No. of person - year individuals
    n_personyearpgs <- all_df %>% count()
  
    # Add calculated statistics to table
    # CAUTION: If a variable is added, moved, or deleted in the following line,
    # the row numbers in the "average cessation rate" just below might also needs to be updated
    column <- c("Mean (SD)",
                sum_age,
                sum_pgs,
                sum_no_waves,
                "\\%",
                sum_female,
                sum_smoking,
                sum_unins,
                sum_cv,
                n_individuals,
                n_personyearpgs)
    table[[i]] <- column
  }
  
  # Change column names for table
  names(table)[1] <- "BOLD"
  names(table)[2] <- "BOLD All"
  names(table)[3] <- "BOLD Low PGS"
  names(table)[4] <- "BOLD High PGS"
  
  return(table)
}

# Function: Creates a descriptive table split by new cv -------------------
descriptive_table_cv <- function(df, name){

  # Create skeleton of table ------------------------------------------------
  table <- tibble(variable = c("", "Age (baseline)", "Smoking PGS", "Years of education", "Income (nominal \\$ 1000)", "No. waves present",
                               "", "Female", "Smoking (baseline)", "Persistently uninsured", "CV health shock",
                               "Avg. cessation rate (baseline smokers)", "No. of individuals"),
                  All = "",
                  NewCV = "",
                  NoNewCV = ""
  )
  
  # Create baseline datasets ------------------------------------------------
  # "All"
  baseline_all <- df %>% group_by(hhidpn) %>%
    arrange(year) %>%
    mutate(no_waves = n(),
           max_cv = max(cv)) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    select(hhidpn, age, si_1, no_waves, female,
           smoken, uninsured, max_cv, unins,
           raedyrs, iearn) %>%
    mutate(iearn = iearn/1000)
  
  # "NewCV"
  baseline_cv <- baseline_all %>%
    filter(max_cv == 1)
  
  # "NoNewCV"
  baseline_nocv <- baseline_all %>%
    filter(max_cv == 0)
  
  # Fill columns ------------------------------------------------------------
  for(i in 2:4){
    # In order of columns
    groups <- list(NA, baseline_all, baseline_cv, baseline_nocv)
    # Select data for this column
    baseline_df <- groups[[i]]
    
    # Mean and standard deviations --------------------------------------------
    # Age (baseline)
    sum_age <- paste_mean_std(baseline_df, "age")
    
    # PGS
    sum_pgs <- paste_mean_std(baseline_df, "si_1")
    
    # Years of education
    # This variable contains missing variables!
    sum_edyrs <- paste_mean_std(baseline_df %>%
                                  filter(!is.na(raedyrs)), "raedyrs")
    
    # Individual income (nominal dollars)
    sum_earnings <- paste_mean_std(baseline_df, "iearn")
    
    # No. waves present
    sum_no_waves <- paste_mean_std(baseline_df, "no_waves")
    
    # Mean in percent ---------------------------------------------------------
    # Female
    sum_female <- mean_perc(baseline_df, "female")
    
    # Smoking (baseline)
    sum_smoking <- mean_perc(baseline_df, "smoken")
    
    # Persistently uninsured
    sum_unins <- mean_perc(baseline_df, "unins")
    
    # CV health shock
    sum_cv <- mean_perc(baseline_df, "max_cv")
    
    # No. of individuals
    n_individuals <- count_individuals(baseline_df, simple = TRUE)
    
    # Add calculated statistics to table
    # CAUTION: If a variable is added, moved, or deleted in the following line,
    # the row numbers in the "average cessation rate" just below might also needs to be updated
    column <- c("Mean (SD)", sum_age, sum_pgs, sum_edyrs, sum_earnings, sum_no_waves, "\\%",
                sum_female, sum_smoking, sum_unins,
                sum_cv, "", n_individuals)
    table[[i]] <- column
  }
  
  
  # Average cessation rate among baseline smokers -----------------------------------------
  # Calculated and added outside of above loop
  
  # COLUMN "All"
  #table[12, 2] <- mean_cess_baseline_smk(df, baseline_all)
  table$All[12] <- mean_cess_baseline_smk(df, baseline_all)
  
  # COLUMN "NewCV"
  newcv_cess <- mean_cess_baseline_smk(df, baseline_cv, return_df = TRUE)
  table$NewCV[12] <- newcv_cess$rate
  
  
  # COLUMN "NoCV"
  nocv_cess <- mean_cess_baseline_smk(df, baseline_nocv, return_df = TRUE)
  table$NoNewCV[12] <- nocv_cess$rate
  
  
  # Change column names for table
  names(table)[1] <- "BOLD "
  names(table)[2] <- "BOLD All"
  names(table)[3] <- "BOLD New Shock"
  names(table)[4] <- "BOLD No new Shock"
  
  # Adding the P values ---------------------------------------------
  
  # T tests
  test_vars <- c("age", "si_1", "raedyrs", "iearn", "no_waves", "female",
                 "smoken", "unins", "max_cv")
  p_2 <- c()
  for (v in test_vars) {
    if (v == "max_cv") {
      p_temp <- "-"
    } else {
      p_temp <- format(round(perform_ttest(baseline_cv, baseline_nocv, v), digits = 2), nsmall = 2)
    }
    p_2 <- c(p_2, p_temp)
  }
  p_2 <- c("",p_2[1:5],"",p_2[6:9], round(t.test(newcv_cess$df$cess, nocv_cess$df$cess)$p.value, digits = 2),"")
  
  # Add to table
  table <-cbind(table,p_2)
  names(table)[5] <- "BOLD P value"
  
  return(table)
}

# Function: Creates a descriptive table split by subgroups ----------------
descriptive_table_subgroups <- function(df,  name, pgs = "si_1"){

  # Create skeleton of table ------------------------------------------------
  # Will later on be filled with values
  if (pgs == "si_1") {
    var_labels <- c("", "Age (baseline)", "Smoking PGS", "Years of education", "Income (nominal \\$ 1000)", "No. waves present",
                  "", "Female", "Smoking (baseline)", "Persistently uninsured",
                  "Avg. cessation rate (baseline smokers)",
                  "No. of individuals","No. of Person-year individuals")
    outcome <- "smoken"
  } else if (pgs == "dpw_1") {
    var_labels <- c("", "Age (baseline)", "Drinking PGS", "Years of education", "Income (nominal \\$ 1000)", "No. waves present",
                    "", "Female", "Drinking (baseline)", "Persistently uninsured",
                    "Avg. cessation rate (baseline drinkers)",
                    "No. of individuals","Person-year individuals")
    outcome <- "drink"
  }
  # } else if (pgs = "bmi_1") {
  #   var_labels <- c("", "Age (baseline)", "BMI PGS", "Years of education", "Income (nominal \\$ 1000)", "No. waves present",
  #                   "", "Female", "Drinking (baseline)", "Persistently uninsured",
  #                   "Avg. cessation rate (baseline drinkers)",
  #                   "No. of individuals")
  # }
  
  table <- tibble(variable = c("Shock at ages 60-64", var_labels,
                               "Shock at ages 67-70", var_labels),
                  lowPGS = "",
                  highPGS = "")
  
  # Create baseline datasets ------------------------------------------------
  # "lowPGS"
  baseline_temp <- df %>%
    group_by(hhidpn) %>%
    arrange(year) %>%
    mutate(no_waves = n(),
           max_cv = max(cv),
           age_at_shock = ifelse(cv == 1, age, 0),
           age_at_shock = max(age_at_shock)) %>%
    filter(max_cv == 1, dplyr::row_number() == 1) %>%
    ungroup() %>%
    select(hhidpn, age, age_at_shock, pgs,
           no_waves, female, outcome,
           high_pgs, unins, raedyrs, iearn) %>%
    mutate(iearn = iearn/1000)
  
  baseline_low <- baseline_temp %>%
    filter(high_pgs == 0)
  
  # "highPGS"
  baseline_high <- baseline_temp %>%
    filter(high_pgs == 1)
  
  # Create pre65 and post65 sub-datasets for "lowPGS"
  # Remember that shocks at 65 and 66 should be excluded in
  # the past in dataset
  baseline_low_pre65 <- baseline_low %>% filter(age_at_shock < 65 & age_at_shock > 0)
  baseline_low_post65 <- baseline_low %>% filter(age_at_shock > 65)
  
  # Create pre65 and post65 sub-datasets for "highPGS"
  baseline_high_pre65 <- baseline_high %>% filter(age_at_shock < 65 & age_at_shock > 0)
  baseline_high_post65 <- baseline_high %>% filter(age_at_shock > 65)
  
  #Person year individuals
  pre65High <- df %>% mutate(age_at_shock = ifelse(cv == 1, age, 0)) %>%
                      group_by(hhidpn) %>% arrange(year) %>%
                      mutate(age_at_shock = max(age_at_shock)) %>% 
                      ungroup() %>% 
                      filter(high_pgs == 1 & age_at_shock < 65 & age_at_shock > 0)  
  
  pre65Low <-  df %>% mutate(age_at_shock = ifelse(cv == 1, age, 0)) %>%
                      group_by(hhidpn) %>% arrange(year) %>%
                      mutate(age_at_shock = max(age_at_shock)) %>% 
                      ungroup() %>% 
                      filter(high_pgs == 0 & age_at_shock < 65 & age_at_shock > 0) 
  
  post65High <-  df %>% mutate(age_at_shock = ifelse(cv == 1, age, 0)) %>%
                      group_by(hhidpn) %>% arrange(year) %>%
                      mutate(age_at_shock = max(age_at_shock)) %>% 
                      ungroup() %>% 
                      filter(high_pgs == 1 & age_at_shock > 65)
  
  post65Low <-  df %>% mutate(age_at_shock = ifelse(cv == 1, age, 0)) %>%
                      group_by(hhidpn) %>% arrange(year) %>%
                      mutate(age_at_shock = max(age_at_shock)) %>% 
                      ungroup() %>% 
                      filter(high_pgs == 0 & age_at_shock > 65) 
    
  # Fill columns ------------------------------------------------------------
  for(i in 2:3){
    # Select correct data for column
    if (i == 2){
      baseline_df_pre65 <- baseline_low_pre65
      baseline_df_post65 <- baseline_low_post65
      py_df_pre65 <- pre65Low
      py_df_post65 <- post65Low
    }
    
    if (i == 3){
      baseline_df_pre65 <- baseline_high_pre65
      baseline_df_post65 <- baseline_high_post65
      py_df_pre65 <- pre65High
      py_df_post65 <- post65High
    }
    
    # Mean and standard deviations --------------------------------------------
    # Age (baseline)
    sum_age_pre65 <- paste_mean_std(baseline_df_pre65, "age")
    sum_age_post65 <- paste_mean_std(baseline_df_post65, "age")
    
    # PGS
    sum_pgs_pre65 <- paste_mean_std(baseline_df_pre65 %>%
                                      filter(!is.na(get(pgs))), pgs)
    sum_pgs_post65 <- paste_mean_std(baseline_df_post65 %>%
                                       filter(!is.na(get(pgs))), pgs)
    
    # Years of education
    # This variable contains missing variables!
    sum_edyrs_pre65 <- paste_mean_std(baseline_df_pre65 %>%
                                        filter(!is.na(raedyrs)), "raedyrs")
    sum_edyrs_post65 <- paste_mean_std(baseline_df_post65 %>%
                                         filter(!is.na(raedyrs)), "raedyrs")
    
    # Income (dollar nominal)
    sum_earnings_pre65 <- paste_mean_std(baseline_df_pre65, "iearn")
    sum_earnings_post65 <- paste_mean_std(baseline_df_post65, "iearn")
    
    # No. waves present
    sum_no_waves_pre65 <- paste_mean_std(baseline_df_pre65, "no_waves")
    sum_no_waves_post65 <- paste_mean_std(baseline_df_post65, "no_waves")
    
    # Mean in percent ---------------------------------------------------------
    # Female
    sum_female_pre65 <- mean_perc(baseline_df_pre65, "female")
    sum_female_post65 <- mean_perc(baseline_df_post65, "female")
    
    # Smoking (baseline)
    sum_smoken_pre65 <- mean_perc(baseline_df_pre65, outcome)
    sum_smoken_post65 <- mean_perc(baseline_df_post65, outcome)
    
    # Persistently uninsured
    sum_mean_unins_pre65 <- mean_perc(baseline_df_pre65, "unins")
    sum_mean_unins_post65 <- mean_perc(baseline_df_post65, "unins")
    
    # No. of individuals
    count_n_pre65 <- count_individuals(baseline_df_pre65, simple = TRUE)
    count_n_post65 <- count_individuals(baseline_df_post65, simple = TRUE)
    
    # No. of  person year individuals
    count_nyear_pre65 <- py_df_pre65 %>%  select(hhidpn) %>% count() %>% as.character()
    count_nyear_post65 <- py_df_post65 %>% select(hhidpn) %>% count() %>% as.character()
    
    
    # Add calculated statistics to table
    # CAUTION: If a variable is added, moved, or deleted in the following line,
    # the row numbers in the "average cessation rate"
    # just below might also needs to be updated
    column <- c("", "Mean (SD)", sum_age_pre65,
                sum_pgs_pre65, sum_edyrs_pre65, sum_earnings_pre65, sum_no_waves_pre65, "\\%",
                sum_female_pre65, sum_smoken_pre65,
                sum_mean_unins_pre65, "", count_n_pre65,count_nyear_pre65,
                "", "Mean (SD)", sum_age_post65,
                sum_pgs_post65, sum_edyrs_post65, sum_earnings_post65, sum_no_waves_post65, "\\%",
                sum_female_post65, sum_smoken_post65,
                sum_mean_unins_post65, "", count_n_post65,count_nyear_post65)
    table[[i]] <- column
  }
  
  # Average cessation rate among baseline smokers ---------------------------------------
  # Calculated and added outside of above loop
  # LOW PRE
  low_pre_cess <- mean_cess_baseline_smk(df, baseline_low_pre65, return_df = TRUE, pgs = pgs)
  table$lowPGS[12] <-low_pre_cess$rate 
 
  # LOW POST
  low_post_cess <- mean_cess_baseline_smk(df, baseline_low_post65, return_df = TRUE, pgs = pgs)
  table$lowPGS[26] <-low_post_cess$rate
  
  # HIGH PRE
  high_pre_cess <- mean_cess_baseline_smk(df, baseline_high_pre65, return_df = TRUE, pgs = pgs)
  table$highPGS[12] <-high_pre_cess$rate
  
  # HIGH POST
  high_post_cess <- mean_cess_baseline_smk(df, baseline_high_post65, return_df = TRUE, pgs = pgs)
  table$highPGS[26] <- high_post_cess$rate
  
  # Change column names for tables
  names(table)[1] <- "BOLD "
  names(table)[2] <- "BOLD Low PGS"
  names(table)[3] <- "BOLD High PGS"
  
  # Adding the P values ---------------------------------------------
  
  # T tests between shock pre-65 & post-65!! What do we want to have here?
  test_vars <- c("age", pgs, "raedyrs", "iearn", "no_waves", "female",
                 outcome, "unins")
  
  p_2 <- c()
  for (v in test_vars) {
    p_2 <- c(p_2, format(round(perform_ttest(baseline_low_pre65, baseline_high_pre65, v), digits = 2), nsmall = 2))
  }
  p_2 <- c("","",p_2[1:5],"",p_2[6:8], format(round(t.test(low_pre_cess$df$cess, high_pre_cess$df$cess)$p.value, digits = 2),nsmall = 2),"","")
  
  p_3 <- c()
  for (v in test_vars) {
    p_3 <- c(p_3, format(round(perform_ttest(baseline_low_post65, baseline_high_post65, v), digits = 2), nsmall = 2))
  }
  p_3 <- c(p_2,"","",p_3[1:5],"",p_3[6:8], format(round(t.test(low_post_cess$df$cess, high_post_cess$df$cess)$p.value, digits = 2),nsmall = 2),"","")
  
  # Add to table
  table <-cbind(table,p_3)
  names(table)[4] <- "BOLD P value"
  
  return(table)
}

# Function: Regression with high_pgs dummy --------------------------------
# PGS3 is TRUE if pgs is divided into three parts, so the respective regression with two PGS dummies is implemented, however, it is FALSE by default
reg_hipgs <- function(df, outcome, shock_var, PGS3 = FALSE, name) {
  # df <- df5575_100
  # outcome <- "smoken"
  # shock_var <- "cv"
  # name  <-  f_identifier
  # PGS3 <-  FALSE

  df_reg <- df
  
  
  df_reg$shock <- df_reg[[shock_var]]
  df_reg$outcome <- df_reg[[outcome]]
  
  if (PGS3 == TRUE) {
    reg_shock <- plm(outcome ~ shock +
                       post65 +
                       shock:post65 +
                       shock:unins +
                       post65:unins +
                       mid_pgs:shock +
                       mid_pgs:post65 +
                       high_pgs:shock +
                       high_pgs:post65 +
                       shock:post65:unins +
                       mid_pgs:shock:unins +
                       mid_pgs:shock:post65 +
                       mid_pgs:post65:unins +
                       high_pgs:shock:unins +
                       high_pgs:shock:post65 +
                       high_pgs:post65:unins +
                       shock:post65:unins:high_pgs +
                       shock:post65:unins:mid_pgs +
                       agem + agem2 + agem3,
                     df_reg,
                     na.action=na.omit,
                     index=c("hhidpn", "year"),
                     model="within",
                     effect="twoway")
    
    # Specify robust or normal standard errors
    # reg_shock_se <- coeftest(reg_shock, vcov. = vcovHC(reg_shock, cluster = "group"))
    # reg_shock_se <- reg_shock
    return(reg_shock)
  } else {
    reg_shock <- plm(outcome ~ shock +
                       post65 +
                       shock:post65 +
                       shock:unins +
                       post65:unins +
                       high_pgs:shock +
                       high_pgs:post65 +
                       shock:post65:unins +
                       high_pgs:shock:unins +
                       high_pgs:shock:post65 +
                       high_pgs:post65:unins +
                       shock:post65:unins:high_pgs +
                       agem + agem2 + agem3,
                     df_reg,
                     na.action= na.omit,
                     index=c("hhidpn", "year"),
                     model="within",
                     effect="twoway")
    
    # Specify robust or normal standard errors
    # reg_shock_se <- coeftest(reg_shock, vcov. = vcovHC(reg_shock, cluster = "group"))
    # reg_shock_se <- reg_shock
    
    # cov <- vcovHC(reg_shock, cluster = "group")
    # stderror <- sqrt(diag(cov))
    
    return(reg_shock)
  }
}

# # Function: Regression to replicate Marti & Richards (2017) ---------------
# reg_replication <- function(df_all, df_pgs, shock_var, name) {
#
#   df_reg_all <- df_all
#   df_reg_all$shock <- df_reg_all[[shock_var]]
#
#   reg_replication <- plm(smoken ~ shock +
#                      post65 +
#                      shock:post65 +
#                      shock:unins +
#                      post65:unins +
#                      shock:post65:unins +
#                      age,
#                    df_reg_all,
#                    na.action=na.omit,
#                    index=c("hhidpn", "year"),
#                    model="within",
#                    effect="twoway")
#
#   # Specify robust or normal standard errors
#   # reg_shock_se <- coeftest(reg_shock, vcov. = vcovHC(reg_shock, cluster = "group"))
#   reg_replication_se <- reg_replication
#
#   cov_replication <- vcovHC(reg_replication, cluster = "group")
#   stderror_replication <- sqrt(diag(cov_replication))
#
#   # Run my PGS regression
#   df_reg_pgs <- df_pgs
#   df_reg_pgs$shock <- df_reg_pgs[[shock_var]]
#
#   reg_pgs <- plm(smoken ~ shock +
#                      post65 +
#                      shock:post65 +
#                      shock:unins +
#                      post65:unins +
#                      high_pgs:shock +
#                      high_pgs:post65 +
#                      shock:post65:unins +
#                      high_pgs:shock:unins +
#                      high_pgs:shock:post65 +
#                      high_pgs:post65:unins +
#                      shock:post65:unins:high_pgs +
#                      age,
#                    df_reg_pgs,
#                    na.action=na.omit,
#                    index=c("hhidpn", "year"),
#                    model="within",
#                    effect="twoway")
#
#   # Specify robust or normal standard errors
#   # reg_shock_se <- coeftest(reg_shock, vcov. = vcovHC(reg_shock, cluster = "group"))
#   reg_pgs_se <- reg_pgs
#
#   cov_pgs <- vcovHC(reg_pgs, cluster = "group")
#   stderror_pgs <- sqrt(diag(cov_pgs))
# }

# Function: Returns a custom ggplot theme ---------------------------------
mytheme <- function(){
  theme(
    text = element_text(
      # family = "Calibri",
      color = "gray25",
      size = 15),
    plot.background =
      element_rect(fill = "white"),
    panel.background = element_rect(fill = "gray87"),
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.text.x = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.background = element_rect(fill = "transparent"),
    #legend.text = element_text(size = 12),
    legend.key = element_blank(),
    panel.grid.major.x = element_blank()
  )
}

# Function: Get significance stars ----------------------------------------
signif_stars <- function(p_value) {
  if (p_value < 0.01) {
    stars <- "***"
  } else if (p_value < 0.05) {
    stars <- "**"
  } else if (p_value < 0.1) {
    stars <- "*"
  } else {
    stars <- ""
  }
  return(stars)
}

# Function: Calculates and plots shock effects for the uninsured----------------------------
# 3rd input is the y achsis of the plot, per default it will be change in smoking
# 4th argument is the amount of PGS dummies: If PGS3 = TRUE, we inlcuded 2 dummies, by default PGS3 = FALSE, so we included 1 dummy
# IMPORTANT: If PGS3 = TRUE, x_achsis label required! (default only for binary PGS variable)
# percent: True if coefficient can be interpreted as percentage
# nb: is the number of breaks in the plot (hogher nb givers more numbers un the y axis). Default is 10
# scale: defalut is TRUE. Allows to modify y axis scale
shock_effects <- function(regression, name, percent = TRUE, scale = TRUE, nb=10, y_achsis = "%-point-change in smoking probability", x_achsis = c("Low PGS", "High PGS"), PGS3 = FALSE){

  # Effect of a health shock for low-PGS people BEFORE 65
  glht_low_pre65 <- summary(glht(regression,
                                 linfct = c("shock + shock:unins = 0"),
                                 vcov = vcovHC(regression, cluster = "group")))
  if (PGS3 == TRUE) {
    # Effect of a health shock for mid-PGS people BEFORE 65
    glht_mid_pre65 <-summary(glht(regression,
                                  linfct = c("shock + shock:unins +
                                            shock:mid_pgs +
                                            shock:unins:mid_pgs = 0"),
                                  vcov = vcovHC(regression, cluster = "group")))
  }
  # Effect of a health shock for high-PGS people BEFORE 65
  glht_high_pre65 <-summary(glht(regression,
                                 linfct = c("shock + shock:unins +
                                            shock:high_pgs +
                                            shock:unins:high_pgs = 0"),
                                 vcov = vcovHC(regression, cluster = "group")))
  
  # Effect of a health shock for low-PGS people AFTER 65
  glht_low_post65 <-summary(glht(regression,
                                 linfct = c("shock + shock:post65 + shock:unins +
                                             shock:post65:unins = 0"),
                                 vcov = vcovHC(regression, cluster = "group")))
  if (PGS3 == TRUE) {
    # Effect of a health shock for mid-PGS people AFTER 65
    glht_mid_post65 <- summary(glht(regression,
                                    linfct = c("shock +
                                              shock:post65 +
                                              shock:unins +
                                              shock:mid_pgs +
                                              shock:post65:unins +
                                              shock:unins:mid_pgs +
                                              shock:post65:mid_pgs +
                                              shock:post65:unins:mid_pgs = 0"),
                                    vcov = vcovHC(regression, cluster = "group")))
  }
  
  # Effect of a health shock for high-PGS people AFTER 65
  glht_high_post65 <- summary(glht(regression,
                                   linfct = c("shock +
								               shock:post65 +
      											   shock:unins +
		      									   shock:high_pgs +
								      			   shock:post65:unins +
											         shock:unins:high_pgs +
			      								   shock:post65:high_pgs +
									      		   shock:post65:unins:high_pgs = 0"),
                                   vcov = vcovHC(regression, cluster = "group")))
  
  
  # Difference in the effect of a health shock for low-PGS people after 65 minus before 65
  ####Check: it should be the same as glht_low_post65 - glht_low_pre65
  glht_prelow_vs_postlow <- summary(glht(regression,
                                         linfct = c("shock:post65 +
                                                    shock:post65:unins = 0"),
                                         vcov = vcovHC(regression, cluster = "group")))
  
  # Difference in the effect of a health shock for mid-PGS people before and after 65
  ####Check: it should be the same as glht_mid_post65 - glht_mid_pre65
  if (PGS3 == TRUE) {
    glht_premid_vs_postmid <- summary(glht(regression,
                                           linfct = c("shock:post65 +
                                                      shock:post65:unins +
                                                      shock:post65:mid_pgs +
                                                      shock:post65:unins:mid_pgs = 0"),
                                           vcov = vcovHC(regression, cluster = "group")))
  }
  
  # Difference in the effect of a health shock for high-PGS people before and after 65
  ####Check: it should be the same as glht_high_post65 - glht_high_pre65
  glht_prehigh_vs_posthigh <- summary(glht(regression,
                                           linfct = c("shock:post65 +
                                                      shock:post65:unins +
                                                      shock:post65:high_pgs +
                                                      shock:post65:unins:high_pgs = 0"),
                                           vcov = vcovHC(regression, cluster = "group")))
  
  # Difference in the effect of a health shock between low- and high-PGS people BEFORE 65
  ####Check: it should be the same as glht_high_pre65 - glht_low_pre65
  glht_prelow_vs_prehigh <- summary(glht(regression,
                                         linfct = c("shock:high_pgs +
                                                    shock:unins:high_pgs = 0"),
                                         vcov = vcovHC(regression, cluster = "group")))
  if (PGS3 == TRUE) {
    # Difference in the effect of a health shock between low- and mid-PGS people BEFORE 65
    ####Check: it should be the same as glht_mid_pre65 - glht_low_pre65
    glht_prelow_vs_premid <- summary(glht(regression,
                                          linfct = c("shock:mid_pgs +
                                                    shock:unins:mid_pgs = 0"),
                                          vcov = vcovHC(regression, cluster = "group")))
    
    # Difference in the effect of a health shock between mid- and high-PGS people BEFORE 65
    ####Check: it should be the same as glht_high_pre65 - glht_mid_pre65
    glht_premid_vs_prehigh <- summary(glht(regression,
                                           linfct = c("shock:high_pgs +
                                                      shock:unins:high_pgs -
                                                      shock:mid_pgs -
                                                      shock:unins:mid_pgs = 0"),
                                           vcov = vcovHC(regression, cluster = "group")))
  }
  
  # Difference in the effect of a health shock between low- and high-PGS people AFTER 65
  ####Check: it should be the same as glht_high_post65 - glht_low_post65
  glht_postlow_vs_posthigh <- summary(glht(regression,
                                           linfct = c("shock:high_pgs +
                                                      post65:high_pgs +
                                                      shock:unins:high_pgs +
                                                      shock:post65:high_pgs +
                                                      post65:unins:high_pgs +
                                                      shock:post65:unins:high_pgs = 0"),
                                           vcov = vcovHC(regression, cluster = "group")))
  if (PGS3 == TRUE) {
    # Difference in the effect of a health shock between low- and mid-PGS people AFTER 65
    ####Check: it should be the same as glht_mid_post65 - glht_low_post65
    glht_postlow_vs_postmid <- summary(glht(regression,
                                            linfct = c("shock:mid_pgs +
                                                        post65:mid_pgs +
                                                        shock:unins:mid_pgs +
                                                        shock:post65:mid_pgs +
                                                        post65:unins:mid_pgs +
                                                        shock:post65:unins:mid_pgs = 0"),
                                            vcov = vcovHC(regression, cluster = "group")))
    
    # Difference in the effect of a health shock between low- and high-PGS people AFTER 65
    ####Check: it should be the same as glht_high_post65 - glht_low_post65
    glht_postmid_vs_posthigh <- summary(glht(regression,
                                             linfct = c("shock:high_pgs +
                                                        post65:high_pgs +
                                                        shock:unins:high_pgs +
                                                        shock:post65:high_pgs +
                                                        post65:unins:high_pgs +
                                                        shock:post65:unins:high_pgs -
                                                        shock:mid_pgs -
                                                        post65:mid_pgs -
                                                        shock:unins:mid_pgs -
                                                        shock:post65:mid_pgs -
                                                        post65:unins:mid_pgs -
                                                        shock:post65:unins:mid_pgs = 0"),
                                             vcov = vcovHC(regression, cluster = "group")))
    
  }
  
  # Difference in the effect of health insurance status on the effect of a health shock
  # between high- and low-PGS people
  ####Check: this should be the triple difference
  glht_slopehigh_vs_slopelow <- summary(glht(regression,
                                             linfct = c("shock:post65:high_pgs +
                                                        shock:post65:unins:high_pgs = 0"),
                                             vcov = vcovHC(regression, cluster = "group")))
  
  if (PGS3 == TRUE) {
    # Difference in the effect of health insurance status on the effect of a health shock
    # between mid- and low-PGS people
    ####Check: this should be the triple difference
    glht_slopemid_vs_slopelow <- summary(glht(regression,
                                              linfct = c("shock:post65:mid_pgs +
                                                          shock:post65:unins:mid_pgs = 0"),
                                              vcov = vcovHC(regression, cluster = "group")))
    
    # Difference in the effect of health insurance status on the effect of a health shock
    # between high- and mid-PGS people
    ####Check: this should be the triple difference
    glht_slopehigh_vs_slopemid <- summary(glht(regression,
                                               linfct = c("shock:post65:high_pgs +
                                                          shock:post65:unins:high_pgs -
                                                          shock:post65:mid_pgs -
                                                          shock:post65:unins:mid_pgs = 0"),
                                               vcov = vcovHC(regression, cluster = "group")))
  }
  
  # Make graph of results
  
  if (PGS3 == FALSE) {
    df_plot <- tibble(time = c(rep("pre65", 2),rep("post65", 2)),
                    group = c(rep(c("lowPGS", "highPGS"), 2)),
                    coeff = c(glht_low_pre65$test$coef,
                              glht_high_pre65$test$coef,
                              glht_low_post65$test$coef,
                              glht_high_post65$test$coef),
                    se    = c(glht_low_pre65$test$sigma,
                              glht_high_pre65$test$sigma,
                              glht_low_post65$test$sigma,
                              glht_high_post65$test$sigma
                    ))
    
    df_plot$time <- as.factor(df_plot$time)
    df_plot$time <- factor(df_plot$time, levels = c("pre65", "post65"))
    df_plot$group <- factor(df_plot$group, levels = c("lowPGS", "highPGS"))
    
    # if beta cannot interpreted as percentage, we have to divide by 100 as it will be mulitplied later
    if (percent==FALSE){
      df_plot <- df_plot %>% mutate(coeff = coeff / 100, se = se / 100)
    }
    
    minx <- round(min(df_plot$coeff*100-1.96*df_plot$se*100))
    maxx <- round(max(df_plot$coeff*100+1.96*df_plot$se*100))
    
    plot <- ggplot(df_plot %>% mutate(coeff = coeff * 100, se = se * 100), aes(x=group, y=coeff, color=time)) +
      geom_point(position = position_dodge(width = 0.1), size = 3) +
      geom_errorbar(aes(ymin = coeff-1.96*se, ymax = coeff+1.96*se),
                    width=.1,
                    position = position_dodge(width = 0.1),
                    size = 1) +
      labs(x = NULL,
           y = y_achsis) +
      scale_x_discrete(labels = x_achsis)  +
      scale_color_manual(labels = c("Pre-65", "Post-65"),
                         name = "Shock timing",
                         values=c("#F36A6A", "#5CACEE")) +
      geom_hline(aes(yintercept = 0), color = "gray40", linetype = "dashed") +
      #ylim(-0.5,0.3) +
      theme_minimal() +
      theme(legend.position = c(0.88, 0.1),
            axis.ticks.x = element_blank())
    
    if(scale == TRUE){
      plot <- plot +  scale_y_continuous(n.breaks = nb)}
    
    ggsave(paste0( name, ".png"),
           plot = plot, width = 7, height = 7)
  } else {
    df_plot <- tibble(time = c(rep("pre65", 3),rep("post65", 3)),
                      group = c(rep(c("lowPGS", "middlePGS", "highPGS"), 2)),
                      coeff = c(glht_low_pre65$test$coef,
                                glht_mid_pre65$test$coef,
                                glht_high_pre65$test$coef,
                                glht_low_post65$test$coef,
                                glht_mid_post65$test$coef,
                                glht_high_post65$test$coef),
                      se    = c(glht_low_pre65$test$sigma,
                                glht_mid_pre65$test$sigma,
                                glht_high_pre65$test$sigma,
                                glht_low_post65$test$sigma,
                                glht_mid_post65$test$sigma,
                                glht_high_post65$test$sigma
                      ))
    
    df_plot$time <- as.factor(df_plot$time)
    df_plot$time <- factor(df_plot$time, levels = c("pre65", "post65"))
    df_plot$group <- factor(df_plot$group, levels = c("lowPGS","middlePGS", "highPGS"))
    
    minx <- round(min(df_plot$coeff*100-1.96*df_plot$se*100))
    maxx <- round(max(df_plot$coeff*100+1.96*df_plot$se*100))
    
    plot <- ggplot(df_plot %>% mutate(coeff = coeff * 100, se = se * 100), aes(x=group, y=coeff, color=time)) +
      geom_point(position = position_dodge(width = 0.1), size = 3) +
      geom_errorbar(aes(ymin = coeff-1.96*se, ymax = coeff+1.96*se),
                    width=.1,
                    position = position_dodge(width = 0.1),
                    size = 1) +
      labs(x = NULL,
           y = y_achsis) +
      scale_x_discrete(labels = x_achsis)  +
      scale_color_manual(labels = c("Pre-65", "Post-65"),
                         name = "Shock timing",
                         values=c("#F36A6A", "#5CACEE")) +
      geom_hline(aes(yintercept = 0), color = "gray40", linetype = "dashed") +
      #ylim(-0.5,0.3) +
      theme_minimal() +
      theme(legend.position = c(0.88, 0.1),
            axis.ticks.x = element_blank())
    
    if(scale == TRUE){
      plot <- plot +  scale_y_continuous(n.breaks = nb)}
    
    ggsave(paste0(name, ".png"),
           plot = plot, width = 7, height = 7)
  }
  
  # Make an overview table for the main paper
  if (PGS3 == FALSE) {
    table <- tibble(HI = c("","", "Pre 65", "", "Post 65", "", "", "Post 65 - Pre 65", "","","Post 65 - Pre 65", ""),
                    lowPGS = c(""),
                    highPGS = c(""))
    
    #table[1, 2] <- "BOLD Effect of health shock on smoking probability"
    
    table[2, 1] <- ""
    table[2, 2] <- "Low PGS"
    table[2, 3] <- "High PGS"
    
    get_coef_str <- function(lin_hyp) {
      return(paste0(round(lin_hyp$test$coef, digits = 3), signif_stars(lin_hyp$test$pvalues[1])))
    }
    
    table[3, 2] <- get_coef_str(glht_low_pre65)
    table[4, 2] <- paste0("(", round(glht_low_pre65$test$sigma, digits = 3), ")")
    table[3, 3] <- get_coef_str(glht_high_pre65)
    table[4, 3] <- paste0("(", round(glht_high_pre65$test$sigma, digits = 3), ")")
    
    table[5, 2] <- get_coef_str(glht_low_post65)
    table[6, 2] <- paste0("(", round(glht_low_post65$test$sigma, digits = 3), ")")
    table[5, 3] <- get_coef_str(glht_high_post65)
    table[6, 3] <- paste0("(", round(glht_high_post65$test$sigma, digits = 3), ")")
    
    #table[7, 2] <- "BOLD Effect of health insurance on effect of health shock"
    table[7, 2] <- "Low PGS"
    table[7, 3] <- "High PGS"
    
    table[8, 2] <- get_coef_str(glht_prelow_vs_postlow)
    table[9, 2] <- paste0("(", round(glht_prelow_vs_postlow$test$sigma, digits = 3), ")")
    table[8, 3] <- get_coef_str(glht_prehigh_vs_posthigh)
    table[9, 3] <- paste0("(", round(glht_prehigh_vs_posthigh$test$sigma, digits = 3), ")")
    
    #table[11, 2] <- "BOLD Differential effect of health insurance by genetic group"
    table[10, 2] <- "High PGS "
    table[10, 3] <- "- low PGS"
    
    table[11, 2] <- get_coef_str(glht_slopehigh_vs_slopelow)
    table[12, 2] <- paste0("(", round(glht_slopehigh_vs_slopelow$test$sigma, digits = 3), ")")
    
    names(table)[1] <- ""
    names(table)[2] <- "Low PGS"
    names(table)[3] <- "High PGS"
    
    # Return
    results_list = list("plot" = plot, "overview_table" = table)
    return(results_list)
    
  } else {
    table <- tibble(HI = c( "","", "Pre 65", "", "Post 65", "", "", "Post 65 - Pre 65", "", "", "Post 65 - Pre 65", ""),
                    lowPGS = c(""),
                    midPGS = c(""),
                    highPGS = c(""))
    
    #table[1, 2] <- "BOLD Effect of health shock on smoking probability"
    table[2, 2] <- "Low PGS"
    table[2, 3] <- "Middle PGS"
    table[2, 4] <- "High PGS"
    
    get_coef_str <- function(lin_hyp) {
      return(paste0(round(lin_hyp$test$coef, digits = 3), signif_stars(lin_hyp$test$pvalues[1])))
    }
    
    table[3, 2] <- get_coef_str(glht_low_pre65)
    table[4, 2] <- paste0("(", round(glht_low_pre65$test$sigma, digits = 3), ")")
    table[3, 3] <- get_coef_str(glht_mid_pre65)
    table[4, 3] <- paste0("(", round(glht_mid_pre65$test$sigma, digits = 3), ")")
    table[3, 4] <- get_coef_str(glht_high_pre65)
    table[4, 4] <- paste0("(", round(glht_high_pre65$test$sigma, digits = 3), ")")
    
    table[5, 2] <- get_coef_str(glht_low_post65)
    table[6, 2] <- paste0("(", round(glht_low_post65$test$sigma, digits = 3), ")")
    table[5, 3] <- get_coef_str(glht_mid_post65)
    table[6, 3] <- paste0("(", round(glht_mid_post65$test$sigma, digits = 3), ")")
    table[5, 4] <- get_coef_str(glht_high_post65)
    table[6, 4] <- paste0("(", round(glht_high_post65$test$sigma, digits = 3), ")")
    
    #table[7, 2] <- "BOLD Effect of health insurance on effect of health shock"
    table[7, 2] <- "Low PGS"
    table[7, 3] <- "Middle PGS"
    table[7, 4] <- "High PGS"
    
    table[8, 2] <- get_coef_str(glht_prelow_vs_postlow)
    table[9, 2] <- paste0("(", round(glht_prelow_vs_postlow$test$sigma, digits = 3), ")")
    table[8, 3] <- get_coef_str(glht_premid_vs_postmid)
    table[9, 3] <- paste0("(", round(glht_premid_vs_postmid$test$sigma, digits = 3), ")")
    table[8, 4] <- get_coef_str(glht_prehigh_vs_posthigh)
    table[9, 4] <- paste0("(", round(glht_prehigh_vs_posthigh$test$sigma, digits = 3), ")")
    
    #table[11, 2] <- "BOLD Differential effect of health insurance by genetic group"
    table[10, 2] <- "High PGS - low PGS"
    table[10, 3] <- "Middle PGS - low PGS"
    table[10, 4] <- "High PGS - middle PGS"
    
    
    table[11, 2] <- get_coef_str(glht_slopehigh_vs_slopelow)
    table[12, 2] <- paste0("(", round(glht_slopehigh_vs_slopelow$test$sigma, digits = 3), ")")
    table[11, 3] <- get_coef_str(glht_slopemid_vs_slopelow)
    table[12, 3] <- paste0("(", round(glht_slopemid_vs_slopelow$test$sigma, digits = 3), ")")
    table[11, 4] <- get_coef_str(glht_slopehigh_vs_slopemid)
    table[12, 4] <- paste0("(", round(glht_slopehigh_vs_slopemid$test$sigma, digits = 3), ")")
    
    names(table)[1] <- ""
    names(table)[2] <- ""
    names(table)[3] <- ""
    names(table)[4] <- ""
    
    # Return
    results_list = list("plot" = plot, "overview_table" = table)
    return(results_list)
    
  }
}

# # For analysis where we run linear hypothesis in stata
# stata_shock_effects <- function(df_coef, name, y_achsis = "%-point-change in smoking probability", x_achsis = c("Low PGS", "High PGS")){
#   df_plot <- tibble(time = c(rep("pre65", 2),rep("post65", 2)),
#                     group = c(rep(c("lowPGS", "highPGS"), 2)),
#                     coeff = c(as.numeric(df_coef[1]),
#                               as.numeric(df_coef[4]),
#                               as.numeric(df_coef[7]),
#                               as.numeric(df_coef[10])),
#                     se    = c(as.numeric(df_coef[2]),
#                               as.numeric(df_coef[5]),
#                               as.numeric(df_coef[8]),
#                               as.numeric(df_coef[11])
#                     ))
# 
# 
#   df_plot$time <- as.factor(df_plot$time)
#   df_plot$time <- factor(df_plot$time, levels = c("pre65", "post65"))
#   df_plot$group <- factor(df_plot$group, levels = c("lowPGS", "highPGS"))
# 
#   plot <- ggplot(df_plot %>% mutate(coeff = coeff * 100, se = se * 100), aes(x=group, y=coeff, color=time)) +
#     geom_point(position = position_dodge(width = 0.1), size = 3) +
#     geom_errorbar(aes(ymin = coeff-1.96*se, ymax = coeff+1.96*se),
#                   width=.1,
#                   position = position_dodge(width = 0.1),
#                   size = 1) +
#     labs(x = NULL,
#          y = y_achsis) +
#     scale_x_discrete(labels = x_achsis)  +
#     scale_color_manual(labels = c("Pre-65", "Post-65"),
#                        name = "Shock timing",
#                        values=c("#F36A6A", "#5CACEE")) +
#     geom_hline(aes(yintercept = 0), color = "gray40", linetype = "dashed") +
#     #ylim(-0.5,0.3) +
#     theme_minimal() +
#     theme(legend.position = c(0.88, 0.1),
#           axis.ticks.x = element_blank())
# 
#   ggsave(paste0("../3_output/shock_effects/", name, "plot.png"),
#          plot = plot, width = 7, height = 7)
# 
#   # Create the table
#   table <- tibble(HI = c("","", "Pre 65", "", "Post 65", "", "", "Post 65 - Pre 65", "","","Post 65 - Pre 65", ""),
#                   lowPGS = c(""),
#                   highPGS = c(""))
# 
#   #table[1, 2] <- "BOLD Effect of health shock on smoking probability"
# 
#   table[2, 1] <- ""
#   table[2, 2] <- "Low PGS"
#   table[2, 3] <- "High PGS"
# 
#   #get_coef_str <- function(lin_hyp) {
#   #  return(paste0(round(lin_hyp$test$coef, digits = 3), signif_stars(lin_hyp$test$pvalues[1])))
#   #}
# 
#   table[3, 2] <- paste0(round(as.numeric(df_coef[1]), digits = 3), signif_stars(pt(as.numeric(df_coef[1])/as.numeric(df_coef[2]), df = as.numeric(df_coef[3]))))
#   table[4, 2] <- round(as.numeric(df_coef[2]), digits = 3)
#   table[3, 3] <- paste0(round(as.numeric(df_coef[4]), digits = 3), signif_stars(pt(as.numeric(df_coef[4])/as.numeric(df_coef[5]), df = as.numeric(df_coef[6]))))
#   table[4, 3] <- round(as.numeric(df_coef[5]), digits = 3)
# 
#   table[5, 2] <- paste0(round(as.numeric(df_coef[7]), digits = 3), signif_stars(pt(as.numeric(df_coef[7])/as.numeric(df_coef[8]), df = as.numeric(df_coef[9]))))
#   table[6, 2] <- round(as.numeric(df_coef[8]), digits = 3)
#   table[5, 3] <- paste0(round(as.numeric(df_coef[10]), digits = 3), signif_stars(pt(as.numeric(df_coef[10])/as.numeric(df_coef[11]), df = as.numeric(df_coef[12]))))
#   table[6, 3] <- round(as.numeric(df_coef[11]), digits = 3)
# 
#   #table[7, 2] <- "BOLD Effect of health insurance on effect of health shock"
#   table[7, 2] <- "Low PGS"
#   table[7, 3] <- "High PGS"
# 
#   table[8, 2] <- paste0(round(as.numeric(df_coef[13]), digits = 3), signif_stars(pt(as.numeric(df_coef[13])/as.numeric(df_coef[14]), df = as.numeric(df_coef[15]))))
#   table[9, 2] <- round(as.numeric(df_coef[14]), digits = 3)
#   table[8, 3] <- paste0(round(as.numeric(df_coef[16]), digits = 3), signif_stars(pt(as.numeric(df_coef[16])/as.numeric(df_coef[17]), df = as.numeric(df_coef[18]))))
#   table[9, 3] <- round(as.numeric(df_coef[17]), digits = 3)
# 
#   #table[11, 2] <- "BOLD Differential effect of health insurance by genetic group"
#   table[10, 2] <- "High PGS "
#   table[10, 3] <- "- low PGS"
# 
#   table[11, 2] <- paste0(round(as.numeric(df_coef[19]), digits = 3), signif_stars(pt(as.numeric(df_coef[19])/as.numeric(df_coef[20]), df = as.numeric(df_coef[21]))))
#   table[12, 2] <- round(as.numeric(df_coef[20]), digits = 3)
# 
#   names(table)[1] <- ""
#   names(table)[2] <- "Low PGS"
#   names(table)[3] <- "High PGS"
# 
#   # Return
#   results_list = list("plot" = plot, "overview_table" = table)
#   return(results_list)
#   }

# Function: Calculates and plots shock effects FOR THE INSURED----------------------------
# 3rd input is the y achsis of the plot, per default it will be change in smoking
# 4th argument is the amount of PGS dummies: If PGS3 = TRUE, we inlcuded 2 dummies, by default PGS3 = FALSE, so we included 1 dummy
# IMPORTANT: If PGS3 = TRUE, x_achsis label required! (default only for binary PGS variable)
# percent: True if coefficient can be interpreted as percentage
shock_effects_ins <- function(regression, name, percent = TRUE, nb = 10,  y_achsis = "%-point-change in smoking probability", x_achsis = c("Low PGS", "High PGS"), PGS3 = FALSE){

    glht_low_pre65 <- summary(glht(regression,
                                 linfct = c("shock = 0"),
                                 vcov = vcov(regression, cluster = "group")))
  if (PGS3 == TRUE) {
    # Effect of a health shock for mid-PGS people BEFORE 65
    glht_mid_pre65 <-summary(glht(regression,
                                  linfct = c("shock +
                                            shock:mid_pgs = 0"),
                                  vcov = vcovHC(regression, cluster = "group")))
  }
  # Effect of a health shock for high-PGS people BEFORE 65
  glht_high_pre65 <-summary(glht(regression,
                                 linfct = c("shock + 
                                            shock:high_pgs = 0"),
                                 vcov = vcovHC(regression, cluster = "group")))
  
  # Effect of a health shock for low-PGS people AFTER 65
  glht_low_post65 <-summary(glht(regression,
                                 linfct = c("shock + shock:post65 = 0"),
                                 vcov = vcovHC(regression, cluster = "group")))
  if (PGS3 == TRUE) {
    # Effect of a health shock for mid-PGS people AFTER 65
    glht_mid_post65 <- summary(glht(regression,
                                    linfct = c("shock +
                                              shock:post65 +
                                              shock:mid_pgs +
                                              shock:post65:mid_pgs = 0"),
                                    vcov = vcovHC(regression, cluster = "group")))
  }
  
  # Effect of a health shock for high-PGS people AFTER 65
  glht_high_post65 <- summary(glht(regression,
                                   linfct = c("shock +
                                               shock:post65 +
                                               shock:high_pgs +
                                               shock:post65:high_pgs = 0"),
                                   vcov = vcovHC(regression, cluster = "group")))
  
  
  # Difference in the effect of a health shock for low-PGS people after 65 minus before 65
  ####Check: it should be the same as glht_low_post65 - glht_low_pre65
  glht_prelow_vs_postlow <- summary(glht(regression,
                                         linfct = c("shock:post65 = 0"),
                                         vcov = vcovHC(regression, cluster = "group")))
  
  # Difference in the effect of a health shock for mid-PGS people before and after 65
  ####Check: it should be the same as glht_mid_post65 - glht_mid_pre65
  if (PGS3 == TRUE) {
    glht_premid_vs_postmid <- summary(glht(regression,
                                           linfct = c("shock:post65 +
                                                      shock:post65:mid_pgs = 0"),
                                           vcov = vcovHC(regression, cluster = "group")))
  }
  
  # Difference in the effect of a health shock for high-PGS people before and after 65
  ####Check: it should be the same as glht_high_post65 - glht_high_pre65
  glht_prehigh_vs_posthigh <- summary(glht(regression,
                                           linfct = c("shock:post65 +
                                                      shock:post65:high_pgs = 0"),
                                           vcov = vcovHC(regression, cluster = "group")))
  
  # Difference in the effect of a health shock between low- and high-PGS people BEFORE 65
  ####Check: it should be the same as glht_high_pre65 - glht_low_pre65
  glht_prelow_vs_prehigh <- summary(glht(regression,
                                         linfct = c("shock:high_pgs = 0"),
                                         vcov = vcovHC(regression, cluster = "group")))
  if (PGS3 == TRUE) {
    # Difference in the effect of a health shock between low- and mid-PGS people BEFORE 65
    ####Check: it should be the same as glht_mid_pre65 - glht_low_pre65
    glht_prelow_vs_premid <- summary(glht(regression,
                                          linfct = c("shock:mid_pgs = 0"),
                                          vcov = vcovHC(regression, cluster = "group")))
    
    # Difference in the effect of a health shock between mid- and high-PGS people BEFORE 65
    ####Check: it should be the same as glht_high_pre65 - glht_mid_pre65
    glht_premid_vs_prehigh <- summary(glht(regression,
                                           linfct = c("shock:high_pgs -
                                                      shock:mid_pgs = 0"),
                                           vcov = vcovHC(regression, cluster = "group")))
  }
  
  # Difference in the effect of a health shock between low- and high-PGS people AFTER 65
  ####Check: it should be the same as glht_high_post65 - glht_low_post65
  glht_postlow_vs_posthigh <- summary(glht(regression,
                                           linfct = c("shock:high_pgs +
                                                      post65:high_pgs +
                                                      shock:post65:high_pgs = 0"),
                                           vcov = vcovHC(regression, cluster = "group")))
  if (PGS3 == TRUE) {
    # Difference in the effect of a health shock between low- and mid-PGS people AFTER 65
    ####Check: it should be the same as glht_mid_post65 - glht_low_post65
    glht_postlow_vs_postmid <- summary(glht(regression,
                                            linfct = c("shock:mid_pgs +
                                                        post65:mid_pgs +
                                                        shock:post65:mid_pgs = 0"),
                                            vcov = vcovHC(regression, cluster = "group")))
    
    # Difference in the effect of a health shock between low- and high-PGS people AFTER 65
    ####Check: it should be the same as glht_high_post65 - glht_low_post65
    glht_postmid_vs_posthigh <- summary(glht(regression,
                                             linfct = c("shock:high_pgs +
                                                        post65:high_pgs +
                                                        shock:post65:high_pgs -
                                                        shock:mid_pgs -
                                                        post65:mid_pgs -
                                                        shock:post65:mid_pgs = 0"),
                                             vcov = vcovHC(regression, cluster = "group")))
    
  }
  
  # Difference in the effect of health insurance status on the effect of a health shock
  # between high- and low-PGS people
  ####Check: this should be the triple difference
  glht_slopehigh_vs_slopelow <- summary(glht(regression,
                                             linfct = c("shock:post65:high_pgs = 0"),
                                             vcov = vcovHC(regression, cluster = "group")))
  
  if (PGS3 == TRUE) {
    # Difference in the effect of health insurance status on the effect of a health shock
    # between mid- and low-PGS people
    ####Check: this should be the triple difference
    glht_slopemid_vs_slopelow <- summary(glht(regression,
                                              linfct = c("shock:post65:mid_pgs = 0"),
                                              vcov = vcovHC(regression, cluster = "group")))
    
    # Difference in the effect of health insurance status on the effect of a health shock
    # between high- and mid-PGS people
    ####Check: this should be the triple difference
    glht_slopehigh_vs_slopemid <- summary(glht(regression,
                                               linfct = c("shock:post65:high_pgs -
                                                          shock:post65:mid_pgs = 0"),
                                               vcov = vcovHC(regression, cluster = "group")))
  }
  
  # Make a graph
  
  if (PGS3 == FALSE) {
    df_plot <- tibble(time = c(rep("pre65", 2),rep("post65", 2)),
                      group = c(rep(c("lowPGS", "highPGS"), 2)),
                      coeff = c(glht_low_pre65$test$coef,
                                glht_high_pre65$test$coef,
                                glht_low_post65$test$coef,
                                glht_high_post65$test$coef),
                      se    = c(glht_low_pre65$test$sigma,
                                glht_high_pre65$test$sigma,
                                glht_low_post65$test$sigma,
                                glht_high_post65$test$sigma
                      ))
    
    df_plot$time <- as.factor(df_plot$time)
    df_plot$time <- factor(df_plot$time, levels = c("pre65", "post65"))
    df_plot$group <- factor(df_plot$group, levels = c("lowPGS", "highPGS"))
    
    # if beta cannot interpreted as percentage, we have to divide by 100 as it will be mulitplied later
    if (percent==FALSE){
      df_plot <- df_plot %>% mutate(coeff = coeff / 100, se = se / 100)
    }

    
    plot <- ggplot(df_plot %>% mutate(coeff = coeff * 100, se = se * 100), aes(x=group, y=coeff, color=time)) +
      geom_point(position = position_dodge(width = 0.1), size = 3) +
      geom_errorbar(aes(ymin = coeff-1.96*se, ymax = coeff+1.96*se),
                    width=.1,
                    position = position_dodge(width = 0.1),
                    size = 1) +
      labs(x = NULL,
           y = y_achsis) +
      scale_x_discrete(labels = x_achsis)  +
      scale_color_manual(labels = c("Pre-65", "Post-65"),
                         name = "Shock timing",
                         values=c("#F36A6A", "#5CACEE")) +
      geom_hline(aes(yintercept = 0), color = "gray40", linetype = "dashed") +
      #ylim(-0.5,0.3) +
      theme_minimal() +
      theme(legend.position = c(0.88, 0.1),
            axis.ticks.x = element_blank())

    plot <- plot +  scale_y_continuous(n.breaks = nb)
    
    ggsave(paste0(name, ".png"),
           plot = plot, width = 7, height = 7)
  } else {
    df_plot <- tibble(time = c(rep("pre65", 3),rep("post65", 3)),
                      group = c(rep(c("lowPGS", "middlePGS", "highPGS"), 2)),
                      coeff = c(glht_low_pre65$test$coef,
                                glht_mid_pre65$test$coef,
                                glht_high_pre65$test$coef,
                                glht_low_post65$test$coef,
                                glht_mid_post65$test$coef,
                                glht_high_post65$test$coef),
                      se    = c(glht_low_pre65$test$sigma,
                                glht_mid_pre65$test$sigma,
                                glht_high_pre65$test$sigma,
                                glht_low_post65$test$sigma,
                                glht_mid_post65$test$sigma,
                                glht_high_post65$test$sigma
                      ))
    
    df_plot$time <- as.factor(df_plot$time)
    df_plot$time <- factor(df_plot$time, levels = c("pre65", "post65"))
    df_plot$group <- factor(df_plot$group, levels = c("lowPGS","middlePGS", "highPGS"))
    

    plot <- ggplot(df_plot %>% mutate(coeff = coeff * 100, se = se * 100), aes(x=group, y=coeff, color=time)) +
      geom_point(position = position_dodge(width = 0.1), size = 3) +
      geom_errorbar(aes(ymin = coeff-1.96*se, ymax = coeff+1.96*se),
                    width=.1,
                    position = position_dodge(width = 0.1),
                    size = 1) +
      labs(x = NULL,
           y = y_achsis) +
      scale_x_discrete(labels = x_achsis)  +
      scale_color_manual(labels = c("Pre-65", "Post-65"),
                         name = "Shock timing",
                         values=c("#F36A6A", "#5CACEE")) +
      geom_hline(aes(yintercept = 0), color = "gray40", linetype = "dashed") +
      #ylim(-0.5,0.3) +
      theme_minimal() +
      theme(legend.position = c(0.88, 0.1),
            axis.ticks.x = element_blank())
    plot <- plot +  scale_y_continuous(n.breaks = nb)
    
    ggsave(paste0(name, ".png"),
           plot = plot, width = 7, height = 7)
  }
  
  # Make an overview table for the main paper
  if (PGS3 == FALSE) {
    table <- tibble(HI = c("","", "Pre 65", "", "Post 65", "", "", "Post 65 - Pre 65", "","","Post 65 - Pre 65", ""),
                    lowPGS = c(""),
                    highPGS = c(""))
    
    #table[1, 2] <- "BOLD Effect of health shock on smoking probability"
    
    table[2, 1] <- ""
    table[2, 2] <- "Low PGS"
    table[2, 3] <- "High PGS"
    
    get_coef_str <- function(lin_hyp) {
      return(paste0(round(lin_hyp$test$coef, digits = 3), signif_stars(lin_hyp$test$pvalues[1])))
    }
    
    table[3, 2] <- get_coef_str(glht_low_pre65)
    table[4, 2] <- paste0("(", round(glht_low_pre65$test$sigma, digits = 3), ")")
    table[3, 3] <- get_coef_str(glht_high_pre65)
    table[4, 3] <- paste0("(", round(glht_high_pre65$test$sigma, digits = 3), ")")
    
    table[5, 2] <- get_coef_str(glht_low_post65)
    table[6, 2] <- paste0("(", round(glht_low_post65$test$sigma, digits = 3), ")")
    table[5, 3] <- get_coef_str(glht_high_post65)
    table[6, 3] <- paste0("(", round(glht_high_post65$test$sigma, digits = 3), ")")
    
    #table[7, 2] <- "BOLD Effect of health insurance on effect of health shock"
    table[7, 2] <- "Low PGS"
    table[7, 3] <- "High PGS"
    
    table[8, 2] <- get_coef_str(glht_prelow_vs_postlow)
    table[9, 2] <- paste0("(", round(glht_prelow_vs_postlow$test$sigma, digits = 3), ")")
    table[8, 3] <- get_coef_str(glht_prehigh_vs_posthigh)
    table[9, 3] <- paste0("(", round(glht_prehigh_vs_posthigh$test$sigma, digits = 3), ")")
    
    #table[11, 2] <- "BOLD Differential effect of health insurance by genetic group"
    table[10, 2] <- "High PGS "
    table[10, 3] <- "- low PGS"
    
    table[11, 2] <- get_coef_str(glht_slopehigh_vs_slopelow)
    table[12, 2] <- paste0("(", round(glht_slopehigh_vs_slopelow$test$sigma, digits = 3), ")")
    
    names(table)[1] <- ""
    names(table)[2] <- "Low PGS"
    names(table)[3] <- "High PGS"
    
    # Return
    results_list = list("plot" = plot, "overview_table" = table)
    return(results_list)
    
  } else {
    table <- tibble(HI = c( "","", "Pre 65", "", "Post 65", "", "", "Post 65 - Pre 65", "", "", "Post 65 - Pre 65", ""),
                    lowPGS = c(""),
                    midPGS = c(""),
                    highPGS = c(""))
    
    #table[1, 2] <- "BOLD Effect of health shock on smoking probability"
    table[2, 2] <- "Low PGS"
    table[2, 3] <- "Middle PGS"
    table[2, 4] <- "High PGS"
    
    get_coef_str <- function(lin_hyp) {
      return(paste0(round(lin_hyp$test$coef, digits = 3), signif_stars(lin_hyp$test$pvalues[1])))
    }
    
    table[3, 2] <- get_coef_str(glht_low_pre65)
    table[4, 2] <- paste0("(", round(glht_low_pre65$test$sigma, digits = 3), ")")
    table[3, 3] <- get_coef_str(glht_mid_pre65)
    table[4, 3] <- paste0("(", round(glht_mid_pre65$test$sigma, digits = 3), ")")
    table[3, 4] <- get_coef_str(glht_high_pre65)
    table[4, 4] <- paste0("(", round(glht_high_pre65$test$sigma, digits = 3), ")")
    
    table[5, 2] <- get_coef_str(glht_low_post65)
    table[6, 2] <- paste0("(", round(glht_low_post65$test$sigma, digits = 3), ")")
    table[5, 3] <- get_coef_str(glht_mid_post65)
    table[6, 3] <- paste0("(", round(glht_mid_post65$test$sigma, digits = 3), ")")
    table[5, 4] <- get_coef_str(glht_high_post65)
    table[6, 4] <- paste0("(", round(glht_high_post65$test$sigma, digits = 3), ")")
    
    #table[7, 2] <- "BOLD Effect of health insurance on effect of health shock"
    table[7, 2] <- "Low PGS"
    table[7, 3] <- "Middle PGS"
    table[7, 4] <- "High PGS"
    
    table[8, 2] <- get_coef_str(glht_prelow_vs_postlow)
    table[9, 2] <- paste0("(", round(glht_prelow_vs_postlow$test$sigma, digits = 3), ")")
    table[8, 3] <- get_coef_str(glht_premid_vs_postmid)
    table[9, 3] <- paste0("(", round(glht_premid_vs_postmid$test$sigma, digits = 3), ")")
    table[8, 4] <- get_coef_str(glht_prehigh_vs_posthigh)
    table[9, 4] <- paste0("(", round(glht_prehigh_vs_posthigh$test$sigma, digits = 3), ")")
    
    #table[11, 2] <- "BOLD Differential effect of health insurance by genetic group"
    table[10, 2] <- "High PGS - low PGS"
    table[10, 3] <- "Middle PGS - low PGS"
    table[10, 4] <- "High PGS - middle PGS"
    
    
    table[11, 2] <- get_coef_str(glht_slopehigh_vs_slopelow)
    table[12, 2] <- paste0("(", round(glht_slopehigh_vs_slopelow$test$sigma, digits = 3), ")")
    table[11, 3] <- get_coef_str(glht_slopemid_vs_slopelow)
    table[12, 3] <- paste0("(", round(glht_slopemid_vs_slopelow$test$sigma, digits = 3), ")")
    table[11, 4] <- get_coef_str(glht_slopehigh_vs_slopemid)
    table[12, 4] <- paste0("(", round(glht_slopehigh_vs_slopemid$test$sigma, digits = 3), ")")
    
    names(table)[1] <- ""
    names(table)[2] <- ""
    names(table)[3] <- ""
    names(table)[4] <- ""
    
    # Return
    results_list = list("plot" = plot, "overview_table" = table)
    return(results_list)
    
  }
}

# # For analysis where we run linear hypothesis in stata
# stata_shock_effects <- function(df_coef, name, y_achsis = "%-point-change in smoking probability", x_achsis = c("Low PGS", "High PGS")){
#   df_plot <- tibble(time = c(rep("pre65", 2),rep("post65", 2)),
#                     group = c(rep(c("lowPGS", "highPGS"), 2)),
#                     coeff = c(as.numeric(df_coef[1]),
#                               as.numeric(df_coef[4]),
#                               as.numeric(df_coef[7]),
#                               as.numeric(df_coef[10])),
#                     se    = c(as.numeric(df_coef[2]),
#                               as.numeric(df_coef[5]),
#                               as.numeric(df_coef[8]),
#                               as.numeric(df_coef[11])
#                     ))
#   
#   
#   df_plot$time <- as.factor(df_plot$time)
#   df_plot$time <- factor(df_plot$time, levels = c("pre65", "post65"))
#   df_plot$group <- factor(df_plot$group, levels = c("lowPGS", "highPGS"))
#   
#   plot <- ggplot(df_plot %>% mutate(coeff = coeff * 100, se = se * 100), aes(x=group, y=coeff, color=time)) +
#     geom_point(position = position_dodge(width = 0.1), size = 3) +
#     geom_errorbar(aes(ymin = coeff-1.96*se, ymax = coeff+1.96*se),
#                   width=.1,
#                   position = position_dodge(width = 0.1),
#                   size = 1) +
#     labs(x = NULL,
#          y = y_achsis) +
#     scale_x_discrete(labels = x_achsis)  +
#     scale_color_manual(labels = c("Pre-65", "Post-65"),
#                        name = "Shock timing",
#                        values=c("#F36A6A", "#5CACEE")) +
#     geom_hline(aes(yintercept = 0), color = "gray40", linetype = "dashed") +
#     #ylim(-0.5,0.3) +
#     theme_minimal() +
#     theme(legend.position = c(0.88, 0.1),
#           axis.ticks.x = element_blank())
#   
#   ggsave(paste0("../3_output/shock_effects/", name, "plot.png"),
#          plot = plot, width = 7, height = 7)
#   
#   # Create the table
#   table <- tibble(HI = c("","", "Pre 65", "", "Post 65", "", "", "Post 65 - Pre 65", "","","Post 65 - Pre 65", ""),
#                   lowPGS = c(""),
#                   highPGS = c(""))
#   
#   #table[1, 2] <- "BOLD Effect of health shock on smoking probability"
#   
#   table[2, 1] <- ""
#   table[2, 2] <- "Low PGS"
#   table[2, 3] <- "High PGS"
#   
#   #get_coef_str <- function(lin_hyp) {
#   #  return(paste0(round(lin_hyp$test$coef, digits = 3), signif_stars(lin_hyp$test$pvalues[1])))
#   #}
#   
#   table[3, 2] <- paste0(round(as.numeric(df_coef[1]), digits = 3), signif_stars(pt(as.numeric(df_coef[1])/as.numeric(df_coef[2]), df = as.numeric(df_coef[3]))))
#   table[4, 2] <- round(as.numeric(df_coef[2]), digits = 3)
#   table[3, 3] <- paste0(round(as.numeric(df_coef[4]), digits = 3), signif_stars(pt(as.numeric(df_coef[4])/as.numeric(df_coef[5]), df = as.numeric(df_coef[6]))))
#   table[4, 3] <- round(as.numeric(df_coef[5]), digits = 3)
#   
#   table[5, 2] <- paste0(round(as.numeric(df_coef[7]), digits = 3), signif_stars(pt(as.numeric(df_coef[7])/as.numeric(df_coef[8]), df = as.numeric(df_coef[9]))))
#   table[6, 2] <- round(as.numeric(df_coef[8]), digits = 3)
#   table[5, 3] <- paste0(round(as.numeric(df_coef[10]), digits = 3), signif_stars(pt(as.numeric(df_coef[10])/as.numeric(df_coef[11]), df = as.numeric(df_coef[12]))))
#   table[6, 3] <- round(as.numeric(df_coef[11]), digits = 3)
#   
#   #table[7, 2] <- "BOLD Effect of health insurance on effect of health shock"
#   table[7, 2] <- "Low PGS"
#   table[7, 3] <- "High PGS"
#   
#   table[8, 2] <- paste0(round(as.numeric(df_coef[13]), digits = 3), signif_stars(pt(as.numeric(df_coef[13])/as.numeric(df_coef[14]), df = as.numeric(df_coef[15]))))
#   table[9, 2] <- round(as.numeric(df_coef[14]), digits = 3)
#   table[8, 3] <- paste0(round(as.numeric(df_coef[16]), digits = 3), signif_stars(pt(as.numeric(df_coef[16])/as.numeric(df_coef[17]), df = as.numeric(df_coef[18]))))
#   table[9, 3] <- round(as.numeric(df_coef[17]), digits = 3)
#   
#   #table[11, 2] <- "BOLD Differential effect of health insurance by genetic group"
#   table[10, 2] <- "High PGS "
#   table[10, 3] <- "- low PGS"
#   
#   table[11, 2] <- paste0(round(as.numeric(df_coef[19]), digits = 3), signif_stars(pt(as.numeric(df_coef[19])/as.numeric(df_coef[20]), df = as.numeric(df_coef[21]))))
#   table[12, 2] <- round(as.numeric(df_coef[20]), digits = 3)
#   
#   names(table)[1] <- ""
#   names(table)[2] <- "Low PGS"
#   names(table)[3] <- "High PGS"
#   
#   # Return
#   results_list = list("plot" = plot, "overview_table" = table)
#   return(results_list)
# }

# Function: Regressing and plotting CV probabilities ----------------------
cv_prob <- function(df, name){
  df_plot <- df
  # Create age trend for post-65 period
  df <- df %>%
    filter(age != 65, age != 66) %>%
    mutate(after_67_age = ifelse(age >= 67, age - 67, 0),
           post_67 = ifelse(age >= 67, 1, 0))
  
  # Continuous age variable regression (not used for plot)
  get_stderror <- function(cv_lm) {
    cov <- vcovHC(cv_lm, type = "HC1")
    return(sqrt(diag(cov)))
  }
  
  cv_prob <- lm("cv ~ age + post_67 + post_67:after_67_age", data = df)
  stderror <- get_stderror(cv_prob)
  
  cv_prob_high <- lm("cv ~ age + post_67 + post_67:after_67_age", data = subset(df, high_pgs == 1))
  stderror_high <- get_stderror(cv_prob_high)
  
  cv_prob_low <- lm("cv ~ age + post_67 + post_67:after_67_age", data = subset(df, high_pgs == 0))
  stderror_low <- get_stderror(cv_prob_low)
  
  cv_prob_male <- lm("cv ~ age + post_67 + post_67:after_67_age", data = subset(df, female == 0))
  stderror_male <- get_stderror(cv_prob_male)
  
  cv_prob_female <- lm("cv ~ age + post_67 + post_67:after_67_age", data = subset(df, female == 1))
  stderror_female <- get_stderror(cv_prob_female)
  
  df <- df_plot
  
  # Dummy age variables
  df$f_age <- as.factor(df$age)
  
  # Regression for basic plot
  cv_prob_d <- lm("cv ~ f_age -1",
                  data = df)
  reg_results_d <- coeftest(cv_prob_d, vcov. = vcovHC(cv_prob_d))
  
  reg_results_d_tidy <- tidy(reg_results_d) %>%
    mutate(age = 59 + row_number())
  
  # Basic plot
  plot_result_basic <- ggplot(reg_results_d_tidy %>% mutate(estimate = estimate * 100, std.error = std.error * 100),
                              aes(x = as.factor(age), y = estimate)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin=estimate - 1.96 * std.error, ymax=estimate + 1.96 * std.error)) +
    labs(x = NULL,
         y = "% of respondents reporting a health shock"
    ) +
    theme_minimal()
  ggsave(paste0(name, ".png"), plot = plot_result_basic)
  
  # Regression for plot divided by PGS
  cv_prob_d_high <- lm("cv ~ f_age - 1",
                       data = subset(df, high_pgs == 1))
  cv_prob_d_low <- lm("cv ~ f_age - 1",
                      data = subset(df, high_pgs == 0))
  
  reg_results_d_high <- coeftest(cv_prob_d_high, vcov. = vcovHC(cv_prob_d_high))
  reg_results_d_low <- coeftest(cv_prob_d_low, vcov. = vcovHC(cv_prob_d_low))
  
  reg_results_d_high_tidy <- tidy(reg_results_d_high) %>%
    mutate(age = 59 + row_number(),
           high_pgs = 1)
  reg_results_d_low_tidy <- tidy(reg_results_d_low) %>%
    mutate(age = 59 + row_number(),
           high_pgs = 0)
  
  reg_results_d_pgs_tidy <- rbind(reg_results_d_high_tidy, reg_results_d_low_tidy)
  
  # Plot divided by PGS
  plot_result_pgs <- ggplot(reg_results_d_pgs_tidy %>% mutate(estimate = estimate * 100,
                                                              std.error = std.error * 100,
                                                              not_included = as.factor(ifelse((age != 65) & (age != 66), 0, 1))),
                            aes(x = as.factor(age), y = estimate, color = as.factor(high_pgs))) +
    geom_point(position = position_dodge(width = 0.4), size = 2) +
    geom_errorbar(aes(ymin = estimate - 1.96 * std.error,
                      ymax = estimate + 1.96 * std.error,
                      linetype = not_included),
                  width=.5,
                  position = position_dodge(width = 0.4),
                  size = 1) +
    labs(x = "Age at interview",
         y = "% of respondents reporting a health shock") +
    scale_color_manual(labels = c("Low PGS", "High PGS"),
                       name = "Genetic group",
                       values = c("#F36A6A", "#5CACEE")) +
    scale_linetype_manual(name = "Included",
                          labels = c("Yes", "No"),
                          values = c("solid", "dotted"),
                          guide = 'none') +
    ylim(0, 7) +
    theme_minimal() +
    theme(panel.grid.minor.x = element_blank(),
          legend.position = c(0.12, 0.88))
  ggsave(paste0("../3_output/cv_prob/", name, "plot_pgs.png"),
         plot = plot_result_pgs, width = 8, height = 6)
  
  # Regression for plot divided by gender
  cv_prob_d_male <- lm("cv ~ f_age - 1",
                       data = subset(df, female == 0))
  cv_prob_d_female <- lm("cv ~ f_age - 1",
                         data = subset(df, female == 1))
  
  reg_results_d_male <- coeftest(cv_prob_d_male, vcov. = vcovHC(cv_prob_d_male))
  reg_results_d_female <- coeftest(cv_prob_d_female, vcov. = vcovHC(cv_prob_d_female))
  
  reg_results_d_male_tidy <- tidy(reg_results_d_male) %>%
    mutate(age = 59 + row_number(),
           female = 0)
  reg_results_d_female_tidy <- tidy(reg_results_d_female) %>%
    mutate(age = 59 + row_number(),
           female = 1)
  
  reg_results_d_gender_tidy <- rbind(reg_results_d_male_tidy, reg_results_d_female_tidy)
  
  # Plot divided by gender
  plot_result_gender <- ggplot(reg_results_d_gender_tidy %>% mutate(estimate = estimate * 100,
                                                                    std.error = std.error * 100,
                                                                    not_included = as.factor(ifelse((age != 65) & (age != 66), 0, 1))),
                               aes(x = as.factor(age), y = estimate, color = as.factor(female))) +
    geom_point(position = position_dodge(width = 0.4), size = 2) +
    geom_errorbar(aes(ymin = estimate - 1.96 * std.error,
                      ymax = estimate + 1.96 * std.error,
                      linetype = not_included),
                  width=.5,
                  position = position_dodge(width = 0.4),
                  size = 1) +
    scale_linetype_manual(name = "Included",
                          labels = c("Yes", "No"),
                          values = c("solid", "dotted"),
                          guide = 'none') +
    labs(x = "Age at interview",
         y = "% of respondents reporting a health shock") +
    scale_color_manual(labels = c("Male", "Female"),
                       name = "Gender",
                       values=c("#F36A6A", "#5CACEE")) +
    ylim(0, 7) +
    theme_minimal() +
    theme(panel.grid.minor.x = element_blank(),
          legend.position = c(0.08, 0.88))
  ggsave(paste0("../3_output/cv_prob/", name, "plot_gender.png"),
         plot = plot_result_gender, width = 8, height = 6)
  
  results <- stargazer(cv_prob, cv_prob_low, cv_prob_high, cv_prob_male, cv_prob_female,
                       keep.stat = "n",
                       dep.var.labels = "Probability of shock",
                       column.labels = c("All", "Low PGS", "High PGS", "Male", "Female"),
                       covariate.labels = c("Age", "Post-67 dummy", "Post-67 slope"), float = FALSE)
  results_list = list("plot_basic" = plot_result_basic,"plot_pgs" = plot_result_pgs, "plot_gender" = plot_result_gender, "results" = results)
  return(results_list)
}

# Function: Plot smoking over time ---------------------------------------
over_time <- function(df, outcome = "smoken", name, filedir){
  
  df$outcome <- df[[outcome]]
  
  # Dummy age variables
  df$agem_e12 <- df$agem_e/12
  df$f_age <- as.factor(df$agem_e12)
  
  # Regression for PGS-divided plot
  reg_high <- lm(outcome ~ f_age - 1 + female,
                 data = subset(df, high_pgs == 1))
  reg_low <- lm(outcome ~ f_age - 1 + female,
                data = subset(df, high_pgs == 0))
  
  reg_results_high <- coeftest(reg_high, vcov. = vcovHC(reg_high))
  reg_results_low <- coeftest(reg_low, vcov. = vcovHC(reg_low))
  
  reg_results_high_tidy <- tidy(reg_results_high) %>%
    mutate(age = 60 + (row_number()-1)/12,
           high_pgs = 1)
  
  reg_results_low_tidy <- tidy(reg_results_low) %>%
    mutate(age = 60 + (row_number()-1)/12,
           high_pgs = 0)
  
  reg_results_tidy <- rbind(reg_results_low_tidy, reg_results_high_tidy)
  
  # y-axis label
  if (outcome == "smoken") {
    lab_y_axis <- "% of respondents who smoke"
  }
  else if (outcome == "retired") {
    lab_y_axis <- "% of respondents who retired"
  }
  else if (outcome == "govmr") {
    lab_y_axis <- "% of respondents enrolled in Medicare"
  }
  else if (outcome == "iearn") {
    lab_y_axis <- "Individual earnings"
  }
  else if (outcome == "cv") {
    lab_y_axis <- "% of respondents rerporting a health shock"
  }
  else {
    lab_y_axis <- paste0(outcome)
  }
  
  # Regression for PGS-divided plot
  plot_by_pgs <- ggplot(reg_results_tidy %>% mutate(estimate = estimate * 100),
                        aes(x = age, y = estimate, color = as.factor(high_pgs))) +
    geom_smooth(method = "auto",
                #formula = y ~ x + I(x^2) + I(x^3),
                size = 1) +
    geom_point() +
    # geom_errorbar(aes(ymin = estimate - 1.96 * std.error, ymax = estimate + 1.96 * std.error),
    #               width=.5,
    #               position = position_dodge(width = 0.4),
    #               size = 1) +
    labs(x = "Age at interview",
         y = lab_y_axis) +
    scale_color_manual(labels = c("Low PGS", "High PGS"),
                       name = "Genetic group",
                       values=c("#F36A6A", "#5CACEE")) +
    theme_minimal() +
    scale_x_continuous(breaks = seq(60, 70, 2), labels = seq(60, 70, 2)) +
    theme(panel.grid.minor.x = element_blank(),
          legend.position = c(0.85, 0.88))
  ggsave(paste0(filedir),
         plot = plot_by_pgs, width = 9, height = 6)
  
  return(plot_by_pgs)
}

# Function: Plot smoking over time with discontinuity at 65----------------------
over_time_discontinuity <- function(df, outcome = "smoken", name, cutoff = "65", filedir){
  
  df$outcome <- df[[outcome]]
  
  # Dummy age variables
  df$agem_e12 <- df$agem_e/12
  df$f_age <- as.factor(df$agem_e12)
  
  # Regression for PGS-divided plot
  reg_high <- lm(outcome ~ f_age - 1 + female,
                 data = subset(df, high_pgs == 1))
  reg_low <- lm(outcome ~ f_age - 1 + female,
                data = subset(df, high_pgs == 0))
  
  reg_results_high <- coeftest(reg_high, vcov. = vcovHC(reg_high))
  reg_results_low <- coeftest(reg_low, vcov. = vcovHC(reg_low))
  
  reg_results_high_tidy <- tidy(reg_results_high) %>%
    mutate(age = 60 + (row_number()-1)/12,
           high_pgs = 1)
  
  reg_results_low_tidy <- tidy(reg_results_low) %>%
    mutate(age = 60 + (row_number()-1)/12,
           high_pgs = 0)
  
  reg_results_tidy <- rbind(reg_results_low_tidy, reg_results_high_tidy)
  
  # y-axis label
  if (outcome == "smoken") {
    lab_y_axis <- "% of respondents who smoke"
  }
  else if (outcome == "retired") {
    lab_y_axis <- "% of respondents who retired"
  }
  else if (outcome == "govmr") {
    lab_y_axis <- "% of respondents enrolled in Medicare"
  }
  else if (outcome == "iearn") {
    lab_y_axis <- "Individual earnings"
  }
  else if (outcome == "cv") {
    lab_y_axis <- "% of respondents rerporting a health shock"
  }
  else {
    lab_y_axis <- paste0(outcome)
  }
  
  # add some final columns to improve plots and allow for age discontinuity
  data_for_discont_plot <- reg_results_tidy %>%
    mutate(estimate = estimate * 100,
           above_64 = age >= cutoff)
  
  # Regression for PGS-divided plot
  plot_by_pgs <- ggplot(mapping=aes(x = age, y = estimate, color = as.factor(high_pgs))) +
    geom_smooth(data = data_for_discont_plot %>% filter(!above_64),
                method = "auto",
                #formula = y ~ x + I(x^2) + I(x^3),
                size = 1) +
    geom_point(data = data_for_discont_plot %>% filter(!above_64)) +
    geom_smooth(data = data_for_discont_plot %>% filter(above_64),
                method = "auto",
                #formula = y ~ x + I(x^2) + I(x^3),
                size = 1) +
    geom_point(data = data_for_discont_plot %>% filter(above_64)) +
    # geom_errorbar(aes(ymin = estimate - 1.96 * std.error, ymax = estimate + 1.96 * std.error),
    #               width=.5,
    #               position = position_dodge(width = 0.4),
    #               size = 1) +
    labs(x = "Age at interview",
         y = lab_y_axis) +
    scale_color_manual(labels = c("Low PGS", "High PGS"),
                       name = "Genetic group",
                       values=c("#F36A6A", "#5CACEE")) +
    theme_minimal() +
    scale_x_continuous(breaks = seq(60, 70, 2), labels = seq(60, 70, 2)) +
    theme(panel.grid.minor.x = element_blank(),
          legend.position = c(0.85, 0.78))
  ggsave(paste0(filedir),
         plot = plot_by_pgs, width = 9, height = 6)
  
  return(plot_by_pgs)
}

# Function: RDD style regressions----------------------
rdd_over_time <- function(df, outcome, cutoff, treatment_name) {
  require("stargazer","plm")
  
  # create age polynomials to avoid multicolinearity
  df <- df %>% mutate(agem = agem_e - mean(agem_e, na.rm = TRUE),
                      agem2 = agem^2,
                      agem3 = agem^3,
                      agem4 = agem^4)
  
  #RDD reg############
  df$outcome <- df[[outcome]]
  
  #new dummy for treatment
  df$cutoff <- ifelse(df$age >= cutoff, 1, 0)
  
  # Fully interacted model
  rdd_reg <- felm(outcome ~ cutoff
                  + cutoff:high_pgs
                  + agem
                  + agem2
                  + agem3
                  + agem4
                  + cutoff:agem
                  + cutoff:agem2
                  + cutoff:agem3
                  + cutoff:agem4
                  + high_pgs:agem
                  + high_pgs:agem2
                  + high_pgs:agem3
                  + high_pgs:agem4
                  + cutoff:high_pgs:agem
                  + cutoff:high_pgs:agem2
                  + cutoff:high_pgs:agem3
                  + cutoff:high_pgs:agem4 
                  | hhidpn + year | #fe
                    0 | #instruments
                    hhidpn, # clustering
                  df,
                  na.action=na.omit)
  
  # stargazer(rdd_reg_high, rdd_reg_low, type = "text", title = paste("RDD for", treatment_name, "years, outcome =", outcome), dep.var.labels   = outcome,
  #           covariate.labels = c(treatment_name, "age", "agepoly2", "agepoly3", "agepoly4", paste(treatment_name,"x age"),
  #                                paste(treatment_name, "x agepoly2"), paste(treatment_name, "x agepoly3"), paste(treatment_name, "x agepoly4")),
  #                                column.labels = c("High PGS", "Low PGS"), out = paste("../3_output/over_time/", outcome,"rdd", cutoff, ".tex"))
  # 
  
  summary(rdd_reg)
  
  #Robustness Check
  
  cutoffb <- as.numeric(cutoff)*12
  
  # HIGH PGS: checking the results using RDrobust
  dfhigh=subset(df, high_pgs==1)
  highrdrobust <- (rdrobust(dfhigh$outcome, dfhigh$agem_e, c = as.numeric(cutoffb), cluster = dfhigh$hhidpn, all = TRUE)) #agem_e is used for input simplicity compared to agem
  
  # LOW PGS: checking the results using RDrobust
  dflow=subset(df, high_pgs==0)
  lowrdrobust <- (rdrobust(dflow$outcome, dflow$agem_e, c = as.numeric(cutoffb), cluster = dflow$hhidpn, all = TRUE)) #agem_e is used for input simplicity compared to agem
  
  outList <- list("rdd_reg" = rdd_reg, "highrdrobust" = highrdrobust, "lowrdrobust" = lowrdrobust)
  return(outList)
}

# Function: Histogram of PGS for smoking by smoking status ----------------
histogram_full <- function(df, name) {
  hist_full <- ggplot(df %>%
                        filter(!is.na(smoken)) %>%
                        group_by(hhidpn) %>%
                        arrange(year) %>%
                        filter(row_number() == 1) %>%
                        ungroup(),
                      aes(x = si_1, y = ..density.., group = as.factor(smoken), fill = as.factor(smoken))) +
    geom_histogram(position = "identity",
                   alpha = 0.5,
                   binwidth = 0.1) +
    labs(x = "PGS for smoking",
         y = "Density") +
    scale_fill_manual(labels = c("Non-smoker", "Smoker"),
                      name = "Smoking status",
                      values=c("#F36A6A", "#5CACEE")) +
    theme_minimal() +
    theme(panel.grid.minor.x = element_blank(),
          legend.position = c(0.88, 0.88))
  ggsave(paste0("../3_output/make_histograms/", name, "smoken.png"),
         plot = hist_full, width = 9, height = 6)
  return(hist_full)
}

# Function: Unlist columns ------------------------------------------------
unlist_columns <- function(df, columns) {
  for (col in columns) {
    df[[col]] <- unlist(df[[col]])
  }
  return(df)
}

# # Function: Convenient wrapper around reg_hipgs and shock_effects ---------
# reg_shock <- function(df, min_age, max_age, ptu_perc, add_info = '') {
#   # Part of file name
#   f_identifier <- paste0("robustness_", min_age, max_age, "_", ptu_perc, "_cv", add_info)
#   reg_temp <- reg_hipgs(df, outcome="smoken", shock_var= "cv", name = f_identifier)
#   sh <- shock_effects(reg_temp, name = f_identifier)
#   return(list(reg = reg_temp, sh = sh))
# }
#


# Cessation rates graph ---------------------------------------------------
cess_rates_graph <- function(df, fileName) {
  
  # Create baseline datasets
  # "lowPGS"
  baseline_temp <- df %>%
    group_by(hhidpn) %>%
    arrange(year) %>%
    mutate(no_waves = n(),
           max_cv = max(cv),
           age_at_shock = ifelse(cv == 1, age, 0),
           age_at_shock = max(age_at_shock)) %>%
    filter(max_cv == 1, row_number() == 1) %>%
    ungroup() %>%
    select(hhidpn, age, age_at_shock, si_1,
           no_waves, female, smoken, uninsured, max_cv,
           high_pgs, unins, raedyrs, iearn)
  
  baseline_low <- baseline_temp %>%
    filter(high_pgs == 0)
  
  # "highPGS"
  baseline_high <- baseline_temp %>%
    filter(high_pgs == 1)
  
  # Create pre65 and post65 sub-datasets for "lowPGS"
  # Remember that shocks at 65 and 66 should be excluded in
  # the past in dataset
  baseline_low_pre65 <- baseline_low %>% filter(age_at_shock < 65 & age_at_shock > 0)
  baseline_low_post65 <- baseline_low %>% filter(age_at_shock > 65)
  
  # Create pre65 and post65 sub-datasets for "highPGS"
  baseline_high_pre65 <- baseline_high %>% filter(age_at_shock < 65 & age_at_shock > 0)
  baseline_high_post65 <- baseline_high %>% filter(age_at_shock > 65)
  
  
  # Calculate mean and sd of cessation rates
  low_pre_cess <- mean_cess_baseline_smk(df, baseline_low_pre65, return_se = TRUE)
  low_post_cess <- mean_cess_baseline_smk(df, baseline_low_post65, return_se = TRUE)
  high_pre_cess <- mean_cess_baseline_smk(df, baseline_high_pre65, return_se = TRUE)
  high_post_cess <- mean_cess_baseline_smk(df, baseline_high_post65, return_se = TRUE)
  
  # Make graph of results
  df_plot <- tibble(time = c(rep("pre65", 2),rep("post65", 2)),
                    group = c(rep(c("lowPGS", "highPGS"), 2)),
                    coeff = c(low_pre_cess$rate,
                              high_pre_cess$rate,
                              low_post_cess$rate,
                              high_post_cess$rate),
                    se = c(low_pre_cess$sd*100/sqrt(low_pre_cess$no),
                           high_pre_cess$sd*100/sqrt(high_pre_cess$no),
                           low_post_cess$sd*100/sqrt(low_post_cess$no),
                           high_post_cess$sd*100/sqrt(high_post_cess$no)
                    ))
  
  df_plot$time <- as.factor(df_plot$time)
  df_plot$time <- factor(df_plot$time, levels = c("pre65", "post65"))
  df_plot$coeff <- as.numeric(df_plot$coeff)
  df_plot$se <- as.numeric(df_plot$se)
  
  cess_rates_plot1 <- ggplot(df_plot, aes(x=time, y=coeff, color=group)) +
    geom_point(position = position_dodge(width = 0.1), size = 3) +
    #geom_line(aes(group = group), size = 1) +
    geom_errorbar(aes(ymin = coeff-1.96*se, ymax = coeff+1.96*se),
                  width=.1,
                  position = position_dodge(width = 0.1),
                  size = 1) +
    labs(x = NULL,
         y = "Avg. cessation rate for baseline smokers") +
    scale_x_discrete(labels = c("Pre-65\n(uninsured)","Post-65\n(eligible for Medicare)"))  +
    scale_color_manual(labels = c("High PGS", "Low PGS"),
                       name = "Genetic group",
                       values=c("#F36A6A", "#5CACEE")) +
    geom_hline(aes(yintercept = 0), color = "gray40", linetype = "dashed") +
    #ylim(-0.5,0.3) +
    theme_minimal() +
    theme(legend.position = c(0.88, 0.1),
          axis.ticks.x = element_blank())
  
  ggsave(paste0(fileName, "plot1.png"),
         plot = cess_rates_plot1, width = 7, height = 7)
  
  
  
  
  df_plot$time <- as.factor(df_plot$time)
  df_plot$group <- factor(df_plot$group, levels = c("lowPGS", "highPGS"))
  
  
  cess_rates_plot2 <- ggplot(df_plot, aes(x=group, y=coeff, color=time)) +
    geom_point(position = position_dodge(width = 0.1), size = 3) +
    geom_errorbar(aes(ymin = coeff-1.96*se, ymax = coeff+1.96*se),
                  width=.1,
                  position = position_dodge(width = 0.1),
                  size = 1) +
    labs(x = NULL,
         y = "Avg. cessation rate for baseline smokers") +
    scale_x_discrete(labels = c("Low PGS","High PGS"))  +
    scale_color_manual(labels = c("Pre-65", "Post-65"),
                       name = "Shock timing",
                       values=c("#F36A6A", "#5CACEE")) +
    geom_hline(aes(yintercept = 0), color = "gray40", linetype = "dashed") +
    theme_minimal() +
    theme(legend.position = c(0.88, 0.1),
          axis.ticks.x = element_blank())
  
  
  ggsave(paste0(fileName, "plot2.png"),
         plot = cess_rates_plot2, width = 7, height = 7)
  
  results_list <- list("plot1" = cess_rates_plot1, "plot2" = cess_rates_plot2)
  return(results_list)
  
}


# Function: Regression with continuous PGS variable --------------------------------
reg_continuous <- function(df, shock_var, name) {
  df_reg <- df
  
  df_reg$shock <- df_reg[[shock_var]]
  reg_shock <- plm(smoken ~ shock +
                     post65 +
                     shock:post65 +
                     shock:unins +
                     post65:unins +
                     si_1:shock +
                     si_1:post65 +
                     shock:post65:unins +
                     si_1:shock:unins +
                     si_1:shock:post65 +
                     si_1:post65:unins +
                     shock:post65:unins:si_1 +
                     agem + agem2 + agem3,
                   df_reg,
                   na.action=na.omit,
                   index=c("hhidpn", "year"),
                   model="within",
                   effect="twoway")
  
  # Specify robust or normal standard errors
  # reg_shock_se <- coeftest(reg_shock, vcov. = vcovHC(reg_shock, cluster = "group"))
  reg_shock_se <- reg_shock
  
  cov <- vcovHC(reg_shock, cluster = "group")
  stderror <- sqrt(diag(cov))
  
  return(reg_shock_se)
}

###################################
# New Functions: AlexanderT. 06.06.2020

reg_hipgs_felm <- function(df, outcome, shock_var, PGS3 = FALSE, name) {
  df_reg <- df
  
  
  df_reg$shock <- df_reg[[shock_var]]
  df_reg$outcome <- df_reg[[outcome]]
  
  if (PGS3 == TRUE) {
    reg_shock <- felm(outcome ~ shock +
                        post65 +
                        shock:post65 +
                        shock:unins +
                        post65:unins +
                        mid_pgs:shock +
                        mid_pgs:post65 +
                        high_pgs:shock +
                        high_pgs:post65 +
                        shock:post65:unins +
                        mid_pgs:shock:unins +
                        mid_pgs:shock:post65 +
                        mid_pgs:post65:unins +
                        high_pgs:shock:unins +
                        high_pgs:shock:post65 +
                        high_pgs:post65:unins +
                        shock:post65:unins:mid_pgs +
                        shock:post65:unins:high_pgs +
                        agem + agem2 + agem3 | hhidpn + year | 0 | hhidpn,
                      data = df_reg, 
                      na.action=na.omit)
    
    
    # Specify robust or normal standard errors
    # reg_shock_se <- coeftest(reg_shock, vcov. = vcov(reg_shock, cluster = "group"))
    # reg_shock_se <- reg_shock
    return(reg_shock)
  } else {
    reg_shock <- felm(outcome ~ shock +
                        post65 +
                        shock:post65 +
                        shock:unins +
                        post65:unins +
                        high_pgs:shock +
                        high_pgs:post65 +
                        shock:post65:unins +
                        high_pgs:shock:unins +
                        high_pgs:shock:post65 +
                        high_pgs:post65:unins +
                        shock:post65:unins:high_pgs +
                        agem + agem2 + agem3 | hhidpn + year | 0 | hhidpn,
                      data = df_reg,
                      na.action=na.omit)
    
    # Specify robust or normal standard errors
    # reg_shock_se <- coeftest(reg_shock, vcov. = vcov(reg_shock, cluster = "group"))
    # reg_shock_se <- reg_shock
    
    # cov <- vcov(reg_shock, cluster = "group")
    # stderror <- sqrt(diag(cov))
    
    return(reg_shock)
  }
}



# new shock effects function

shock_effects_felm <- function(regression, name, percent = TRUE, y_achsis = "%-point-change in smoking probability", x_achsis = c("Low PGS", "High PGS"), PGS3 = FALSE){
  # Effect of a health shock for low-PGS people BEFORE 65
  glht_low_pre65 <- summary(glht(regression,
                                 linfct = c("shock + shock:unins = 0"),
                                 vcov = vcov(regression, cluster = "group")))
  if (PGS3 == TRUE) {
    # Effect of a health shock for mid-PGS people BEFORE 65
    glht_mid_pre65 <-summary(glht(regression,
                                  linfct = c("shock + shock:unins +
                                            shock:mid_pgs +
                                            shock:unins:mid_pgs = 0"),
                                  vcov = vcov(regression, cluster = "group")))
  }
  # Effect of a health shock for high-PGS people BEFORE 65
  glht_high_pre65 <-summary(glht(regression,
                                 linfct = c("shock + shock:unins +
                                            shock:high_pgs +
                                            shock:unins:high_pgs = 0"),
                                 vcov = vcov(regression, cluster = "group")))
  
  # Effect of a health shock for low-PGS people AFTER 65
  glht_low_post65 <-summary(glht(regression,
                                 linfct = c("shock + shock:post65 + shock:unins +
                                             shock:post65:unins = 0"),
                                 vcov = vcov(regression, cluster = "group")))
  if (PGS3 == TRUE) {
    # Effect of a health shock for mid-PGS people AFTER 65
    glht_mid_post65 <- summary(glht(regression,
                                    linfct = c("shock +
                                              shock:post65 +
                                              shock:unins +
                                              shock:mid_pgs +
                                              shock:post65:unins +
                                              shock:unins:mid_pgs +
                                              shock:post65:mid_pgs +
                                              shock:post65:unins:mid_pgs = 0"),
                                    vcov = vcov(regression, cluster = "group")))
  }
  
  # Effect of a health shock for high-PGS people AFTER 65
  glht_high_post65 <- summary(glht(regression,
                                   linfct = c("shock +
								               shock:post65 +
      											   shock:unins +
		      									   shock:high_pgs +
								      			   shock:post65:unins +
											         shock:unins:high_pgs +
			      								   shock:post65:high_pgs +
									      		   shock:post65:unins:high_pgs = 0"),
                                   vcov = vcov(regression, cluster = "group")))
  
  
  # Difference in the effect of a health shock for low-PGS people after 65 minus before 65
  ####Check: it should be the same as glht_low_post65 - glht_low_pre65
  glht_prelow_vs_postlow <- summary(glht(regression,
                                         linfct = c("shock:post65 +
                                                    shock:post65:unins = 0"),
                                         vcov = vcov(regression, cluster = "group")))
  
  # Difference in the effect of a health shock for mid-PGS people before and after 65
  ####Check: it should be the same as glht_mid_post65 - glht_mid_pre65
  if (PGS3 == TRUE) {
    glht_premid_vs_postmid <- summary(glht(regression,
                                           linfct = c("shock:post65 +
                                                      shock:post65:unins +
                                                      shock:post65:mid_pgs +
                                                      shock:post65:unins:mid_pgs = 0"),
                                           vcov = vcov(regression, cluster = "group")))
  }
  
  # Difference in the effect of a health shock for high-PGS people before and after 65
  ####Check: it should be the same as glht_high_post65 - glht_high_pre65
  glht_prehigh_vs_posthigh <- summary(glht(regression,
                                           linfct = c("shock:post65 +
                                                      shock:post65:unins +
                                                      shock:post65:high_pgs +
                                                      shock:post65:unins:high_pgs = 0"),
                                           vcov = vcov(regression, cluster = "group")))
  
  # Difference in the effect of a health shock between low- and high-PGS people BEFORE 65
  ####Check: it should be the same as glht_high_pre65 - glht_low_pre65
  glht_prelow_vs_prehigh <- summary(glht(regression,
                                         linfct = c("shock:high_pgs +
                                                    shock:unins:high_pgs = 0"),
                                         vcov = vcov(regression, cluster = "group")))
  if (PGS3 == TRUE) {
    # Difference in the effect of a health shock between low- and mid-PGS people BEFORE 65
    ####Check: it should be the same as glht_mid_pre65 - glht_low_pre65
    glht_prelow_vs_premid <- summary(glht(regression,
                                          linfct = c("shock:mid_pgs +
                                                    shock:unins:mid_pgs = 0"),
                                          vcov = vcov(regression, cluster = "group")))
    
    # Difference in the effect of a health shock between mid- and high-PGS people BEFORE 65
    ####Check: it should be the same as glht_high_pre65 - glht_mid_pre65
    glht_premid_vs_prehigh <- summary(glht(regression,
                                           linfct = c("shock:high_pgs +
                                                      shock:unins:high_pgs -
                                                      shock:mid_pgs -
                                                      shock:unins:mid_pgs = 0"),
                                           vcov = vcov(regression, cluster = "group")))
  }
  
  # Difference in the effect of a health shock between low- and high-PGS people AFTER 65
  ####Check: it should be the same as glht_high_post65 - glht_low_post65
  glht_postlow_vs_posthigh <- summary(glht(regression,
                                           linfct = c("shock:high_pgs +
                                                      post65:high_pgs +
                                                      shock:unins:high_pgs +
                                                      shock:post65:high_pgs +
                                                      post65:unins:high_pgs +
                                                      shock:post65:unins:high_pgs = 0"),
                                           vcov = vcov(regression, cluster = "group")))
  if (PGS3 == TRUE) {
    # Difference in the effect of a health shock between low- and mid-PGS people AFTER 65
    ####Check: it should be the same as glht_mid_post65 - glht_low_post65
    glht_postlow_vs_postmid <- summary(glht(regression,
                                            linfct = c("shock:mid_pgs +
                                                        post65:mid_pgs +
                                                        shock:unins:mid_pgs +
                                                        shock:post65:mid_pgs +
                                                        post65:unins:mid_pgs +
                                                        shock:post65:unins:mid_pgs = 0"),
                                            vcov = vcov(regression, cluster = "group")))
    
    # Difference in the effect of a health shock between low- and high-PGS people AFTER 65
    ####Check: it should be the same as glht_high_post65 - glht_low_post65
    glht_postmid_vs_posthigh <- summary(glht(regression,
                                             linfct = c("shock:high_pgs +
                                                        post65:high_pgs +
                                                        shock:unins:high_pgs +
                                                        shock:post65:high_pgs +
                                                        post65:unins:high_pgs +
                                                        shock:post65:unins:high_pgs -
                                                        shock:mid_pgs -
                                                        post65:mid_pgs -
                                                        shock:unins:mid_pgs -
                                                        shock:post65:mid_pgs -
                                                        post65:unins:mid_pgs -
                                                        shock:post65:unins:mid_pgs = 0"),
                                             vcov = vcov(regression, cluster = "group")))
    
  }
  
  # Difference in the effect of health insurance status on the effect of a health shock
  # between high- and low-PGS people
  ####Check: this should be the triple difference
  glht_slopehigh_vs_slopelow <- summary(glht(regression,
                                             linfct = c("shock:post65:high_pgs +
                                                        shock:post65:unins:high_pgs = 0"),
                                             vcov = vcov(regression, cluster = "group")))
  
  if (PGS3 == TRUE) {
    # Difference in the effect of health insurance status on the effect of a health shock
    # between mid- and low-PGS people
    ####Check: this should be the triple difference
    glht_slopemid_vs_slopelow <- summary(glht(regression,
                                              linfct = c("shock:post65:mid_pgs +
                                                          shock:post65:unins:mid_pgs = 0"),
                                              vcov = vcov(regression, cluster = "group")))
    
    # Difference in the effect of health insurance status on the effect of a health shock
    # between high- and mid-PGS people
    ####Check: this should be the triple difference
    glht_slopehigh_vs_slopemid <- summary(glht(regression,
                                               linfct = c("shock:post65:high_pgs +
                                                          shock:post65:unins:high_pgs -
                                                          shock:post65:mid_pgs -
                                                          shock:post65:unins:mid_pgs = 0"),
                                               vcov = vcov(regression, cluster = "group")))
  }
  
  # Make graph of results
  
  if (PGS3 == FALSE) {
    df_plot <- tibble(time = c(rep("pre65", 2),rep("post65", 2)),
                      group = c(rep(c("lowPGS", "highPGS"), 2)),
                      coeff = c(glht_low_pre65$test$coef,
                                glht_high_pre65$test$coef,
                                glht_low_post65$test$coef,
                                glht_high_post65$test$coef),
                      se    = c(glht_low_pre65$test$sigma,
                                glht_high_pre65$test$sigma,
                                glht_low_post65$test$sigma,
                                glht_high_post65$test$sigma
                      ))
    
    df_plot$time <- as.factor(df_plot$time)
    df_plot$time <- factor(df_plot$time, levels = c("pre65", "post65"))
    df_plot$group <- factor(df_plot$group, levels = c("lowPGS", "highPGS"))
    
    # if beta cannot interpreted as percentage, we have to divide by 100 as it will be mulitplied later
    if (percent==FALSE){
      df_plot <- df_plot %>% mutate(coeff = coeff / 100, se = se / 100)
    }
    
    
    plot <- ggplot(df_plot %>% mutate(coeff = coeff * 100, se = se * 100), aes(x=group, y=coeff, color=time)) +
      geom_point(position = position_dodge(width = 0.1), size = 3) +
      geom_errorbar(aes(ymin = coeff-1.96*se, ymax = coeff+1.96*se),
                    width=.1,
                    position = position_dodge(width = 0.1),
                    size = 1) +
      labs(x = NULL,
           y = y_achsis) +
      scale_x_discrete(labels = x_achsis)  +
      scale_color_manual(labels = c("Pre-65", "Post-65"),
                         name = "Shock timing",
                         values=c("#F36A6A", "#5CACEE")) +
      geom_hline(aes(yintercept = 0), color = "gray40", linetype = "dashed") +
      #ylim(-0.5,0.3) +
      theme_minimal() +
      theme(legend.position = c(0.88, 0.1),
            axis.ticks.x = element_blank())
    
    ggsave(paste0("../3_output/shock_effects/", name, "plot.png"),
           plot = plot, width = 7, height = 7)
  } else {
    df_plot <- tibble(time = c(rep("pre65", 3),rep("post65", 3)),
                      group = c(rep(c("lowPGS", "middlePGS", "highPGS"), 2)),
                      coeff = c(glht_low_pre65$test$coef,
                                glht_mid_pre65$test$coef,
                                glht_high_pre65$test$coef,
                                glht_low_post65$test$coef,
                                glht_mid_post65$test$coef,
                                glht_high_post65$test$coef),
                      se    = c(glht_low_pre65$test$sigma,
                                glht_mid_pre65$test$sigma,
                                glht_high_pre65$test$sigma,
                                glht_low_post65$test$sigma,
                                glht_mid_post65$test$sigma,
                                glht_high_post65$test$sigma
                      ))
    
    df_plot$time <- as.factor(df_plot$time)
    df_plot$time <- factor(df_plot$time, levels = c("pre65", "post65"))
    df_plot$group <- factor(df_plot$group, levels = c("lowPGS","middlePGS", "highPGS"))
    
    plot <- ggplot(df_plot %>% mutate(coeff = coeff * 100, se = se * 100), aes(x=group, y=coeff, color=time)) +
      geom_point(position = position_dodge(width = 0.1), size = 3) +
      geom_errorbar(aes(ymin = coeff-1.96*se, ymax = coeff+1.96*se),
                    width=.1,
                    position = position_dodge(width = 0.1),
                    size = 1) +
      labs(x = NULL,
           y = y_achsis) +
      scale_x_discrete(labels = x_achsis)  +
      scale_color_manual(labels = c("Pre-65", "Post-65"),
                         name = "Shock timing",
                         values=c("#F36A6A", "#5CACEE")) +
      geom_hline(aes(yintercept = 0), color = "gray40", linetype = "dashed") +
      #ylim(-0.5,0.3) +
      theme_minimal() +
      theme(legend.position = c(0.88, 0.1),
            axis.ticks.x = element_blank())
    
    ggsave(paste0("../3_output/shock_effects/", name, "plot.png"),
           plot = plot, width = 7, height = 7)
  }
  
  # Make an overview table for the main paper
  if (PGS3 == FALSE) {
    table <- tibble(HI = c("","", "Pre 65", "", "Post 65", "", "", "Post 65 - Pre 65", "","","Post 65 - Pre 65", ""),
                    lowPGS = c(""),
                    highPGS = c(""))
    
    #table[1, 2] <- "BOLD Effect of health shock on smoking probability"
    
    table[2, 1] <- ""
    table[2, 2] <- "Low PGS"
    table[2, 3] <- "High PGS"
    
    get_coef_str <- function(lin_hyp) {
      return(paste0(round(lin_hyp$test$coef, digits = 3), signif_stars(lin_hyp$test$pvalues[1])))
    }
    
    table[3, 2] <- get_coef_str(glht_low_pre65)
    table[4, 2] <- paste0("(", round(glht_low_pre65$test$sigma, digits = 3), ")")
    table[3, 3] <- get_coef_str(glht_high_pre65)
    table[4, 3] <- paste0("(", round(glht_high_pre65$test$sigma, digits = 3), ")")
    
    table[5, 2] <- get_coef_str(glht_low_post65)
    table[6, 2] <- paste0("(", round(glht_low_post65$test$sigma, digits = 3), ")")
    table[5, 3] <- get_coef_str(glht_high_post65)
    table[6, 3] <- paste0("(", round(glht_high_post65$test$sigma, digits = 3), ")")
    
    #table[7, 2] <- "BOLD Effect of health insurance on effect of health shock"
    table[7, 2] <- "Low PGS"
    table[7, 3] <- "High PGS"
    
    table[8, 2] <- get_coef_str(glht_prelow_vs_postlow)
    table[9, 2] <- paste0("(", round(glht_prelow_vs_postlow$test$sigma, digits = 3), ")")
    table[8, 3] <- get_coef_str(glht_prehigh_vs_posthigh)
    table[9, 3] <- paste0("(", round(glht_prehigh_vs_posthigh$test$sigma, digits = 3), ")")
    
    #table[11, 2] <- "BOLD Differential effect of health insurance by genetic group"
    table[10, 2] <- "High PGS "
    table[10, 3] <- "- low PGS"
    
    table[11, 2] <- get_coef_str(glht_slopehigh_vs_slopelow)
    table[12, 2] <- paste0("(", round(glht_slopehigh_vs_slopelow$test$sigma, digits = 3), ")")
    
    names(table)[1] <- ""
    names(table)[2] <- "Low PGS"
    names(table)[3] <- "High PGS"
    
    # Return
    results_list = list("plot" = plot, "overview_table" = table)
    return(results_list)
    
  } else {
    table <- tibble(HI = c( "","", "Pre 65", "", "Post 65", "", "", "Post 65 - Pre 65", "", "", "Post 65 - Pre 65", ""),
                    lowPGS = c(""),
                    midPGS = c(""),
                    highPGS = c(""))
    
    #table[1, 2] <- "BOLD Effect of health shock on smoking probability"
    table[2, 2] <- "Low PGS"
    table[2, 3] <- "Middle PGS"
    table[2, 4] <- "High PGS"
    
    get_coef_str <- function(lin_hyp) {
      return(paste0(round(lin_hyp$test$coef, digits = 3), signif_stars(lin_hyp$test$pvalues[1])))
    }
    
    table[3, 2] <- get_coef_str(glht_low_pre65)
    table[4, 2] <- paste0("(", round(glht_low_pre65$test$sigma, digits = 3), ")")
    table[3, 3] <- get_coef_str(glht_mid_pre65)
    table[4, 3] <- paste0("(", round(glht_mid_pre65$test$sigma, digits = 3), ")")
    table[3, 4] <- get_coef_str(glht_high_pre65)
    table[4, 4] <- paste0("(", round(glht_high_pre65$test$sigma, digits = 3), ")")
    
    table[5, 2] <- get_coef_str(glht_low_post65)
    table[6, 2] <- paste0("(", round(glht_low_post65$test$sigma, digits = 3), ")")
    table[5, 3] <- get_coef_str(glht_mid_post65)
    table[6, 3] <- paste0("(", round(glht_mid_post65$test$sigma, digits = 3), ")")
    table[5, 4] <- get_coef_str(glht_high_post65)
    table[6, 4] <- paste0("(", round(glht_high_post65$test$sigma, digits = 3), ")")
    
    #table[7, 2] <- "BOLD Effect of health insurance on effect of health shock"
    table[7, 2] <- "Low PGS"
    table[7, 3] <- "Middle PGS"
    table[7, 4] <- "High PGS"
    
    table[8, 2] <- get_coef_str(glht_prelow_vs_postlow)
    table[9, 2] <- paste0("(", round(glht_prelow_vs_postlow$test$sigma, digits = 3), ")")
    table[8, 3] <- get_coef_str(glht_premid_vs_postmid)
    table[9, 3] <- paste0("(", round(glht_premid_vs_postmid$test$sigma, digits = 3), ")")
    table[8, 4] <- get_coef_str(glht_prehigh_vs_posthigh)
    table[9, 4] <- paste0("(", round(glht_prehigh_vs_posthigh$test$sigma, digits = 3), ")")
    
    #table[11, 2] <- "BOLD Differential effect of health insurance by genetic group"
    table[10, 2] <- "High PGS - low PGS"
    table[10, 3] <- "Middle PGS - low PGS"
    table[10, 4] <- "High PGS - middle PGS"
    
    
    table[11, 2] <- get_coef_str(glht_slopehigh_vs_slopelow)
    table[12, 2] <- paste0("(", round(glht_slopehigh_vs_slopelow$test$sigma, digits = 3), ")")
    table[11, 3] <- get_coef_str(glht_slopemid_vs_slopelow)
    table[12, 3] <- paste0("(", round(glht_slopemid_vs_slopelow$test$sigma, digits = 3), ")")
    table[11, 4] <- get_coef_str(glht_slopehigh_vs_slopemid)
    table[12, 4] <- paste0("(", round(glht_slopehigh_vs_slopemid$test$sigma, digits = 3), ")")
    
    names(table)[1] <- ""
    names(table)[2] <- ""
    names(table)[3] <- ""
    names(table)[4] <- ""
    
    # Return
    results_list = list("plot" = plot, "overview_table" = table)
    return(results_list)
    
  }
}


## ADDED BY PIA 
# Function: Returning descriptive table for 2 diferent group defined by "vargroup"

descriptive_table_var <- function(df, name, vargroup, group1name, group2name ){
  
  # Create skeleton of table ------------------------------------------------
  table <- tibble(variable = c("", "Age (baseline)", "Smoking PGS", "Years of education", 
                               "Income (nominal \\$ 1000)", "No. waves present",
                               "", "Female", "Smoking (baseline)", "Persistently uninsured", 
                               "CV health shock",
                               "Avg. cessation rate (baseline smokers)", "No. of individuals", 
                               "No. of person-year individuals"),
                  All = "",
                  Uninsured = "",
                  Insured = "" )
  # Create baseline datasets ------------------------------------------------
  # "All"
  baseline_all <- df %>% group_by(hhidpn) %>%
    arrange(year) %>%
    mutate(no_waves = n(),
           max_cv = max(cv)) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    select(hhidpn, age, si_1, no_waves, female,
           smoken, uninsured, max_cv, unins,
           raedyrs, iearn) %>%
    mutate(iearn = iearn/1000)
  
  # "Uninsured"
  baseline_u <- baseline_all %>%
    filter(UQ(sym(vargroup))   == 1)
  
  # "Insured"
  baseline_i <- baseline_all %>%
    filter(UQ(sym(vargroup) )  == 0)
  
  #All sample for n of person year individuals
  dfAll <- df
  dfuns <- df %>% filter(unins == 0)
  dfins <- df %>% filter(unins == 1)
  # Fill columns ------------------------------------------------------------
  groups <- list(NA, baseline_all, baseline_u, baseline_i)
  groupsU <- list(NA, dfAll, dfuns, dfins)
  
  for(i in 2:4){
    # In order of columns
    # Select data for this column
    baseline_df <- groups[[i]]
    dfU <- groupsU[[i]]
    
    # Mean and standard deviations --------------------------------------------
    # Age (baseline)
    sum_age <- paste_mean_std(baseline_df, "age")
    
    # PGS
    sum_pgs <- paste_mean_std(baseline_df, "si_1")
    
    # Years of education
    # This variable contains missing variables!
    sum_edyrs <- paste_mean_std(baseline_df %>%
                                  filter(!is.na(raedyrs)), "raedyrs")
    
    # Individual income (nominal dollars)
    sum_earnings <- paste_mean_std(baseline_df, "iearn")
    
    # No. waves present
    sum_no_waves <- paste_mean_std(baseline_df, "no_waves")
    
    # Mean in percent ---------------------------------------------------------
    # Female
    sum_female <- mean_perc(baseline_df, "female")
    
    # Smoking (baseline)
    sum_smoking <- mean_perc(baseline_df, "smoken")
    
    # Persistently uninsured
    sum_unins <- mean_perc(baseline_df, "unins")
    
    # CV health shock
    sum_cv <- mean_perc(baseline_df, "max_cv")
    
    # No. of individuals
    n_individuals <- count_individuals(baseline_df, simple = TRUE)
    
    #No of person - year individulas
    n_personyearU <- dfU %>% count()
    
    # Add calculated statistics to table
    # CAUTION: If a variable is added, moved, or deleted in the following line,
    # the row numbers in the "average cessation rate" just below might also needs to be updated
    column <- c("Mean (SD)", sum_age, sum_pgs, sum_edyrs, sum_earnings, sum_no_waves, "\\%",
                sum_female, sum_smoking, sum_unins,
                sum_cv, " ", n_individuals, n_personyearU)
    table[[i]] <- column
  }
  # Average cessation rate among baseline smokers -----------------------------------------
  # Calculated and added outside of above loop
  
  # COLUMN "All"
  #table[12, 2] <- mean_cess_baseline_smk(df, baseline_all)
  table$All[12] <- mean_cess_baseline_smk(df, baseline_all)
  
  #table1 <- as.data.frame(table)
  #all_cess <- mean_cess_baseline_smk(df, baseline_all)
  #table1[12, 2] <- all_cess
  # COLUMN "Uninsured"
  uninsured_cess <- mean_cess_baseline_smk(df, baseline_u, return_df = TRUE)
  table$Uninsured[12] <- uninsured_cess$rate
  
  
  # COLUMN "Insured"
  insured_cess <- mean_cess_baseline_smk(df, baseline_i, return_df = TRUE)
  table$Insured[12] <- insured_cess$rate
  
  
  # Change column names for table
  names(table)[1] <- "BOLD "
  names(table)[2] <- "BOLD All"
  names(table)[3] <- group1name
  names(table)[4] <- group2name
  # Adding the P values ---------------------------------------------
  
  # T tests
  test_vars <- c("age", "si_1", "raedyrs", "iearn", "no_waves", "female",
                 "smoken", "unins", "max_cv")
  p_2 <- c()
  for (v in test_vars) {
    if (v == vargroup) {
      p_temp <- "-"
    } else {
      p_temp <- format(round(perform_ttest(baseline_u, baseline_i, v), digits = 2), nsmall = 2)
    }
    p_2 <- c(p_2, p_temp)
  }
  p_2 <- c("",p_2[1:5],"",p_2[6:9], round(t.test(uninsured_cess$df$cess, insured_cess$df$cess)$p.value, digits = 2),"","")
  
  # Add to table
  table <-cbind(table,p_2)
  names(table)[5] <- "BOLD P value"
  
  return(table)
}
