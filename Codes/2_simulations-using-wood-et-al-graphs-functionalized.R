# This script analyzes/summarizes all the simulation results based on Wood et al. 
# application that were ran in Azure.
# Do not need to re-run the simulations!
#----------------------------------------------------------------------
#Load packages
library(dplyr)
library(here)
library(stargazer)
library(gt)
library(glue)
#----------------------------------------------------------------------
## More Systemized tables 
#----------------------------------------------------------------------
# Coverage Table
create_coverage_table <- function(simResults, 
                                  CS = FALSE){
  
  #Rename thetahat to estimate if we don't have a varible named estimate
  if(is.null(simResults$estimate)){
    simResults <- simResults %>% 
    rename(estimate = thetahat)}
  
  #Change "group" to "cohort", since they mean the same thing and we use both
  simResults <- simResults %>% 
    mutate(estimand = ifelse(estimand == "group", 
                             "cohort",
                             estimand))
  
  
  coverageResults <-
    simResults %>% 
    mutate(ub = estimate + se * 1.96,
           lb = estimate- se*1.96) %>%
    group_by(estimand) %>%
    summarise(mean = mean(estimate), 
              covers0 = mean(lb <= 0 & 0<=ub), 
              mean_se = mean(se), 
              sd = sd(estimate),
              fisher_size = mean(fisher_pval <= 0.05))
  
  return(coverageResults)
}
#----------------------------------------------------------------------
# Long Coverage Table
create_long_coverage_table <- function(simResultsDiDA0,
                                       simResultsCS, 
                                       simResultsSA = NULL,
                                       rescaleFactor = 100, 
                                       renameColumns = FALSE, 
                                       roundDigits = NULL){
  long_coverage_table <-
    bind_rows(  
      create_coverage_table(simResultsDiDA0, CS = FALSE) %>% 
        mutate(estimator ="Efficient"),
      create_coverage_table(simResultsCS, CS = TRUE) %>% 
        mutate(estimator = "CS") %>%
        mutate(estimand = ifelse(estimand == "group", 
                                 "cohort", 
                                 estimand))
    )
  
  if(!is.null(simResultsSA)){
    long_coverage_table <- bind_rows(long_coverage_table,
                                     create_coverage_table(simResultsSA, 
                                                           CS = TRUE) %>%
                                       mutate(estimator = "SA") %>% 
                                       mutate(estimand = ifelse(estimand == "group",
                                                                "cohort",
                                                                estimand)))
  }  
  #Multiply means and sds/ses by rescale factor
  long_coverage_table <- long_coverage_table %>% 
    mutate_at(c("mean", 
                "mean_se",
                "sd"), ~.*100)
  
  #Order columns more intuitively
  long_coverage_table <- long_coverage_table %>% 
    select(estimator, estimand, mean, covers0, fisher_size, mean_se, sd)
  
  #Rename columns for table formatting if renameColumns = T
  if(renameColumns){
    long_coverage_table <- long_coverage_table %>% 
      rename(Estimator = estimator, 
             Estimand = estimand, 
             Bias = mean, 
             Coverage = covers0,
             `FRT Size` = fisher_size,
             `Mean SE` = mean_se, 
             `SD` =sd)
  }
  
  #Round to roundDigits if it is not null
  if(!is.null(roundDigits)){
    long_coverage_table <- long_coverage_table %>%
      mutate_if(is.numeric, ~round(.,digits = roundDigits))
  }
  return(long_coverage_table)  
}
#----------------------------------------------------------------------
# Table for Ratio of MC standard deviations 

create_ratio_of_sds_table <- function(simResultsDiDA0, 
                                      simResultsCS, 
                                      simResultsSA = NULL){
  long_coverage_table <- create_long_coverage_table(simResultsDiDA0,
                                                    simResultsCS,
                                                    simResultsSA = simResultsSA)
  
  #Wide table comparing them (each row is estimand)
  wide_table <-
    long_coverage_table %>%
    tidyr::pivot_wider(id_cols = estimand,
                       names_from = estimator, 
                       values_from = -c(estimand,estimator)) %>%
    select(estimand, 
           contains("sd")) %>%
    mutate(ratio_CS = sd_CS / sd_Efficient)
  
  if(!is.null(simResultsSA)){
    wide_table <- wide_table %>% 
      mutate(ratio_SA = sd_SA / sd_Efficient)
  }
  
  return(wide_table)
}
#----------------------------------------------------------------------
# Function to make the tables we want

makeGraphs <- function(outcome,
                       annualData = FALSE, 
                       hetEffects = FALSE){
  
  baseDir <- glue("Temp/Wood-et-al-sims/{outcome}")
  annualDataSuffix <- ifelse(annualData,
                             "-annual-t-and-g", 
                             "")
  hetEffectsSuffix <- ifelse(hetEffects, 
                             "-het-effects",
                             "")
  baseDir <- paste0(baseDir,
                    annualDataSuffix)
  
  
  simResultsDiDA0 <- readRDS(here(glue("{baseDir}{hetEffectsSuffix}/wood-et-al-based-sims-efficient-DiDA0-T1-Inf.rds")))
  simResultsCS <- readRDS(here(glue("{baseDir}{hetEffectsSuffix}/wood-et-al-based-sims-CS2-T1-Inf.rds")))
  simResultsLongA0 <- readRDS(here(glue("{baseDir}{hetEffectsSuffix}/wood-et-al-based-sims-efficient-T1-Inf.rds")))
  simResultsSA <- readRDS(here(glue("{baseDir}{hetEffectsSuffix}/wood-et-al-based-sims-SA-T1-Inf.rds")))
  
  long_coverage_table <- create_long_coverage_table(simResultsDiDA0,
                                                    simResultsCS,
                                                    simResultsSA,
                                                    renameColumns = TRUE,
                                                    roundDigits = 2)
  
  long_coverage_table %>% 
    mutate(Estimator =ifelse(Estimator == "Efficient", 
                             "PlugIn", 
                             Estimator)) %>%
    mutate(Estimator = ifelse(Estimator == "CS" & Estimand == "ES0", 
                              "CS/dCDH", 
                              Estimator)) %>% #add dCDH to estimator name for ES0
    gt() %>% gt::gtsave(here(glue("Tables/Wood-et-al-sims/coverage-table-{outcome}{hetEffectsSuffix}{annualDataSuffix}.tex")))
  
  
  create_ratio_of_sds_table(simResultsDiDA0,
                            simResultsCS,
                            simResultsSA) %>%
    select(estimand, contains("ratio")) %>%
    rename(Estimand = estimand, `CS`= ratio_CS, 
           `SA` = ratio_SA) %>%
    gt() %>% fmt_number(c("CS", 
                          "SA"), 
                        decimals =2) %>%
    tab_spanner(label = "Ratio of SD to Plug-In",
                columns = c("CS","SA")) %>% 
    gt::gtsave(here(glue("Tables/Wood-et-al-sims/sd-ratio-table-{outcome}{hetEffectsSuffix}{annualDataSuffix}.tex")))
  
  
  
  #Create table using the full A0
  create_long_coverage_table(simResultsLongA0, 
                             simResultsDiDA0, 
                             renameColumns = TRUE,
                             roundDigits = 2) %>%
    mutate(Estimator = ifelse(Estimator == "Efficient", 
                              "PlugIn - Long X",
                              "PlugIn")) %>%
    gt() %>% gt::gtsave(here(glue("Tables/Wood-et-al-sims/coverage-table-long-A0-{outcome}{hetEffectsSuffix}{annualDataSuffix}.tex")))
  
  #Create a table showing the number of successful sims for each spec
  observations_table <- 
    bind_rows(
      data.frame(estimator = c("Plug-In"), 
                 n = NROW(simResultsDiDA0), 
                 min_fisher_permutations = min(simResultsDiDA0$num_fisher_permutations)  ),
      data.frame(estimator = c("Plug-In - Long X"), 
                 n = NROW(simResultsLongA0),
                 min_fisher_permutations = min(simResultsLongA0$num_fisher_permutations)),
      data.frame(estimator = c("CS"), 
                 n = NROW(simResultsCS),
                 min_fisher_permutations = min(simResultsCS$num_fisher_permutations)),
      data.frame(estimator = c("SA"), 
                 n = NROW(simResultsSA),
                 min_fisher_permutations = min(simResultsSA$num_fisher_permutations))
    )
  
  observations_table$outcome <- outcome
  return(observations_table)
}
#----------------------------------------------------------------------
# Produce the tables and graphs
makeGraphs(outcome = "complaints", 
           annualData = FALSE)
makeGraphs(outcome = "sustained",
           annualData = FALSE)
makeGraphs(outcome = "force", 
           annualData = FALSE)

makeGraphs(outcome = "complaints",
           annualData = TRUE)
makeGraphs(outcome = "sustained",
           annualData = TRUE)
makeGraphs(outcome = "force", 
           annualData = TRUE)

makeGraphs(outcome = "complaints", 
           annualData = FALSE, 
           hetEffects = TRUE)
#----------------------------------------------------------------------