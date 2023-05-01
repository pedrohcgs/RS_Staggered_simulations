# Codes to run the simulations based on Wood et al application
# This script runs all simulations in Azure, and store all results into
# the Temp folder
#----------------------------------------------------------------------
#Load packages
library(here)
library(dplyr)
library(ggplot2)
library(purrr)
library(doAzureParallel)
library(staggered)
library(glue)
#----------------------------------------------------------------------
## Azure setup 
# Make sure you setup the credential correctly before running this!
setCredentials("credentials2_anonymous.json")
#cluster <- makeCluster("cluster_small.json")
cluster <- makeCluster("Replication/cluster.json")
registerDoAzureParallel(cluster)
getDoParWorkers()
setVerbose(TRUE)
setAutoDeleteJob(FALSE)
#----------------------------------------------------------------------
# Simulation function that submits Azure jobs
runSims <- function(baseDir = "Temp", 
                    outcome = "complaints", 
                    numSims = 1000, 
                    skipDiDPackage = TRUE, 
                    collapseToAnnual = FALSE, 
                    addHeterogeneity = FALSE,
                    estimators = c("efficient-longX", 
                                   "efficient", 
                                   "SDIM",
                                   "CS2",
                                   "SA"),
                    default_num_fisher_permutations = 500,
                    errorhandling = "pass",
                    repairSims = FALSE){
  
  #This function submits Azure jobs to run simulations
  #Basedir = directory where results will be saved
  #outcome = outcomes used in simulations
  #numSims = numbe of sims
  #skipDiDPackage = if true, uses our own computation of Callaway &SantAnna rather than Pedro's package
  #collapseToAnnual = if true, collapse data to annual level from monthly level
  #addHeterogeneity = add heterogeneous treatment effects
  # estimators = list of which estimators to use
  
  #----------------------------------------------------------------------
  #This sub-function simulates the data for a given permutation and see
  #first_trained table = table with (i,g) pairs from the data
  #outcome_table = panel dataframe with the (i,t,y)
  simulate_data_fn <- function(seed,
                               first_trained_table, 
                               outcome_table, 
                               N_trained){
    set.seed(seed)
    
    #Pick N_trained officer ids (w/o replacement) from the user ids of the untrained
    #We will then match these officer ids to the trained officers
    first_trained_table$simID <- sample(x = unique(outcome_table$uid),
                                        size = N_trained,
                                        replace = F )
    sim_data <- left_join(outcome_table,
                          first_trained_table, 
                          by = c("uid" = "simID"))
    
    sim_data <- sim_data %>% 
      mutate(first_trained = ifelse(is.na(first_trained),
                                    Inf,
                                    first_trained)) #set first_trained to inf for those without a match
    sim_data <-sim_data %>% 
      filter(t_min <= period & period <= t_max)
    
    #If add heterogeneity, make Y = Y(infty) + 1[period>=first_trained]*u*sd, where sd is the sd of the outcome
    if(addHeterogeneity){
      sim_data <- 
        sim_data %>%
        ungroup() %>%
        mutate(complaints = complaints +
                 u*sd(complaints)*(period >= first_trained)) %>%
        mutate(sustained = sustained +
                 u*sd(sustained)*(period >= first_trained)) %>%
        mutate(force = force + 
                 u*sd(force)*(period >= first_trained))
    }
    
    return(sim_data)
  }
  #----------------------------------------------------------------------
  #Function used to call staggered for a given simulation
  compute_efficient_estimator_for_simulation <-
    function(seed,
             df, 
             compute_fisher = TRUE, 
             num_fisher_permutations = default_num_fisher_permutations,
             ...){
      #df <- simulate_data_fn(seed)
      df <- df %>% rename(t = period,
                          y = !!outcome,
                          g = first_trained,
                          i = uid)
      
      simple_results <- staggered(df = df, 
                                  estimand = "simple",
                                  compute_fisher = compute_fisher, 
                                  num_fisher_permutations = num_fisher_permutations,
                                  ...) %>% 
        mutate(estimand = "simple")
      calendar_results <- staggered(df = df, 
                                    estimand = "calendar",
                                    compute_fisher = compute_fisher, 
                                    num_fisher_permutations = num_fisher_permutations,
                                    ...) %>% 
        mutate(estimand = "calendar")
      cohort_results <- staggered(df = df,
                                  estimand = "cohort",
                                  compute_fisher = compute_fisher, 
                                  num_fisher_permutations = num_fisher_permutations,
                                  ...) %>% 
        mutate(estimand = "cohort")
      eventtime_results <- staggered(df = df,
                                     estimand = "eventstudy",
                                     eventTime = 0, 
                                     compute_fisher = compute_fisher, 
                                     num_fisher_permutations = num_fisher_permutations,
                                     ...) %>% 
        mutate(estimand = "ES0")
      
      results <- bind_rows(simple_results, 
                           calendar_results, 
                           cohort_results, 
                           eventtime_results) %>% 
        mutate(seed = seed)
      
      return(results)
    }
  #----------------------------------------------------------------------
  ## Prepare the data ##
  
  #Create the base directory if it doesn't exist (suppress warnings if it exists)  
  dir.create(here(baseDir),
             showWarnings = FALSE, 
             recursive = TRUE)
  
  load(here("Data/procedural_justice_revisited/products/rdata/3_officer_level_data.RData"))
  
  ##Remove special units from the data
  assignment <- readRDS(here("Data/procedural_justice-master/products/rdata/assignment.rds"))
  pj_officer_level_balanced <-
    pj_officer_level_balanced %>%
    left_join(assignment %>%
                select(uid,type,unit),
              by = "uid") %>%
    filter(!is.na(type) & type != "special") %>%
    filter(first_trained  > 14) #exclude the pilot as well
  
  
  #Create a dataframe with the marginal distribution of when officers were first trained
  first_trained_table <- pj_officer_level_balanced %>% 
    filter(period==1) %>% 
    select(first_trained) %>% 
    group_by()
  
  #Create a dataset with the marginal distribution of outcomes from the data (without first trained date)
  outcome_table <- pj_officer_level_balanced %>% 
    select(-first_trained)
  N_trained <- nrow(first_trained_table) #Number of officers trained
  
  #Create a fixed mean-zero officer level shock to be used for heterogeneous treatment effects (if specified)
  set.seed(0)
  zdraws <- rnorm(n = nrow(first_trained_table))
  zdraws <- zdraws - mean(zdraws) #center to be mean zero
  shock_table <- data.frame(uid = unique(outcome_table$uid), 
                            u = zdraws)
  outcome_table <- left_join(outcome_table, 
                             shock_table, 
                             by = "uid")
  
  
  #If collapseToAnnual = TRUE, collapse table so that months 1-12 correspond with t=1, months 13-24 correspond with t=2, etc
  if(collapseToAnnual == TRUE){
    first_trained_table <- first_trained_table %>% 
      mutate(first_trained = floor((first_trained-1)/12)+1) ##modify to cohort year
    outcome_table <- outcome_table %>% 
      mutate(period = floor((period-1)/12)+1) %>%
      group_by(uid,period) %>% 
      summarise(complaints = sum(complaints),
                force = sum(force), 
                sustained = sum(sustained), 
                u = sum(u)) ##modify to cohort year
  }
  
  t_min <- 1
  t_max <- Inf
  
  #Set foreach to not wait after submitting a job
  foreach_opts <- list(wait = FALSE)
  
  
  ##Create lists of the seeds to run
  #If repair_sims = TRUE, we only runs the seeds saved in missing
  #If repair_sims = FALSE, we run 1:numSims
  
  if(repairSims == TRUE){
    seedsLongX <- readRDS(here(glue("{baseDir}/wood-et-al-based-sims-efficient-T{t_min}-{t_max}-missing-seeds.rds")))
    seedsDiDA0 <- readRDS(here(glue("{baseDir}/wood-et-al-based-sims-efficient-DiDA0-T{t_min}-{t_max}-missing-seeds.rds")))
    seedsSDIM <- readRDS(here(glue("{baseDir}/wood-et-al-based-sims-SDIM-T{t_min}-{t_max}-missing-seeds.rds")))
    seedsCS2 <- readRDS(here(glue("{baseDir}/wood-et-al-based-sims-CS2-T{t_min}-{t_max}-missing-seeds.rds")))
    seedsSA <- readRDS(here(glue("{baseDir}/wood-et-al-based-sims-SA-T{t_min}-{t_max}-missing-seeds.rds")))
    #seedsCS <- readRDS(here(glue("{baseDir}/wood-et-al-based-sims-CS-T{t_min}-{t_max}-missing-seeds.rds")))
    repairSuffix <- "-repair" #this will be added to the filenames
  }else{
    seedsLongX <- 1:numSims
    seedsDiDA0 <- 1:numSims
    seedsSDIM <- 1:numSims
    seedsCS2 <- 1:numSims
    seedsSA <- 1:numSims
    seedsCS <- 1:numSims
    repairSuffix <- ""
  }
  
  ## Simulation results for efficient estimator (using full A0)----
  
  if("efficient-longX" %in% estimators & length(seedsLongX) > 0){
    simResults <- foreach(seed = seedsLongX, 
                          .options.azure = foreach_opts,
                          .errorhandling = errorhandling ) %dopar% {
                            library(dplyr)
                            library(purrr)
                            library(staggered)
                            
                            df <- simulate_data_fn(seed,
                                                   first_trained_table, 
                                                   outcome_table, 
                                                   N_trained)
                            compute_efficient_estimator_for_simulation(seed,df, 
                                                                       use_DiD_A0 = FALSE,
                                                                       num_fisher_permutations = 50)
                          }
    
    #simResults <- reduce(.x = simResults, .f = bind_rows)
    #saveRDS(simResults, here(glue("{baseDir}/wood-et-al-based-sims-efficient-T{t_min}-{t_max}.rds")))
    saveRDS(simResults, 
            here(glue("{baseDir}/wood-et-al-based-sims-efficient-T{t_min}-{t_max}-jobID{repairSuffix}.rds")))
  }
  
  ## Simulation results for efficient estimator using only the DiD A0 ----
  if("efficient" %in% estimators & length(seedsDiDA0) > 0){
    simResultsDiDA0 <- foreach(seed = seedsDiDA0, 
                               .options.azure = foreach_opts, 
                               .errorhandling = errorhandling ) %dopar% {
      library(dplyr)
      library(purrr)
      library(staggered)
      #source(here::here("Code/calculation_functions.R"))
      df <- simulate_data_fn(seed,first_trained_table, 
                             outcome_table,
                             N_trained)
      
      compute_efficient_estimator_for_simulation(seed,
                                                 df, 
                                                 use_DiD_A0 =T)
    }
    
    #simResultsDiDA0 <- reduce(.x = simResultsDiDA0, .f = bind_rows)
    #saveRDS(simResultsDiDA0, here(glue("{baseDir}/wood-et-al-based-sims-efficient-DiDA0-T{t_min}-{t_max}.rds")))
    saveRDS(simResultsDiDA0, 
            here(glue("{baseDir}/wood-et-al-based-sims-efficient-DiDA0-T{t_min}-{t_max}-jobID{repairSuffix}.rds")))
  }
  
  ## Simulation results for SDIM estimator ----
  
  if("SDIM" %in% estimators & length(seedsSDIM) > 0){
    simResultsSDIM <- foreach(seed = seedsSDIM,
                              .options.azure = foreach_opts,
                              .errorhandling = errorhandling ) %dopar% {
      library(dplyr)
      library(purrr)
      library(staggered)
      #source(here::here("Code/calculation_functions.R"))
      df <- simulate_data_fn(seed,
                             first_trained_table, 
                             outcome_table,
                             N_trained)
      compute_efficient_estimator_for_simulation(seed,
                                                 df, 
                                                 beta =0) %>% 
        mutate(betaType = "SDIM")
    }
    
    #simResultsSDIM <- reduce(.x = simResultsSDIM, .f = bind_rows)
    #saveRDS(simResultsSDIM, here(glue("{baseDir}/wood-et-al-based-sims-SDIM-T{t_min}-{t_max}.rds")))
    saveRDS(simResultsSDIM, 
            here(glue("{baseDir}/wood-et-al-based-sims-SDIM-T{t_min}-{t_max}-jobID{repairSuffix}.rds")))
    
  }
  
  ## Our version of CS ----
  if("CS2" %in% estimators & length(seedsCS2) > 0 ){
    simResultsCS2 <- foreach(seed = seedsCS2,
                             .options.azure = foreach_opts,
                             .errorhandling = errorhandling ) %dopar% {
      library(dplyr)
      library(purrr)
      library(staggered)
      #source(here::here("Code/calculation_functions.R"))
      df <- simulate_data_fn(seed,
                             first_trained_table, 
                             outcome_table,
                             N_trained)
      compute_efficient_estimator_for_simulation(seed,
                                                 df, 
                                                 beta =1, 
                                                 use_DiD_A0=T) %>%
        mutate(betaType = "CS2")
    }
    
    # simResultsCS2 <- reduce(.x = simResultsCS2, .f = bind_rows)
    # saveRDS(simResultsCS2, here(glue("{baseDir}/wood-et-al-based-sims-CS2-T{t_min}-{t_max}.rds")))
    saveRDS(simResultsCS2, 
            here(glue("{baseDir}/wood-et-al-based-sims-CS2-T{t_min}-{t_max}-jobID{repairSuffix}.rds")))
  }
  
  
  ## Our version of Sun & Abraham
  if("SA" %in% estimators & length(seedsSA) > 0){
    simResultsSA <- foreach(seed = seedsSA,
                            .options.azure = foreach_opts,
                            .errorhandling = errorhandling ) %dopar% {
      library(dplyr)
      library(purrr)
      library(staggered)
      
      df <- simulate_data_fn(seed,
                             first_trained_table,
                             outcome_table, 
                             N_trained)
      compute_efficient_estimator_for_simulation(seed, 
                                                 df, 
                                                 beta =1, 
                                                 use_DiD_A0= TRUE, 
                                                 use_last_treated_only = TRUE) %>%
        mutate(betaType = "SA")
    }
    
    # simResultsSA <- reduce(.x = simResultsSA, .f = bind_rows)
    # saveRDS(simResultsSA, here(glue("{baseDir}/wood-et-al-based-sims-SA-T{t_min}-{t_max}.rds")))
    saveRDS(simResultsSA, 
            here(glue("{baseDir}/wood-et-al-based-sims-SA-T{t_min}-{t_max}-jobID{repairSuffix}.rds")))
    
  }
  
  ## Use oracle betastar ----
  #Note that this assumes the sharp null so would need to be updated if you want to do this for heterogeneous effects
  #Note also that this uses the oracle for the *simple* estimand
  
  # df <- simulate_data_fn(1,first_trained_table, outcome_table, N_trained)
  # df <- df %>% rename(t = period, y = !!outcome, g = first_trained, i = uid)
  # 
  # df_wide <-
  #   df  %>% 
  #   reshape2::dcast(i ~ t, value.var = "y") %>%
  #   dplyr::select(-i)
  # 
  # sigma_true <- cov(as.matrix(df_wide))
  # g_level_summaries <- compute_g_level_summaries(df, refine_S_g = FALSE)
  # S_g_true <- rep(list(sigma_true), length(g_level_summaries$S_g_List))
  # A_theta_list <- staggered::create_Atheta_list_for_simple_average_ATE(g_level_summaries$g_list,g_level_summaries$t_list, g_level_summaries$N_g_List)
  # A_0_list <- staggered::create_A0_list(g_level_summaries$g_list,g_level_summaries$t_list)
  # betastar_oracle <- compute_Betastar(Ybar_g_list = g_level_summaries$Ybar_g_List,A_theta_list = A_theta_list, A_0_list = A_0_list,S_g_list = S_g_true,N_g_list = g_level_summaries$N_g_List)
  # betastar_df <- compute_Betastar(Ybar_g_list = g_level_summaries$Ybar_g_List,A_theta_list = A_theta_list, A_0_list = A_0_list,S_g_list = g_level_summaries$S_g_List,N_g_list = g_level_summaries$N_g_List)
  # 
  # 
  # simResultsOracle <- foreach(seed = 1:numSims) %dopar% {
  #   library(dplyr)
  #   library(purrr)
  #   library(staggered)
  #   #source(here::here("Code/calculation_functions.R"))
  #   df <- simulate_data_fn(seed,first_trained_table, outcome_table, N_trained)
  # 
  #   compute_efficient_estimator_for_simulation(seed,df, use_DiD_A0 = FALSE, beta = betastar_oracle)
  # }
  # 
  # simResultsOracle <- reduce(.x = simResultsOracle, .f = bind_rows)
  # 
  # saveRDS(simResultsOracle, here(glue("{baseDir}/wood-et-al-based-sims-efficient-oracle-T{t_min}-{t_max}.rds")))
  
  
  
  ## Callaway & sant'anna sims using DiD package ----
  
  if(!skipDiDPackage){
    simResultsCS <- foreach(seed = seedsCS,
                            .options.azure = foreach_opts,
                            .errorhandling = errorhandling ) %dopar% {
      library(dplyr)
      library(purrr)
      library(did)
      #library(tictoc)
      
      df <- simulate_data_fn(seed,
                             first_trained_table, 
                             outcome_table, 
                             N_trained)
      df <- df %>% rename(t = period,
                          y = !!outcome, 
                          g = first_trained, 
                          i = uid)
      df$g[df$g==Inf] <- 0 #Modify never treated to g=0 to comply with Pedro's silly numbering convention
      
      #tic()
      cs_results <-
        att_gt(yname="y",
               gname="g",
               idname="i",
               tname="t",
               xformla=~1,
               data=df,
               est_method="reg",
               print_details=FALSE,
               control_group = "notyettreated",
               bstrap =F
        )
      simple_results <- data.frame(estimate = aggte(cs_results,
                                                    type = "simple")$overall.att,
                                   se = aggte(cs_results,
                                              type = "simple")$overall.se)
      calendar_results <- data.frame(estimate = aggte(cs_results,
                                                      type = "calendar")$overall.att,
                                     se = aggte(cs_results,
                                                type = "calendar")$overall.se)
      group_results <- data.frame(estimate = aggte(cs_results,
                                                   type = "group")$overall.att,
                                  se = aggte(cs_results,
                                             type = "group")$overall.se)
      
      
      dynamic_object <- aggte(cs_results, 
                              type = "dynamic")
      dynamic_results <- data.frame(estimate = dynamic_object$att.egt[dynamic_object$egt == 0],
                                    se = dynamic_object$se.egt[dynamic_object$egt == 0])
      
      bind_rows(simple_results %>% 
                  mutate(estimand = "simple"),
                calendar_results %>% 
                  mutate(estimand = "calendar"),
                group_results %>% 
                  mutate(estimand = "group"),
                dynamic_results %>% 
                  mutate(estimand = "ES0")) %>%
        mutate(seed = seed)
      
      #toc()
      
    }
    
    # simResultsCS <- reduce(.x = simResultsCS, .f = bind_rows)
    # saveRDS(simResultsCS, here(glue("{baseDir}/wood-et-al-based-sims-CS-T{t_min}-{t_max}.rds")))
    saveRDS(simResultsCS, 
            here(glue("{baseDir}/wood-et-al-based-sims-CS-T{t_min}-{t_max}-jobID{repairSuffix}.rds")))
  }
  
}
#----------------------------------------------------------------------
# Actually run the simulations


runSims(baseDir = "Temp/Wood-et-al-sims/force", 
        numSims = 1000, 
        outcome = "force")
runSims(baseDir = "Temp/Wood-et-al-sims/complaints",
        numSims = 1000, 
        outcome = "complaints")
runSims(baseDir = "Temp/Wood-et-al-sims/sustained", 
        numSims = 1000,
        outcome = "sustained")

runSims(baseDir = "Temp/Wood-et-al-sims/force-annual-t-and-g", 
        numSims = 1000, 
        outcome = "force", 
        collapseToAnnual = TRUE)
runSims(baseDir = "Temp/Wood-et-al-sims/complaints-annual-t-and-g", 
        numSims = 1000, 
        outcome = "complaints", 
        collapseToAnnual = TRUE)
runSims(baseDir = "Temp/Wood-et-al-sims/sustained-annual-t-and-g",
        numSims = 1000, 
        outcome = "sustained",
        collapseToAnnual = TRUE)

runSims(baseDir = "Temp/Wood-et-al-sims/complaints-het-effects", 
        numSims = 1000, 
        outcome = "complaints",
        addHeterogeneity = TRUE)

#stopCluster(cluster)
#----------------------------------------------------------------------
# Function to load sim results from Azure

loadSims <- function(baseDir = "Temp", 
                     outcome = "complaints", 
                     numSims = 1000, 
                     skipDiDPackage = TRUE, 
                     collapseToAnnual = FALSE, 
                     addHeterogeneity = FALSE,
                     estimators = c("efficient-longX", 
                                    "efficient", 
                                    "SDIM", 
                                    "CS2", 
                                    "SA"),
                     ...){
  
  #This function loads the results from Azure from the runSims() run with the same arguments
  
  
  t_min <- 1
  t_max <- Inf
  
  #Create a fn that converts errors to warnings
  errorHandler = function(e){warning(e)}
  
  #Function for binding rows if df's, where we set each input to an empty df if it's not a data.frame
  #This is used for combining results from Azure, where we get an error message in the list instead of a df if there was an Azure error
  bind_rows_if_dfs <- function(df1,df2){
    if(!is.data.frame(df1)){df1 <- data.frame()}
    if(!is.data.frame(df2)){df2 <- data.frame()}
    
    return(bind_rows(df1,df2))
  }
  
  ## Simulation results for efficient estimator (using full A0)----
  
  if("efficient-longX" %in% estimators){
    tryCatch({
      simResults <- getJobResult(readRDS(here(glue("{baseDir}/wood-et-al-based-sims-efficient-T{t_min}-{t_max}-jobID.rds"))))
      simResults <- reduce(.x = simResults, 
                           .f = bind_rows_if_dfs)
      saveRDS(simResults,
              here(glue("{baseDir}/wood-et-al-based-sims-efficient-T{t_min}-{t_max}.rds")))
    }, error = errorHandler)
    
    missingSeeds<- which(! (1:numSims) %in% unique(simResults$seed) )
    saveRDS(missingSeeds, 
            here(glue("{baseDir}/wood-et-al-based-sims-efficient-T{t_min}-{t_max}-missing-seeds.rds")))
  }
  
  ## Simulation results for efficient estimator using only the DiD A0 ----
  if("efficient" %in% estimators){
    
    tryCatch({
      simResultsDiDA0 <-getJobResult(readRDS(here(glue("{baseDir}/wood-et-al-based-sims-efficient-DiDA0-T{t_min}-{t_max}-jobID.rds"))))
      simResultsDiDA0 <- reduce(.x = simResultsDiDA0,
                                .f = bind_rows_if_dfs)
      saveRDS(simResultsDiDA0,
              here(glue("{baseDir}/wood-et-al-based-sims-efficient-DiDA0-T{t_min}-{t_max}.rds")))
    }, error = errorHandler)
    missingSeedsDiDA0 <- which(! (1:numSims) %in% unique(simResultsDiDA0$seed) )
    saveRDS(missingSeedsDiDA0, 
            here(glue("{baseDir}/wood-et-al-based-sims-efficient-DiDA0-T{t_min}-{t_max}-missing-seeds.rds")))
  }
  
  ## Simulation results for SDIM estimator ----
  
  if("SDIM" %in% estimators){
    tryCatch({
      simResultsSDIM <- getJobResult(readRDS(here(glue("{baseDir}/wood-et-al-based-sims-SDIM-T{t_min}-{t_max}-jobID.rds"))))
      simResultsSDIM <- reduce(.x = simResultsSDIM, 
                               .f = bind_rows_if_dfs)
      saveRDS(simResultsSDIM,
              here(glue("{baseDir}/wood-et-al-based-sims-SDIM-T{t_min}-{t_max}.rds")))
      
      #save the missing indices (from Azure errors)
      missingSeedsSDIM <- which(! (1:numSims) %in% unique(simResultsSDIM$seed) )
      saveRDS(missingSeedsSDIM,
              here(glue("{baseDir}/wood-et-al-based-sims-SDIM-T{t_min}-{t_max}-missing-seeds.rds")))
      
    }, error = errorHandler)
  }
  ## Our version of CS ----
  if("CS2" %in% estimators){
    tryCatch({
      simResultsCS2 <- getJobResult(readRDS(here(glue("{baseDir}/wood-et-al-based-sims-CS2-T{t_min}-{t_max}-jobID.rds"))))
      simResultsCS2 <- reduce(.x = simResultsCS2, 
                              .f = bind_rows_if_dfs)
      saveRDS(simResultsCS2,
              here(glue("{baseDir}/wood-et-al-based-sims-CS2-T{t_min}-{t_max}.rds")))
    }, error = errorHandler)
    
    missingSeedsCS2 <- which(! (1:numSims) %in% unique(simResultsCS2$seed) )
    saveRDS(missingSeedsCS2, 
            here(glue("{baseDir}/wood-et-al-based-sims-CS2-T{t_min}-{t_max}-missing-seeds.rds")))
    
  }
  
  
  ## Our version of Sun & Abraham
  if("SA" %in% estimators){
    tryCatch({
      simResultsSA <- getJobResult(readRDS(here(glue("{baseDir}/wood-et-al-based-sims-SA-T{t_min}-{t_max}-jobID.rds"))))
      simResultsSA <- reduce(.x = simResultsSA, 
                             .f = bind_rows_if_dfs)
      saveRDS(simResultsSA,
              here(glue("{baseDir}/wood-et-al-based-sims-SA-T{t_min}-{t_max}.rds")))
    }, error = errorHandler)
    
    missingSeedsSA <- which(! (1:numSims) %in% unique(simResultsSA$seed) )
    saveRDS(missingSeedsSA,
            here(glue("{baseDir}/wood-et-al-based-sims-SA-T{t_min}-{t_max}-missing-seeds.rds")))
  }
  
  
  ## Callaway & sant'anna sims using DiD package ----
  
  if(!skipDiDPackage){
    tryCatch({
      simResultsCS <- getJobResult(readRDS(here(glue("{baseDir}/wood-et-al-based-sims-CS-T{t_min}-{t_max}-jobID.rds"))))
      simResultsCS <- reduce(.x = simResultsCS, 
                             .f = bind_rows_if_dfs)
      saveRDS(simResultsCS,
              here(glue("{baseDir}/wood-et-al-based-sims-CS-T{t_min}-{t_max}.rds")))
    }, error = errorHandler)
    
    missingSeedsCS <- which(! (1:numSims) %in% unique(simResultsCS$seed) )
    saveRDS(missingSeedsCS,
            here(glue("{baseDir}/wood-et-al-based-sims-CS-T{t_min}-{t_max}-missing-seeds.rds")))
  }
  
}
#----------------------------------------------------------------------
#Actually load the simulation results

loadSims(baseDir = "Temp/Wood-et-al-sims/force",
         numSims = 1000,
         outcome = "force")
loadSims(baseDir = "Temp/Wood-et-al-sims/complaints", 
         numSims = 1000, 
         outcome = "complaints")
loadSims(baseDir = "Temp/Wood-et-al-sims/sustained",
         numSims = 1000, 
         outcome = "sustained")

loadSims(baseDir = "Temp/Wood-et-al-sims/force-annual-t-and-g",
         numSims = 1000, 
         outcome = "force", 
         collapseToAnnual = TRUE)
loadSims(baseDir = "Temp/Wood-et-al-sims/complaints-annual-t-and-g",
         numSims = 1000,
         outcome = "complaints", 
         collapseToAnnual = TRUE)
loadSims(baseDir = "Temp/Wood-et-al-sims/sustained-annual-t-and-g", 
         numSims = 1000, 
         outcome = "sustained",
         collapseToAnnual = TRUE)

loadSims(baseDir = "Temp/Wood-et-al-sims/complaints-het-effects", 
         numSims = 1000,
         outcome = "complaints",
         addHeterogeneity = TRUE)
#----------------------------------------------------------------------


### Repair sims from azure failures 
# runSims(baseDir = "Temp/Wood-et-al-sims/force", numSims = 1000, outcome = "force", repairSims = TRUE)
# runSims(baseDir = "Temp/Wood-et-al-sims/complaints", numSims = 1000, outcome = "complaints", repairSims = TRUE)
# runSims(baseDir = "Temp/Wood-et-al-sims/sustained", numSims = 1000, outcome = "sustained", repairSims = TRUE)
