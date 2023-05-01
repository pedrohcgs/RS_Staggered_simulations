library(doAzureParallel)
createPODataset <- function(poSeed, 
                            N_treated,
                            N_control, 
                            rhoY0 = 0, 
                            hetParam = 0){
  
  createObservation <- function(i){
    set.seed(i)
    
    #Draw y0 from a normal distribution with covariance rhoY0
    y0 <- MASS::mvrnorm(n=1, 
                        mu = c(0,0),
                        Sigma = matrix(
                          c(1,rhoY0,rhoY0,1),
                          nrow =2))

    df <- data.frame(i = rep(i,2),
                     t = c(1,2),
                     y0 = y0)
    return(df)
  }
  
  df <- map_dfr(.x = poSeed + 1:(N_treated + N_control), 
                ~createObservation(.x))
  
  #Create y1 so that there are a treatment effects proportional to (y0-y0Mean) within each period
    #Note that this ensures that the ATE within each period is 0
  df <- df %>% group_by(t) %>% 
    mutate(y1 = y0 + hetParam * (y0 - mean(y0)))
  
  return(df) 
}

assignTreatments <- function(assignmentSeed, 
                             df,
                             N_treated, 
                             N_control){
  set.seed(assignmentSeed)
  treatedIndices <- sample(x = 1:(N_treated + N_control),
                           size = N_treated, replace = FALSE)
  
  df$g <- Inf
  df$g[which(df$i %in% treatedIndices)] <- 2
  
  
  #Create observed y 
  df <- df %>% mutate(y = ifelse(g == Inf, 
                                 y0, 
                                 (t<g) * y0 + (t >= g) * y1) )
  
  return(df)
}



simulationWrapperFN <- function(N_treated, 
                                N_control,
                                rhoY0,
                                hetParam,
                                N_sims = 1000){
  POdf <- createPODataset(poSeed = 1, 
                          N_treated = N_treated,
                          N_control = N_control,
                          rhoY0 = rhoY0,
                          hetParam = hetParam)
  
  #Compute the oracle betastar
    #First reshape the data in wide format so each row is an observeration
    wide_df <-
    POdf %>% 
      pivot_wider(id_cols = "i", 
                  values_from = c(y0,y1), 
                  names_from = "t") 
    
    beta1 <- lm(y1_2 ~ y0_1, wide_df)$coefficients["y0_1"]
    beta0 <- lm(y0_2 ~ y0_1, wide_df)$coefficients["y0_1"]
    betastar <- (N_control/(N_control + N_treated))*beta1 + (N_treated/(N_control + N_treated))*beta0
    
  A_theta_list <- list(matrix(c(0,1),nrow =1), 
                       matrix(c(0,-1),nrow =1 ))
  A_0_list <- list(matrix(c(1,0),nrow =1),
                   matrix(c(-1,0),nrow =1 ))
  
  # betaStarResultsDF <- map_dfr(.x = 1:N_sims, .f = ~assignTreatments(assignmentSeed = .x, POdf, N_treated = N_treated, N_control = N_control) %>% 
  #                                staggered::staggered(A_theta_list = A_theta_list, A_0_list = A_0_list, compute_fisher = TRUE) %>% mutate(betaType = "betaStar")  )
  
  betaStarResults <- foreach(seed = 1:N_sims, 
                             .export = c("assignTreatments") ) %dopar% {
    library(dplyr)
    library(purrr)
    library(staggered)
    
    assignTreatments(assignmentSeed = seed,
                     POdf, 
                     N_treated = N_treated,
                     N_control = N_control) %>% 
      staggered::staggered(A_theta_list = A_theta_list,
                           A_0_list = A_0_list,
                           compute_fisher = TRUE) %>% 
      mutate(betaType = "betaStar")  
    
  }
  
  
  DiDResults <- foreach(seed = 1:N_sims, 
                        .export = c("assignTreatments") ) %dopar% {
    library(dplyr)
    library(purrr)
    library(staggered)
    
    assignTreatments(assignmentSeed = seed,
                     POdf, 
                     N_treated = N_treated,
                     N_control = N_control) %>% 
      staggered::staggered(A_theta_list = A_theta_list, 
                           A_0_list = A_0_list,
                           compute_fisher = TRUE,
                           beta = 1) %>% 
      mutate(betaType = "DiD")  
    
  }
  
  DiMResults <- foreach(seed = 1:N_sims, 
                        .export = c("assignTreatments") ) %dopar% {
    library(dplyr)
    library(purrr)
    library(staggered)
    
    assignTreatments(assignmentSeed = seed, 
                     POdf, N_treated = N_treated, 
                     N_control = N_control) %>% 
      staggered::staggered(A_theta_list = A_theta_list,
                           A_0_list = A_0_list,
                           compute_fisher = TRUE,
                           beta = 0) %>% 
      mutate(betaType = "DiM")  
    
  }
  
  betaStarResultsDF <- purrr::reduce(betaStarResults, bind_rows)
  DiDResultsDF <- purrr::reduce(DiDResults, bind_rows)
  DiMResultsDF <- purrr::reduce(DiMResults, bind_rows)
  
  coversTheta <- function(thetahat, 
                          se, 
                          thetaTrue = 0,
                          alpha = 0.05){
    c_alpha <- qnorm(1-alpha/2)
    coversThetaTrue <- abs(thetahat - thetaTrue) < c_alpha * se
    return(coversThetaTrue)
  }
  
  
  combinedResults <- bind_rows(betaStarResultsDF, 
                               DiDResultsDF, 
                               DiMResultsDF)
  
  #This assumes that true theta is zero
  MCSummary <- 
    combinedResults %>% group_by(betaType) %>% 
    summarise(bias = mean(estimate), sdThetahat = sd(estimate), 
              avgSE = mean(se), avgSE_conservative = mean(se_neyman),
              coverage_SE = mean( coversTheta(estimate, se) ),
              coverage_SE_conservative = mean( coversTheta(estimate, se_neyman) ),
              fisher_size = mean(fisher_pval < 0.05),
              fisher_size_neyman = mean(fisher_pval_se_neyman < 0.05))
  
  MCSummary <- MCSummary %>% mutate(N_treated = N_treated,
                                    N_control = N_control, 
                                    rhoY0 = rhoY0, 
                                    hetParam = hetParam, 
                                    betastar = betastar)
  return(MCSummary)
}