# This script analyzes/summarizes all the simulation results based on 2x2 DiD 
# setup that were ran in Azure.
# Do not need to re-run the simulations!
#----------------------------------------------------------------------
# Load libraries
library(dplyr)
library(tidyr)
library(purrr)
library(lfe)
library(here)
library(gt) 
library(latex2exp)
library(staggered)
library(tidyr)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
## Make tables for the 2x2 simulations

simResults <- readRDS(here("Temp/2period-sim-results.rds"))
resultsTable <-
  simResults %>% 
  select(N_treated, 
         N_control,
         rhoY0, 
         hetParam, 
         betastar, 
         betaType, 
         bias, 
         sdThetahat, 
         coverage_SE,
         fisher_size) %>%
  tidyr::pivot_wider(id_cols = c(N_treated, N_control,rhoY0, hetParam, betastar), 
                     names_from = betaType, 
                     values_from = c(bias,
                                     sdThetahat, 
                                     coverage_SE, 
                                     fisher_size) ) %>%
  arrange(-N_treated,-rhoY0,hetParam)

resultsTable %>%
  select(-betastar) %>%
  gt() %>%
  fmt_number(columns = contains("bias"), decimals = 2) %>%
  fmt_number(columns = contains("sd"), decimals = 2) %>%
  fmt_number(columns = contains("coverage"), decimals = 2) %>%
  fmt_number(columns = contains("fisher"), decimals = 2) %>%
  tab_spanner(label = "Bias", columns = contains("bias")) %>%
  tab_spanner(label = "SD", columns = contains("sd")) %>%
  tab_spanner(label = "Coverage", columns = contains("Coverage")) %>%
  tab_spanner(label = "FRT Size", columns = contains("fisher")) %>%
  cols_label(bias_betaStar = "PlugIn", bias_DiD = "DiD", bias_DiM = "DiM") %>%
  cols_label(sdThetahat_betaStar = "PlugIn", sdThetahat_DiD = "DiD", 
             sdThetahat_DiM = "DiM") %>%
  cols_label(coverage_SE_betaStar = "PlugIn", coverage_SE_DiD = "DiD",
             coverage_SE_DiM = "DiM")%>%
  cols_label(fisher_size_betaStar = "PlugIn", fisher_size_DiD = "DiD",
             fisher_size_DiM = "DiM")%>%
  cols_label(N_treated = "N_1", N_control = "N_0", rhoY0 = "rho",
             hetParam = "gamma") %>%
  gtsave(filename = here("Tables/2-period-sims.tex"))

resultsTable %>% 
  select(N_treated,N_control,rhoY0, hetParam, betastar, contains("sdThetahat") ) %>%
  mutate_at(c("sdThetahat_DiD", "sdThetahat_DiM", "sdThetahat_betaStar"),
            ~./sdThetahat_betaStar) %>%
  gt() %>%
  cols_label(sdThetahat_betaStar = "PlugIn", sdThetahat_DiD = "DiD", 
             sdThetahat_DiM = "DiM") %>%
  tab_spanner(label = "SD Relative to Plug-In", columns = contains("sd")) %>%
  fmt_number(columns = contains("sd"), decimals = 2) %>%
  fmt_number(columns = "betastar", decimals = 2) %>%
  cols_label(N_treated = "N_1", N_control = "N_0", rhoY0 = "rho", 
             hetParam = "gamma") %>%
  gtsave(filename = here("Tables/2-period-sims-sd-ratios.tex"))


