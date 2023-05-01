# Codes to run the simulations for the 2x2 Setup in Azure and 
# This script runs all simulations in Azure, and store all results into
# the Temp folder
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
# Load auxiliary function
source(here("Codes/aux_simulation_functions.R"))
#----------------------------------------------------------------------
## Azure setup 
# Make sure you setup the credential correctly before running this!
setCredentials("credentials2_anonymous.json")
#cluster <- makeCluster("cluster_small.json")
cluster <- makeCluster("cluster.json")
registerDoAzureParallel(cluster)
getDoParWorkers()
setVerbose(TRUE)
setAutoDeleteJob(FALSE)

# Set numnber of simulations
N_sims <- 1000

#----------------------------------------------------------------------
# Run the simulations in Azure
simResults <- bind_rows(
  simulationWrapperFN(N_treated = 1000,
                      N_control = 1000,
                      rhoY0 = 0, 
                      hetParam = 0,
                      N_sims = N_sims),
  
  simulationWrapperFN(N_treated = 1000,
                      N_control = 1000, 
                      rhoY0 = 0.5,
                      hetParam = 0,
                      N_sims = N_sims),
  
  simulationWrapperFN(N_treated = 1000,
                      N_control = 1000,
                      rhoY0 = 0.99, 
                      hetParam = 0,
                      N_sims = N_sims),
  
  simulationWrapperFN(N_treated = 1000, 
                      N_control = 1000, 
                      rhoY0 = 0, 
                      hetParam = 0.5,
                      N_sims = N_sims),
  
  simulationWrapperFN(N_treated = 1000,
                      N_control = 1000,
                      rhoY0 = 0.5,
                      hetParam = 0.5, 
                      N_sims = N_sims),
  
  simulationWrapperFN(N_treated = 1000,
                      N_control = 1000, 
                      rhoY0 = 0.99, 
                      hetParam = 0.5, 
                      N_sims = N_sims),
  
  simulationWrapperFN(N_treated = 25, 
                      N_control = 25, 
                      rhoY0 = 0,
                      hetParam = 0, 
                      N_sims = N_sims),
  
  simulationWrapperFN(N_treated = 25, 
                      N_control = 25,
                      rhoY0 = 0.5, 
                      hetParam = 0, 
                      N_sims = N_sims),
  
  simulationWrapperFN(N_treated = 25,
                      N_control = 25, 
                      rhoY0 = 0.99, 
                      hetParam = 0, 
                      N_sims = N_sims),
  
  simulationWrapperFN(N_treated = 25, 
                      N_control = 25,
                      rhoY0 = 0, 
                      hetParam = 0.5,
                      N_sims = N_sims),
  
  simulationWrapperFN(N_treated = 25, 
                      N_control = 25, 
                      rhoY0 = 0.5,
                      hetParam = 0.5, 
                      N_sims = N_sims),
  
  simulationWrapperFN(N_treated = 25,
                      N_control = 25,
                      rhoY0 = 0.99,
                      hetParam = 0.5,
                      N_sims = N_sims)
)
#----------------------------------------------------------------------
# Save all results
saveRDS(simResults, file= here("Temp/2period-sim-results.rds"))
#----------------------------------------------------------------------

