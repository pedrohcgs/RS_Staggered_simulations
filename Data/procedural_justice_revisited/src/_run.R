source("src/_packages.R")
source("src/_functions.R")

# import and wrangle data
source("src/analysis/1_import.R")
source("src/analysis/2_training.R")
source("src/analysis/3_panel.R")

# run analysis
source("src/analysis/4_did.R")

# plots and output for memo
source("src/plots/plot_did.R")
source("src/plots/plot_truncated.R")
source("src/plots/plot_other.R")

# plots and output for correction
source("src/plots/plot_fig2.R")
