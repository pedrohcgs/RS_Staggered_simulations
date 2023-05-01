##########################################################################
# File: 3_estimation.R
# Description: Estimate ATT using IFE model
#              Estimate ATT on early- and late-adopters
##########################################################################

source(here::here("Data/procedural_justice-master/src/aux_files/packages.R"))
source(here("Data/procedural_justice-master/src/aux_files/functions.R"))
load(here("Data/procedural_justice-master/products/rdata/2_join.RData"))

# set constants
nboots <- 2000

# build training matrix
y <-
  generate_tscs(
    training,
    complaints = complaints,
    sustained  = filter(complaints, !is.na(settlement) | sustained == 1),
    force      = force,
    cluster_by = "day"
  ) %>%
  mutate(coverage = ifelse(month <= settlement_terminus, "sfc", "fc")) %>%
  filter(month <= discontinuity)

# view design
panelView(complaints ~ trained, data = y,
          index = c("cluster", "period"), pre.post = TRUE, by.timing = TRUE)

# estimate ATT
models <-
  expand_grid(
    formula = c(complaints ~ trained,
                sustained  ~ trained,
                force      ~ trained)
  ) %>%
  mutate(
    outcome = map_chr(formula, ~ as_string(lhs(.))),
    fit     = map(formula,
                  ~ estimate_model(y, formula = .x, method = "ife",
                                            nboots = nboots, seed = 29590))
  )

# estimate possibly heterogenous effects for early- and late-adopters
early_adopters <- pull_adopters(y, ymd("2012-05-01"), ymd("2013-09-01"))
late_adopters  <- pull_adopters(y, ymd("2014-04-01"), ymd("2016-03-01"))

models <-
  mutate(
    models,
    early = map(formula,
                ~ estimate_model(
                  filter(y, !cluster %in% late_adopters, cluster != 1),
                  formula = .x, method = "ife",
                  nboots = nboots, seed = 29590)),
    late  = map(formula,
                ~ estimate_model(
                  filter(y, !cluster %in% early_adopters, cluster != 1),
                   formula = .x, method = "ife",
                  nboots = nboots, seed = 29590))
  )

# save
save.image(here("Data/procedural_justice-master/products/rdata/3_estimation.RData"))


