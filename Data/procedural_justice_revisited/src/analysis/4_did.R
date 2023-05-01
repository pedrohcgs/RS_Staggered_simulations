##########################################################################
# Description: Estimate DID model using Callaway and Sant'Anna (2020)
#              procedure
#              Use balanced panel (excluding officers who retired during study)
##########################################################################

# model for balanced panel
balanced_model <-
  function(yname = "complaints",
           data = pj_balanced,
           truncated = FALSE) {
    
    if (yname == "sustained") {
      data <- filter(data, month <= settlement_terminus)
    }
    
    if (truncated) {
      data <- filter(data, month < ymd("2016-03-01"))
    }
    
    att_gt(
      yname            = yname,
      tname            = "period", 
      idname           = "uid",
      first.treat.name = "first_trained",
      data             = data,
      bstrap           = FALSE,
      panel            = TRUE,
      control.group    = "notyettreated"
    )
    
  }

# estimate models and aggregate effects
did <-
  tibble(outcome = c("complaints", "sustained", "force")) %>%
  rowwise() %>%
  mutate(
    fit        = list(balanced_model(yname = outcome)),
    simple     = list(aggte(fit, type = "simple")),
    dynamic    = list(aggte(fit, type = "dynamic")),
    calendar   = list(aggte(fit, type = "calendar")),
    selective  = list(aggte(fit, type = "selective")),
    dynamic_se = list(bootstrap_dynamic(fit, dynamic))
  )

did_truncated <-
  tibble(outcome = c("complaints", "sustained",  "force")) %>%
  rowwise() %>%
  mutate(
    fit        = list(balanced_model(yname = outcome, truncated = TRUE)),
    simple     = list(aggte(fit, type = "simple")),
    dynamic    = list(aggte(fit, type = "dynamic")),
    calendar   = list(aggte(fit, type = "calendar")),
    selective  = list(aggte(fit, type = "selective")),
    dynamic_se = list(bootstrap_dynamic(fit, dynamic))
  )

# save
save.image("output/rdata/did.RData")

