# clean suffixes
remove_suffix <- 
  function(x) {
    suffix <- c("SR", "JR", "II", "III", "IV")
    split  <- stri_split(x, regex = "\\s|\\.|\\,|\\(|\\)")
    x <-
      sapply(split,
             function(z)
               paste(z[!(z %in% suffix)], collapse = " "),
             USE.NAMES = FALSE)
    x <- trimws(x)
    return(x)
  }

# read csv or excel files
read_data <-
  function(x) {
    if (tools::file_ext(x) == "csv") {
      read_csv(x)
    } else {
      read_excel(x)
    }
  }

# filter new officers
remove_new_officers <-
  function(x, end_date) {
    x %>%
      filter(min(assigned) %m-% months(12) >= appointed)
  }

# identify uid for officers receiving a complaint using names, then resolve
# duplicates using birth year followed by middle initial
resolve_join <-
  function(x) {
    r1 <-
      group_by(x, cr_id, last_name, first_name, date, sustained) %>%
      filter(n() == 1)
    
    r2 <-
      get_dupes(x, cr_id, date, last_name, first_name) %>%
      filter(birth_year.x == birth_year.y) %>%
      group_by(cr_id, date, last_name, first_name, sustained) %>%
      filter(n() == 1)
    
    r3 <-
      get_dupes(x, cr_id, date, last_name, first_name) %>%
      filter(mi.x == mi.y) %>%
      group_by(cr_id, date, last_name, first_name, sustained) %>%
      filter(n() == 1)
    
    r4 <-
      get_dupes(x, cr_id, date, last_name, first_name) %>%
      filter(year(date) - birth_year.x <= 65,
             year(date) - birth_year.y <= 65) %>%
      group_by(cr_id, date, last_name, first_name, sustained) %>%
      filter(n() == 1)
    
    bind_rows(r1, r2, r3, r4) %>%
      ungroup() %>%
      distinct(uid, cr_id, date, sustained) %>%
      na.omit()
  }

# bootstrap inference for dynamic effects
bootstrap_dynamic <-
  function(fit, dynamic, time = c(-36, 36)) {
    e_select        <- between(dynamic$egt, time[1], time[2])
    egt             <- dynamic$egt[e_select]
    ES_att          <- dynamic$att.egt[e_select]
    ES_inf_function <- dynamic$inf.function$dynamic.inf.func.e[, e_select]
    boot_es         <- did::mboot(ES_inf_function, fit$DIDparams)
    return(boot_es)
  }

# plot dynamic effects
plot_dyn <-
  function(dynamic, dynamic_se, time = c(-36, 36), ylab = NULL,
           mean) {
    tibble(
      t    = dynamic$egt[between(dynamic$egt, time[1], time[2])],
      att  = dynamic$att.egt[between(dynamic$egt, time[1], time[2])],
      se   = dynamic_se$se,
      crit = dynamic_se$crit.val
    ) %>%
      mutate(
        lower      = att - qnorm(0.975) * se,
        upper      = att + qnorm(0.975) * se,
        lower_crit = att - crit * se,
        upper_crit = att + crit * se
      ) %>%
      ggplot(aes(t, att)) +
      geom_hline(yintercept = 0, color = grey(0.4), linetype = "dashed") +
      geom_vline(xintercept = 0, color = grey(0.4), linetype = "dashed") +
      geom_ribbon(aes(ymin = lower_crit, ymax = upper_crit),
                  fill = grey(0.2), alpha = 0.15) +
      geom_ribbon(aes(ymin = lower, ymax = upper),
                  fill = "blue", alpha = 0.15) +
      geom_line(size = 0.2, color = grey(0.2)) +
      geom_point(size = 0.8, color = grey(0.2)) +
      annotate(
        "text", x = Inf, y = Inf, parse = TRUE,
        hjust = 1.1, vjust = 1.2, size = 3,
        label = as.expression(bquote(bar(y)["[-12,-1]"]~"="~.(round(mean, 3))))
      ) +
      scale_x_continuous("Month relative to training",
                         breaks = seq(-36, 36, by = 12)) +
      scale_y_continuous(glue("Dynamic ATT: ",
                              str_to_title(ylab), " per officer")) +
      theme(panel.grid.major = element_blank())
  }

# plot calendar and selective effects
plot_cs <-
  function(calendar, time = c(-Inf, Inf), ylab = NULL, mean) {
    
    tibble(
      t    = calendar$egt,
      att  = calendar$att.egt,
      se   = calendar$se.egt,
      crit = calendar$crit.val.egt
    ) %>%
      mutate(
        lower      = att - qnorm(0.975) * se,
        upper      = att + qnorm(0.975) * se,
        lower_crit = att - crit * se,
        upper_crit = att + crit * se
      ) %>%
      filter(t >= time[1], t <= time[2]) %>%
      ggplot(aes(t, att)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_errorbar(aes(ymin = lower_crit, ymax = upper_crit),
                    width = 0, color = grey(0.6), size = 0.3) +
      geom_errorbar(aes(ymin = lower, ymax = upper),
                    width = 0, color = "blue", size = 0.3) +
      geom_point(size = 0.9, color = grey(0.2)) +
      annotate(
        "text", x = Inf, y = Inf, parse = TRUE,
        hjust = 1.1, vjust = 1.2, size = 3,
        label = as.expression(bquote(bar(y)["[-12,-1]"]~"="~.(round(mean, 3))))
      ) +
      scale_x_continuous("Month of training roll-out",
                         expand = expansion(add = c(3, 3)),
                         breaks = scales::pretty_breaks()) +
      scale_y_continuous(ylab,
                         breaks = scales::pretty_breaks())
  }

# plot theme
theme_set(
  theme_light(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = 0.2, linetype = "dashed"))
)



### LEGACY
# identify uid for officers receiving a complaint using names, then resolve
# duplicates using birth year followed by middle initial
resolve_join <-
  function(x) {
    r1 <-
      group_by(x, cr_id, last_name, first_name, date) %>%
      filter(n() == 1)
    
    r2 <-
      get_dupes(x, cr_id, date, last_name, first_name) %>%
      filter(birth_year.x == birth_year.y) %>%
      group_by(cr_id, date, last_name, first_name) %>%
      filter(n() == 1)
    
    r3 <-
      get_dupes(x, cr_id, date, last_name, first_name) %>%
      filter(mi.x == mi.y) %>%
      group_by(cr_id, date, last_name, first_name) %>%
      filter(n() == 1)
    
    r4 <-
      get_dupes(x, cr_id, date, last_name, first_name) %>%
      filter(year(date) - birth_year.x <= 65,
             year(date) - birth_year.y <= 65) %>%
      group_by(cr_id, date, last_name, first_name) %>%
      filter(n() == 1)
    
    bind_rows(r1, r2, r3, r4) %>%
      ungroup() %>%
      distinct(uid, cr_id, date) %>%
      na.omit()
  }
