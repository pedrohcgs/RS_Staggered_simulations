load("output/rdata/did.RData")

# outcome labels
outcome_labels <-
  as_labeller(c("complaints" = "Complaints",
                "sustained"  = "Sustained & Settled",
                "force"      = "Use of Force"))

# average in 12 months before training
did_truncated$mean <-
  pj_balanced %>%
  mutate(relative = period - first_trained) %>%
  filter(between(relative, -12, -1), month <= ymd("2016-03-01")) %>%
  summarize_at(vars(complaints, sustained, force), mean) %>%
  pivot_longer(cols = everything()) %>%
  pluck("value")

# tidy average effects
average_att <-
  transmute(
    did_truncated,
    outcome,
    weighted_att = simple$overall.att,
    weighted_se  = simple$overall.se,
    dynamic_att  = dynamic$overall.att,
    dynamic_se   = dynamic$overall.se,
    calendar_att = calendar$overall.att,
    calendar_se  = calendar$overall.se,
    cohort_att   = selective$overall.att,
    cohort_se    = selective$overall.se
  ) %>%
  pivot_longer(cols = weighted_att:cohort_se) %>%
  separate(name, sep = "_", into = c("aggregation", "estimand")) %>%
  pivot_wider(names_from = estimand, values_from = value) %>%
  left_join(transmute(did_truncated, outcome, mean = round(mean, 3))) %>%
  mutate(
    outcome     = fct_relevel(outcome, "complaints", "sustained"),
    aggregation = fct_relevel(aggregation, "weighted", "dynamic"),
    lower       = att - qnorm(0.975) * se,
    upper       = att + qnorm(0.975) * se
  )

average_text <- distinct(average_att, outcome, mean)

# plot average effects
average_att %>%
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(aes(aggregation, ymin = lower, ymax = upper),
                width = 0.05, color = grey(0.4)) +
  geom_label(aes(aggregation, att, label = sprintf("%0.3f", round(att, 3))),
             size = 3) +
  geom_text(
    data = average_text, hjust = 1.1, vjust = 1.5, parse = TRUE, size = 3,
    aes(x = Inf, y = Inf,
        label = paste0("bar(y)['[-12,-1]']~'='~", mean))
  ) +
  scale_x_discrete("Aggregation", labels = function(x) str_to_title(x)) +
  scale_y_continuous("Average of ATTs", breaks = scales::pretty_breaks()) +
  facet_wrap(~ outcome, labeller = outcome_labels, scales = "free_y") +
  ggsave("output/figures/truncated_average_effects.pdf",
         height = 4, width = 9)

# plot effects
truncated_figures <-
  transmute(
    did_truncated,
    dyn      = list(
      plot_dyn(dynamic, dynamic_se,
               ylab = word(outcome, 2, sep = "_"), mean = mean)
    )
  )

# dynamic effects
wrap_plots(truncated_figures$dyn, ncol = 3) +
  ggsave("output/figures/truncated_dynamic_effects.pdf",
         height = 4, width = 10)
