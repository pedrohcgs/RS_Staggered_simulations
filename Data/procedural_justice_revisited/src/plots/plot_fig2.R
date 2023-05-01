load("output/rdata/did.RData")

# plot theme
theme_set(
  theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank())
)

# functions
tidy_dynamic <-
  function(dynamic, dynamic_se, outcome, time = c(-36, 24)) {
    tibble(
      t       = dynamic$egt[between(dynamic$egt, time[1], time[2])],
      att     = dynamic$att.egt[between(dynamic$egt, time[1], time[2])],
      se      = dynamic_se$se[1:length(time[1]:time[2])],
      crit    = dynamic_se$crit.val,
      outcome = fct_relevel(str_to_title(outcome), "Complaints", "Sustained"),
    ) %>%
      mutate(
        lower      = att - qnorm(0.975) * se,
        upper      = att + qnorm(0.975) * se,
        lower_crit = att - crit * se,
        upper_crit = att + crit * se
      ) %>%
      mutate_at(vars(att, lower, upper, lower_crit, upper_crit), ~ . * 100) %>%
      group_by(outcome) %>%
      mutate(fill_att = (att - mean(att)) / sd(att)) %>%
      ungroup()
  }

tidy_breaks <-
  function(x) {
    l <- max(abs(x))
    if (l >= 1) {
      l <- round_half_up(l, digits = 0)
      seq(-l, l, by = 1)
    } else {
      round(seq(-0.6, 0.6, by = 0.2), digits = 1)
    }
  }

# extract data for Figure 2
figure_2_data <-
  mutate(did,
         for_plot = list(tidy_dynamic(dynamic, dynamic_se, outcome))) %>%
  pluck("for_plot") %>%
  bind_rows()

# plot Figure 2A under new estimator
figure_2_top <-
  left_join(
    figure_2_data,
    pj_balanced %>%
      group_by(t = period - first_trained) %>%
      summarize(
        complaints = sum(complaints) / n_distinct(uid) * 100,
        sustained  = sum(sustained) / n_distinct(uid) * 100,
        force      = sum(force) / n_distinct(uid) * 100,
      ) %>%
      pivot_longer(cols = c(complaints:force),
                   names_to = "outcome",
                   values_to = "observed") %>%
      mutate(outcome = fct_relevel(str_to_title(outcome),
                                   "Complaints", "Sustained"))
  ) %>%
  transmute(outcome, t, att, observed, counterfactual = observed - att) %>%
  pivot_longer(cols = c(observed, counterfactual)) %>%
  ggplot(aes(t, value, color = str_to_title(name),
             linetype = str_to_title(name))) +
  geom_vline(xintercept = 0, color = "grey40", linetype = "dashed") +
  geom_line() +
  geom_point(alpha = 0.2) +
  scale_x_continuous(NULL,
                     expand = expansion(mult = 0.05),
                     breaks = seq(-48, 36, by = 12)) +
  scale_y_continuous("Frequency per 100 Officers",
                     limits = c(0, NA),
                     expand = expansion(mult = c(0.02, 0.07)),
                     breaks = scales::pretty_breaks()) +
  scale_color_manual(NULL, values = c("grey20", "#3498db")) +
  scale_linetype_manual(NULL, values = c("longdash", "solid")) +
  facet_wrap(~ outcome, scales = "free_y") +
  theme(axis.line        = element_line(size = 0.3),
        axis.text        = element_text(color = "black"),
        axis.ticks       = element_line(size = 0.3), 
        axis.title.y     = element_text(margin = margin(r = 0.8, unit = "lines")),
        strip.text       = element_text(size = 9),
        legend.direction = "vertical",
        legend.justification = c(0, 0), 
        legend.position = c(0.005, 0.05))

# plot Figure 2B under new estimator
figure_2_bottom <-
  ggplot(figure_2_data, aes(t, att, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, color = "grey40", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey40", linetype = "dashed") +
  geom_ribbon(aes(ymin = lower_crit, ymax = upper_crit), alpha = 0.15) +
  geom_ribbon(fill = "#1f77b4", alpha = 0.15) +
  geom_line(color = "grey60", size = 0.4) +
  geom_point(aes(fill = fill_att, group = outcome),
             shape = 21, stroke = 0.4, color = "grey10", show.legend = FALSE) +
  scale_x_continuous(
    "Month Relative to Training",
    expand = expansion(mult = 0.05),
    breaks = seq(-48, 36, by = 12)
  ) +
  scale_y_continuous(
    "ATT on Frequency",
    limits = function(x) c(-max(abs(x)), max(abs(x))),
    expand = expansion(mult = 0.05),
    breaks = function(x) tidy_breaks(x)
  ) +
  colorspace::scale_fill_continuous_diverging(palette = "Purple-Brown",
                                              rev = FALSE, mid = 0.63) +
  facet_wrap(~ outcome, scales = "free_y") +
  theme(panel.grid = element_blank(),
        axis.line  = element_line(size = 0.3),
        axis.text  = element_text(color = "black"),
        axis.ticks = element_line(size = 0.3),
        strip.text = element_blank())

# compose figure and write
figure_2 <-
  ggarrange(figure_2_top, figure_2_bottom,
            ncol = 1, common.legend = FALSE)

ggsave(figure_2, filename = "output/figures/figure_2.pdf",
       width = 10, height = 5, device = cairo_pdf)


