##########################################################################
# Description: Plot officers per cluster and mean outcomes by month
##########################################################################

load("output/rdata/pj.RData")

# outcome labels
outcome_labels <-
  as_labeller(c("complaints" = "Complaints",
                "sustained"  = "Sustained & Settled",
                "force"      = "Use of Force"))

# plot officers per cluster
training %>%
  group_by(assigned_exact) %>%
  summarize(officers = n_distinct(uid)) %>%
  ggplot(aes(assigned_exact, officers)) +
  geom_point(shape = 21) +
  scale_x_date("Assigned to training") +
  scale_y_continuous("Officers in cluster") +
  ggsave(filename = "output/figures/cluster_size.pdf",
         width = 6, height = 5)

# plot mean outcomes by month
pj_balanced %>%
  group_by(month) %>%
  summarize_at(vars(complaints, sustained, force), mean) %>%
  mutate(sustained = ifelse(month >= settlement_terminus, NA, sustained)) %>%
  pivot_longer(cols = c(complaints:force),
               names_to = "outcome", values_to = "mean") %>%
  mutate(data = ifelse(month >= ymd("2016-03-01"), "Updated", "Original"),
         outcome = fct_relevel(outcome, "complaints", "sustained")) %>%
  ggplot(aes(month, mean, color = data)) +
  geom_point() +
  scale_x_date(NULL) +
  scale_y_continuous("Mean per officer",
                     breaks = scales::pretty_breaks()) +
  scale_color_brewer("Data", palette = "Set2", direction = -1) +
  facet_wrap(~ outcome, labeller = outcome_labels) +
  theme(legend.position = "bottom") +
  ggsave(filename = "output/figures/mean_outcomes.pdf",
         width = 8, height = 5)

