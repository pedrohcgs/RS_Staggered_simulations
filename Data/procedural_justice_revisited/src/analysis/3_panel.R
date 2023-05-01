##########################################################################
# Description: Generate outcomes at officer-level per month
#              Create unbalanced and balanced panel, where the latter
#              excludes officers who resigned or retired during the study
##########################################################################

# generate training frame
training <- 
  mutate(training,
         assigned_exact = assigned,
         assigned = round_date(assigned, "month"),
         resigned = floor_date(resigned, "month")) %>%
  remove_new_officers() %>%
  filter(assigned < resigned | is.na(resigned))

# outcome at officer-level, per month
y_complaints <-
  left_join(training, complaints, by = "uid") %>%
  group_by(uid, month = floor_date(date, "month")) %>%
  summarize(complaints = n_distinct(cr_id, na.rm = TRUE))

y_sustained <-
  left_join(training, filter(complaints, sustained == 1), # sustained_settled == 1
            by = "uid") %>%
  group_by(uid, month = floor_date(date, "month")) %>%
  summarize(sustained = n_distinct(cr_id, na.rm = TRUE))

y_force <-
  left_join(training, force, by = "uid") %>%
  group_by(uid, month = floor_date(date, "month")) %>%
  summarize(force = n_distinct(trr_id, na.rm = TRUE))

# join outcomes
pj <-
  expand_grid(
    uid   = unique(training$uid),
    month = seq(ymd("2011-01-01"), ymd("2016-12-01"), "1 month")
  ) %>%
  left_join(training,     by = "uid") %>%
  left_join(y_complaints, by = c("uid", "month")) %>%
  left_join(y_sustained,  by = c("uid", "month")) %>%
  left_join(y_force,      by = c("uid", "month"))

# first, replace NA entries (i.e. no complaints) with zeroes
pj <-
  replace_na(pj, list(complaints = 0, sustained = 0, force = 0))

# second, ensure that outcomes are NA in months after officers resign
pj <-
  mutate_at(
    pj,
    vars(complaints, sustained, force),
    ~ ifelse(month >= resigned & !is.na(resigned), NA, .)
  )

# set first trained period (first_trained = 0 for untrained officers)
pj <-
  mutate(pj, period = frank(month, ties.method = "dense")) %>%
  group_by(uid) %>%
  mutate(
    first_trained = min(period[month >= assigned]),
    first_trained = ifelse(is.infinite(first_trained), 0, first_trained)
  ) %>%
  ungroup()

# balanced panel, excluding all resigning officers
pj_balanced <-
  filter(pj, resigned > max(month) | is.na(resigned))

nrow(pj_balanced) ==
  length(unique(pj_balanced$uid)) * length(unique(pj_balanced$month))

# save data
save.image("output/rdata/pj.RData")

