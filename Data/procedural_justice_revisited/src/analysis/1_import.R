##########################################################################
# Description: Import officer, complaints, settlements, use of force,
#              and training data
##########################################################################

# read and tidy officer data
officers <-
  read_csv("data/officers/officers.csv") %>%
  transmute(
    uid = UID, last_name = remove_suffix(last_name),
    first_name, mi = middle_initial, appointed = appointed_date,
    resigned = resignation_date, birth_year
  ) %>%
  filter(!is.na(appointed)) %>%
  distinct(
    uid, last_name, first_name, birth_year, appointed, resigned,
    .keep_all = TRUE
  )

# import outcome data
source("src/analysis/_complaints.R")
source("src/analysis/_force.R")

# read training data and filter to first training assignment for each officer
training <-
  read_csv("data/training/training.csv") %>%
  group_by(scrambled_id, last_name, first_name) %>%
  filter(assigned == min(assigned)) %>%
  ungroup()

# the data in the original study ended at 2016-03
original_end <- ymd("2016-03-01")

complaints %>%
  filter(date >= "2011-01-01", date <= ymd("2016-12-31")) %>%
  group_by(month = floor_date(date, "month")) %>%
  summarize(n = n_distinct(cr_id)) %>%
  ggplot(aes(month, n)) +
  geom_vline(xintercept = original_end, linetype = "dashed") +
  geom_line() +
  geom_point() +
  scale_y_continuous("Complaints", limits = c(0, 400))

force %>%
  filter(date >= "2011-01-01") %>%
  group_by(month = floor_date(date, "month")) %>%
  summarize(n = n_distinct(trr_id)) %>%
  ggplot(aes(month, n)) +
  geom_vline(xintercept = original_end, linetype = "dashed") +
  geom_line() +
  geom_point() +
  scale_y_continuous("Force", limits = c(0, 700))

# truncate data at 2016-12-31, all officers are trained by this date
complaints <-
  filter(
    complaints,
    between(date,
            min(training$assigned) %>% floor_date("month") %m-% months(12),
            ymd("2016-12-31"))
  )

force <-
  filter(
    force,
    between(date,
            min(training$assigned) %>% floor_date("month") %m-% months(12),
            ymd("2016-12-31")) 
  )

# settlement data ends at 2015-11
settlement_terminus <-
  max(complaints$date[!is.na(complaints$settlement)], na.rm = TRUE) %>%
  floor_date("month") %m-% days(1)

