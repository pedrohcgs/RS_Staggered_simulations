# import all complaints data
case <-
  tibble(filename = list.files("data/complaints",
                               pattern = "complaints|case",
                               full.names = TRUE)) %>%
  transmute(
    data = map(
      filename, 
      function(x) 
        read_data(x) %>% 
        clean_names() %>%
        select(one_of("log_no", "cr_id", "incident_fromdate",
                      "incident_date", "incident_start_date")) %>%
        set_names(c("cr_id", "date")) %>%
        mutate(
          date = as_date(parse_date_time(date, orders = c("dmy", "mdy", "ymd")))
        )
    )
  ) %>%
  unnest() %>%
  distinct()

accused <-
  tibble(filename = list.files("data/complaints",
                               pattern = "accused",
                               full.names = TRUE)) %>%
  mutate(data = map(
    filename, function(x)
      x %>%
      read_data() %>%
      clean_names() %>%
      select(one_of("log_no", "cr_id", "uid", "officer_last_name", "last_name",
                    "officer_first_name", "first_name", "middle_initial",
                    "birth_year", "final_finding")) %>%
      rename(cr_id = 1) %>%
      left_join(case, by = "cr_id")
  ))

allegation <-
  left_join(
    clean_names(read_excel("data/complaints/allegation.xlsx", sheet = 2)),
    clean_names(read_excel("data/complaints/allegation.xlsx", sheet = 5)),
    by = "officer_id"
  ) %>%
  transmute(
    cr_id = as.numeric(crid), date = as_date(ymd_hms(incident_date)),
    officer_id, last_name = toupper(officer_last.y),
    first_name = toupper(officer_first),
    appointed = as_date(appt_date, origin = "1899-12-30")
  )

# tidy pre March 2016 data
complaints_1 <-
  bind_rows(accused$data[[1]],
            accused$data[[2]]) %>%
  group_by(cr_id, uid, date) %>%
  transmute(finding = ifelse(any(final_finding %in% c("SUSTAINED", "SU")),
                             "sustained", "other_outcome")) %>%
  ungroup()

# tidy post March 2016 data and join uid
complaints_2 <-
  bind_rows(accused$data[[3]],
            accused$data[[4]]) %>%
  transmute(
    cr_id, last_name = remove_suffix(coalesce(officer_last_name, last_name)),
    first_name = coalesce(officer_first_name, first_name),
    mi = middle_initial, birth_year
  ) %>%
  distinct() %>%
  left_join(case, by = "cr_id") %>%
  left_join(transmute(officers, last_name, first_name, mi, birth_year, uid),
            by = c("last_name", "first_name"))  %>%
  resolve_join()

# tidy supplemental allegation data and join uid
complaints_3 <-
  allegation %>%
  left_join(transmute(officers, last_name, first_name, appointed, uid),
            by = c("last_name", "first_name", "appointed")) %>%
  transmute(uid, cr_id, date)

# read settlements data
settlements <-
  read_csv("data/settlements/settlements_1952_2016.csv") %>%
  transmute(uid = UID, date = incident_date,
            settlement = parse_number(settlement))

# bind complaints and join settlements data 
complaints <-
  bind_rows(
    complaints_1,
    complaints_2,
    complaints_3
  ) %>%
  distinct() %>%
  left_join(settlements, by = c("uid", "date")) %>%
  mutate(sustained = ifelse(finding == "sustained" | !is.na(settlement), 1, 0))
