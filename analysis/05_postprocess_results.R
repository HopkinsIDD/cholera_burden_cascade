# This script computes aggregates of inferred incidence and sero-exposures


# Preamble ----------------------------------------------------------------

library(tidyverse)
library(posterior)
library(cmdstanr)
library(here)

source("analysis/utils.R")

# Lazy parameters
n_cores <- 2    # number of cores for processing
chain_subset <- NULL      # chains to subset for computations, NULL for all
draw_subset <- NULL    # subset of draws, NULL for all
do_plots <- F    # make intermediate plots
t_min_plot <- as.Date("2021-01-24") # tmin for plotting
t_max_plot <- getTmax()   # tmin for plotting
datelim <- c(t_min_plot, t_max_plot)
date_subset <- seq.Date(t_min_plot, t_max_plot, by = "1 days")
diarrhea_severity <- "diarrhea2"

save(
  list = c("chain_subset", "draw_subset"),
  file = here(str_glue("generated_data/options_bundle_for_figures_{diarrhea_severity}.rdata"))
)

# Load data ---------------------------------------------------------------

# Serosurvey round dates
dates_wide <- readRDS(here("data/dates_full_wide_agegrpall.rds"))
study_period_dates <- seq.Date(min(dates_wide$R0), max(dates_wide$R2), by = "1 days")

# Full clinical data
case_data_subset_full <- readRDS(here("generated_data/case_data_subset_full.rds"))
case_data_subset_filtered <- readRDS(here("generated_data/case_data_subset_filtered.rds"))

# Clinical data
lab_data <- readRDS(here("data/lab_data.rds"))

# Aggregated weeklyclinical data
period_data <- readRDS(here("data/period_data.rds"))

# Age categories in rest of computations
age_categories <- c(unique(period_data$age_cat), "overall")

# Add rectangles for rounds
round_dates <- dates_wide %>% 
  distinct(R0, R1, R2) %>% 
  pivot_longer(cols = everything(),
               names_to = "round",
               values_to = "date") %>% 
  mutate(round = case_when(date < "2021-05-01" & round == "R0" ~ "R0A",
                           date > "2021-05-01" & round == "R0" ~ "R0B",
                           T ~ round)) %>% 
  group_by(round) %>% 
  summarise(tl = min(date),
            tr = max(date)) %>% 
  rowwise() %>% 
  mutate(tmid = mean.Date(c(tl, tr)))

# Plot for totals
lab_data <- lab_data %>% 
  filter(date %in% case_data_subset_filtered$date) %>% 
  addEpiWeek() %>% 
  mutate(obs_id = row_number()) %>% 
  group_by(epiweek) %>% 
  arrange(epiweek, desc(rdt_res), desc(pcr_res), desc(cult_res)) %>% 
  mutate(obs_id_week = row_number(),
         tl = min(date),
         tr = max(date)) %>% 
  rowwise() %>% 
  mutate(tmid = mean.Date(c(tl, tr))) %>% 
  ungroup()


# Get population data in 2021
pop_dat <- getSitakundatData(age_categories = age_categories) %>% 
  { 
    x <- .
    bind_rows(x, 
              x %>% summarise(pop = sum(pop)) %>% mutate(age_cat= "overall"))
  }


# Titer data
vc_titer_data <- readRDS(here("data/vc_ogawa_titers_full_agegrpall.rds"))

# Dictionnary for age categories
age_cat_dict <- getAgeCatDict()

# Population by age-category for post-stratification
sample_pop_age <- vc_titer_data %>% 
  distinct(id, age_grp) %>%
  group_by(age_grp) %>% 
  summarise(n = n()) %>% 
  mutate(age_cat = names(age_cat_dict)[age_grp]) %>% 
  ungroup()

# Periods for outputs
output_periods_overall <- getOutputPeriodOverall()
output_periods <- getWeeklyOutputPeriods()

# Subset of dates to compute the output statistics
ouptut_dates_subset <- seq.Date(output_periods_overall$tl, output_periods_overall$tr, by = "1 days")


# Save data for figures
save(
  list = c("lab_data", "round_dates", "vc_titer_data", "sample_pop_age",
           "pop_dat", "period_data", "age_cat_dict", "age_categories",
           "output_periods_overall", "output_periods", "ouptut_dates_subset"),
  file = here("generated_data/data_bundle_for_figures.rdata")
)

# Load estimates ----------------------------------------------------------

# Clinical incidence estimates
res_incid <- readRDS(here("generated_data/incid_prior_survonly_stan_output_updated_multiage.rds"))
# prior_incid <- readRDS(here("generated_data/incid_prior_survonly_stan_output_updated_multiage_prior.rds"))

# Seroincidence model estimates
res_sero <- readRDS(here("generated_data/stan_output_final_envlambda_agegrpall_J100_ogawa.rds"))
# prior_sero <- readRDS(here("generated_data/stan_output_jointfullpol_sitakunda_agegrpall_J100_ogawa_prior.rds"))

# Exposure classes based on serology
prob_case <-  readRDS(here("generated_data/stan_output_stats_final_envlambda_agegrpall_J100_ogawa.rds"))$prob_case

# Estimates of exposure using serology
genquant <- readRDS("generated_data/genquant_output_final_envlambda_agegrpall_J100.rds")
genquant_cumulative <- readRDS("generated_data/genquant_output_final_envlambda_agegrpall_J100_cumulative.rds")
genquant_overall <- readRDS("generated_data/genquant_output_final_envlambda_agegrpall_J100_overall.rds")


# A. Clinical incidence aggregates -------------------------------------------

# Draws of daily cholera incidence
I_draws <- res_incid$draws("I_chol") %>%
  as_draws() %>%
  as_draws_df() %>%
  as_tibble() %>%
  pivot_longer(contains("I_chol")) %>%
  mutate(i = str_extract(name, "(?<=\\[)[0-9]+") %>% as.numeric(),
         age_cat = str_extract(name, "(?<=\\,)[0-9]+") %>% as.numeric()) %>%
  mutate(date = case_data_subset_full$date[i])


tot_I_sum_draws <- aggregateIncid(I_draws,
                                  chain_subset = chain_subset,
                                  draw_subset = draw_subset,
                                  date_subset = date_subset,
                                  by_age = FALSE,
                                  weekly = FALSE,
                                  cumulative = FALSE,
                                  do_stats = F,
                                  age_categories = age_categories)

tot_I_sum_age_draws <- aggregateIncid(I_draws,
                                      chain_subset = chain_subset,
                                      draw_subset = draw_subset,
                                      date_subset = date_subset,
                                      by_age = TRUE,
                                      weekly = FALSE,
                                      cumulative = FALSE,
                                      do_stats = F,
                                      age_categories = age_categories)

tot_I_sum <- aggregateIncid(I_draws,
                            chain_subset = chain_subset,
                            draw_subset = draw_subset,
                            date_subset = date_subset,
                            by_age = FALSE,
                            weekly = FALSE,
                            cumulative = FALSE)

tot_I_sum_age <- aggregateIncid(I_draws,
                                chain_subset = chain_subset,
                                draw_subset = draw_subset,
                                date_subset = date_subset,
                                by_age = TRUE,
                                weekly = FALSE,
                                cumulative = FALSE)

weekly_I_sum <- aggregateIncid(I_draws,
                               chain_subset = chain_subset,
                               draw_subset = draw_subset,
                               date_subset = date_subset,
                               by_age = FALSE,
                               weekly = TRUE,
                               cumulative = FALSE)

# Save for use in sketch
saveRDS(weekly_I_sum, here("generated_data/weekly_I_sum.rds"))


weekly_I_cumsum <- aggregateIncid(I_draws,
                                  chain_subset = chain_subset,
                                  draw_subset = draw_subset,
                                  date_subset = date_subset,
                                  by_age = FALSE,
                                  weekly = TRUE,
                                  cumulative = TRUE)

weekly_I_sum_age <- aggregateIncid(I_draws,
                                   chain_subset = chain_subset,
                                   draw_subset = draw_subset,
                                   date_subset = date_subset,
                                   weekly = TRUE,
                                   cumulative = FALSE,
                                   by_age = TRUE,
                                   age_categories = age_categories)

weekly_I_cumsum_age <- aggregateIncid(I_draws,
                                      chain_subset = chain_subset,
                                      draw_subset = draw_subset,
                                      date_subset = date_subset,
                                      by_age = TRUE,
                                      weekly = TRUE,
                                      cumulative = TRUE,
                                      age_categories = age_categories)

# Save for figures
save(
  list = ls() %>% str_subset("weekly|sum") %>% str_subset("draws", negate = T),
  file = here("generated_data/clinical_incidence_bundle_for_figures.rdata")
)

# B. Health seeking probability ----------------------------------------------

# Account for health seeking
set.seed(235908)
# main analysis on "diarrhea2" severity
healthprob_data <- getHealthSeekingData(severity = diarrhea_severity)
healthprob_draws <- getHealthSeekingDraws(hs_data = healthprob_data)

# Save for further use
healthprob_stats <- healthprob_draws %>% 
  group_by(age_cat) %>% 
  summarise(q025 = quantile(p, 0.025),
            q975 = quantile(p, .975),
            mean = mean(p),
            median = median(p)) %>% 
  inner_join(healthprob_data %>% 
               rename(age_cat = age_group,
                      seek_health = x,
                      n_sample = n)) %>% 
  mutate(age_cat = age_cat_dict[age_cat])

saveRDS(healthprob_stats, file = here("generated_data/healthprob_stats.rds"))

save(
  list = c("healthprob_stats", "healthprob_draws"),
  file = here(str_glue("generated_data/healthprob_bundle_for_figures_{diarrhea_severity}.rdata")) 
)

# C. Community incidence aggregates ------------------------------------------

I_draws_true <- I_draws %>%
  inner_join(healthprob_draws) %>% 
  # Increase by reporting probability
  mutate(value = value/p)

tot_I_sum_true <- aggregateIncid(I_draws = I_draws_true,
                                 chain_subset = chain_subset,
                                 draw_subset = draw_subset,
                                 date_subset = date_subset,
                                 by_age = FALSE,
                                 weekly = FALSE,
                                 cumulative = FALSE)

tot_I_sum_true_draws <- aggregateIncid(I_draws_true,
                                       chain_subset = chain_subset,
                                       draw_subset = draw_subset,
                                       date_subset = date_subset,
                                       by_age = FALSE,
                                       weekly = FALSE,
                                       cumulative = FALSE,
                                       do_stats = F,
                                       age_categories = age_categories)

tot_I_sum_true_age_draws <- aggregateIncid(I_draws_true,
                                           chain_subset = chain_subset,
                                           draw_subset = draw_subset,
                                           date_subset = date_subset,
                                           by_age = TRUE,
                                           weekly = FALSE,
                                           cumulative = FALSE,
                                           do_stats = F,
                                           age_categories = age_categories)

tot_I_sum_age_true <- aggregateIncid(I_draws = I_draws_true,
                                     chain_subset = chain_subset,
                                     draw_subset = draw_subset,
                                     date_subset = date_subset,
                                     by_age = TRUE,
                                     weekly = FALSE,
                                     cumulative = FALSE)

weekly_I_sum_true <- aggregateIncid(I_draws = I_draws_true,
                                    chain_subset = chain_subset,
                                    draw_subset = draw_subset,
                                    date_subset = date_subset,
                                    weekly = TRUE,
                                    cumulative = FALSE,
                                    by_age = FALSE)

# Save for use in sketch
saveRDS(weekly_I_sum_true, here("generated_data/weekly_I_sum_true.rds"))

weekly_I_sum_age_true <- aggregateIncid(I_draws = I_draws_true,
                                        chain_subset = chain_subset,
                                        draw_subset = draw_subset,
                                        date_subset = date_subset,
                                        weekly = TRUE,
                                        cumulative = FALSE,
                                        by_age = TRUE,
                                        age_categories = age_categories)

weekly_I_cumsum_true <- aggregateIncid(I_draws = I_draws_true,
                                       chain_subset = chain_subset,
                                       draw_subset = draw_subset,
                                       date_subset = date_subset,
                                       weekly = TRUE,
                                       cumulative = TRUE,
                                       by_age = FALSE)

weekly_I_cumsum_age_true <- aggregateIncid(I_draws = I_draws_true,
                                           chain_subset = chain_subset,
                                           draw_subset = draw_subset,
                                           date_subset = date_subset,
                                           weekly = TRUE,
                                           cumulative = TRUE,
                                           by_age = TRUE,
                                           age_categories = age_categories)


# Save for figures
save(
  list = ls() %>% str_subset("true") %>% str_subset("draws", negate = T),
  file = here(str_glue("generated_data/community_incidence_bundle_for_figures_{diarrhea_severity}.rdata"))
)

# D. Seroinfection aggregates ------------------------------------------------

prob_output_overall_draws <- genquant_overall %>%
  tidybayes::spread_draws(mu_output_total[i,j]) %>%
  mutate(age_cat = j,
         age_cat = age_categories[as.numeric(age_cat)],
         period_id = i,
         tl = output_periods_overall$tl[period_id],
         tr = output_periods_overall$tr[period_id]) %>% 
  ungroup() %>% 
  select(-i, -j, -period_id)

prob_output_stats <- genquant$summary(variables = "mu_output_total",
                                      .cores = n_cores) %>% 
  mutate(age_cat = str_extract(variable, "(?<=\\,)[0-9]"),
         age_cat = age_cat_dict[age_categories[as.numeric(age_cat)]],
         period_id = str_extract(variable, "(?<=\\[)[0-9]+") %>% as.numeric(),
         tl = output_periods$tl[period_id],
         tr = output_periods$tr[period_id]) %>% 
  rowwise() %>% 
  mutate(tmid = mean.Date(c(tl, tr)))

# Post-stratified estimats
prob_output_overall_stats <- genquant %>%
  tidybayes::spread_draws(mu_output_total[i,j]) %>%
  mutate(age_cat = j,
         age_cat = age_categories[as.numeric(age_cat)],
         period_id = i) %>%
  inner_join(sample_pop_age) %>%
  group_by(.draw, period_id) %>% 
  summarise(mu_output_total = weighted.mean(mu_output_total, n)) %>% 
  group_by(period_id) %>% 
  summarise(mean = mean(mu_output_total),
            q5 = quantile(mu_output_total, .05),
            q95 = quantile(mu_output_total, .95)) %>% 
  mutate(tl = output_periods$tl[period_id],
         tr = output_periods$tr[period_id])  %>% 
  rowwise() %>% 
  mutate(tmid = mean.Date(c(tl, tr))) %>% 
  ungroup()


prob_output_cumulative_stats <- genquant_cumulative$summary(variables = "mu_output_total",
                                                            .cores = n_cores) %>% 
  mutate(age_cat = str_extract(variable, "(?<=\\,)[0-9]"),
         age_cat =  age_cat_dict[age_categories[as.numeric(age_cat)]],
         period_id = str_extract(variable, "(?<=\\[)[0-9]+") %>% as.numeric(),
         tl = output_periods$tl[period_id],
         tr = output_periods$tr[period_id]) %>% 
  rowwise() %>% 
  mutate(tmid = mean.Date(c(tl, tr))) %>% 
  ungroup()

prob_output_cumulative_overall_stats <- genquant_cumulative %>%
  tidybayes::spread_draws(mu_output_total[i,j]) %>%
  mutate(age_cat = j,
         age_cat = age_categories[as.numeric(age_cat)],
         period_id = i) %>%
  inner_join(sample_pop_age) %>%
  group_by(.draw, period_id) %>% 
  summarise(mu_output_total = weighted.mean(mu_output_total, n)) %>% 
  group_by(period_id) %>% 
  summarise(mean = mean(mu_output_total),
            q5 = quantile(mu_output_total, .05),
            q95 = quantile(mu_output_total, .95)) %>% 
  mutate(tl = output_periods$tl[period_id],
         tr = output_periods$tr[period_id]) %>% 
  rowwise() %>% 
  mutate(tmid = mean.Date(c(tl, tr))) %>% 
  ungroup()

# Save for future use in report
saveRDS(prob_output_cumulative_overall_stats, here("generated_data/prob_output_cumulative_overall_stats.rds"))


# Marginal effects of time-varying vs. constant FOI
marginal_ratio_draws <- genquant$draws("marginal_ratio") %>% 
  as_draws() %>% 
  as_draws_df() %>% 
  as_tibble() %>% 
  pivot_longer(cols = everything(),
               names_to = "variable",
               values_to = "draw") %>% 
  mutate(age_cat = str_extract(variable, "(?<=\\,)[0-9]"),
         period_id = str_extract(variable, "(?<=\\[)[0-9]+") %>% as.numeric(),
         tl = output_periods$tl[period_id],
         tr = output_periods$tr[period_id],
         age_cat = age_categories[as.numeric(age_cat)]) 

marginal_ratio_stats <- genquant$summary("marginal_ratio") %>% 
  mutate(age_cat = str_extract(variable, "(?<=\\,)[0-9]"),
         period_id = str_extract(variable, "(?<=\\[)[0-9]+") %>% as.numeric(),
         tl = output_periods$tl[period_id],
         tr = output_periods$tr[period_id],
         age_cat = age_categories[as.numeric(age_cat)])  %>% 
  rowwise() %>% 
  mutate(tmid = mean.Date(c(tl, tr))) %>% 
  ungroup()

# Save for figures
save(
  list = ls() %>% str_subset("prob_output") %>% c("marginal_ratio_draws", "marginal_ratio_stats"),
  file = here("generated_data/prob_output_bundle_for_figures.rdata") 
)

# E. Number exposed ----------------------------------------------------------

# Total number exposed during time period by epiweek
num_exposed <- tidybayes::spread_draws(genquant, 
                                       mu_output_total[i, j]) %>% 
  ungroup() %>% 
  rename(age_cat = j,
         period_id = i) %>% 
  mutate(age_cat = age_categories[as.numeric(age_cat)]) %>%
  {
    x <- .
    bind_rows(x,
              # Add post-stratified estimate
              x %>% 
                inner_join(sample_pop_age) %>% 
                group_by(period_id, .chain, .iteration, .draw) %>% 
                summarise(mu_output_total = weighted.mean(mu_output_total, n)) %>% 
                mutate(age_cat = "overall"))
  } %>% 
  inner_join(pop_dat) %>% 
  mutate(num_exposed = mu_output_total * pop)

# Compute statistics overal
tot_exposed <- num_exposed %>% 
  filter(age_cat == "overall") %>% 
  group_by(period_id, .draw) %>% 
  summarise(val = sum(num_exposed)) %>% 
  group_by(period_id) %>% 
  summarise(mean = mean(val),
            q5 = quantile(val, .05),
            q95 = quantile(val, .95)) %>% 
  ungroup() %>% 
  mutate(tl = output_periods$tl[period_id],
         tr = output_periods$tr[period_id]) %>%
  rowwise() %>%
  mutate(tmid = mean.Date(c(tl, tr))) %>%
  ungroup()

# Compute statistics by age
tot_exposed_age <- num_exposed %>% 
  group_by(period_id, .draw, age_cat) %>% 
  summarise(val = sum(mu_output_total*pop)) %>% 
  group_by(period_id, age_cat) %>% 
  summarise(mean = mean(val),
            q5 = quantile(val, .05),
            q95 = quantile(val, .95)) %>% 
  ungroup() %>% 
  mutate(tl = output_periods$tl[period_id],
         tr = output_periods$tr[period_id]) %>%
  rowwise() %>%
  mutate(tmid = mean.Date(c(tl, tr))) %>%
  ungroup()


# Total number exposed during time period overall
num_exposed_overall <- tidybayes::spread_draws(genquant_overall, 
                                               mu_output_total[i, j]) %>% 
  ungroup() %>% 
  rename(age_cat = j,
         period_id = i) %>% 
  mutate(age_cat = age_categories[as.numeric(age_cat)]) %>% {
    x <- .
    bind_rows(x,
              # Add post-stratified estimate
              x %>% 
                inner_join(sample_pop_age) %>% 
                group_by(period_id, .chain, .iteration, .draw) %>% 
                summarise(mu_output_total = weighted.mean(mu_output_total, n)) %>% 
                mutate(age_cat = "overall"))
  } %>% 
  inner_join(pop_dat) %>% 
  mutate(num_exposed = mu_output_total * pop)

# Compute statistics
num_exposed_overall_stats <- num_exposed_overall %>% 
  group_by(age_cat) %>%
  mutate(val = num_exposed) %>% 
  summarise(mean = mean(val),
            q5 = quantile(val, .05),
            q95 = quantile(val, .95))

# Length of modeling period
poi_dt <- getOutputPeriodOverall() %>% 
  unlist() %>% 
  diff()

# Compute annualized statistics
tot_exposed_age_cum_annualized <- num_exposed_overall %>% 
  group_by(age_cat) %>% 
  mutate(val = 365.25*(mu_output_total)/poi_dt) %>% 
  summarise(mean = mean(val),
            q5 = quantile(val, .05),
            q95 = quantile(val, .95)) %>% 
  ungroup()

# Save for figures
save(
  list = ls() %>% str_subset("exposed") %>% str_subset("num_exposed$", negate = T),
  file = here("generated_data/num_exposed_bundle_for_figures.rdata") 
)

# F. Summary table -----------------------------------------------------------

clinical_dt <- diff(range(date_subset)) %>% as.numeric()

# Surveillance data
surv_data <- lab_data %>% 
  # filter(date %in% study_period_dates) %>%
  {
    x <- .
    bind_rows(x, 
              x %>% mutate(age_cat = "overall"))
  } %>% 
  mutate(age_cat = age_cat_dict[age_cat]) %>% 
  group_by(age_cat) %>% 
  summarise(tot_cases = n(),
            tot_rtd_pos = sum(rdt_res)) %>% 
  pivot_longer(cols = contains("tot"),
               names_to = "what") %>% 
  inner_join(pop_dat %>% mutate(age_cat = age_cat_dict[age_cat])) %>% 
  mutate(value = 365.25*value/pop*1e3/clinical_dt) %>% 
  mutate(where = "Clinics only") %>% 
  select(age_cat, where, what, value)

# Estimated data
est_data <-
  bind_rows(
    # Incidence data
    bind_rows(tot_I_sum_age,
              tot_I_sum) %>% 
      mutate(what = "tot_clin_est",
             age_cat = age_cat_dict,
             where = "Clinics only") %>% 
      bind_rows(
        bind_rows(tot_I_sum_age_true,
                  tot_I_sum_true) %>% 
          mutate(what = "tot_comm_est",
                 age_cat = age_cat_dict,
                 where = "Full community")
      ),
    # Seroinfections data
    num_exposed_overall_stats %>% 
      mutate(age_cat = age_cat_dict[age_cat],
             what = "infections",
             where = "Full community")
  ) %>% 
  inner_join(pop_dat %>% mutate(age_cat = age_cat_dict[age_cat])) %>% 
  mutate(dt = ifelse(str_detect(what, "_est"), clinical_dt, poi_dt)) %>% 
  # Annualize
  mutate(across(c("mean", "q5", "q95"), ~ .x *365.25/pop*1e3/dt)) %>% 
  mutate(
    value = str_c(
      formatC(mean, format = "f", digits = 1),
      " (",
      formatC(q5, format = "f", digits = 1),
      "-",
      formatC(q95, format = "f", digits = 1),
      ")"
    ),
    age_cat = factor(age_cat, levels = age_cat_dict)) %>% 
  select(age_cat, where, what, value)

# Incidence calculation data
incidence_data <- lab_data %>% 
  {
    x <- .
    bind_rows(x, 
              x %>% mutate(age_cat = "overall"))
  } %>% 
  mutate(age_cat = age_cat_dict[age_cat]) %>% 
  group_by(age_cat) %>% 
  summarise(tot_cases = n(),
            tot_rtd_pos = sum(rdt_res)) %>% 
  pivot_longer(cols = contains("tot"),
               names_to = "what") %>% 
  inner_join(pop_dat %>% mutate(age_cat = age_cat_dict[age_cat])) %>% 
  mutate(value_annual = 365.25*value/clinical_dt) %>% 
  mutate(where = "Clinics only",
         annual_fraction = 365.25/clinical_dt,
         value = as.character(value),
         value_annual = as.character(value_annual)) %>%
  bind_rows(., bind_rows(
    # Incidence data
    bind_rows(tot_I_sum_age,
              tot_I_sum) %>% 
      mutate(what = "tot_clin_est",
             age_cat = age_cat_dict,
             where = "Clinics only") %>% 
      bind_rows(
        bind_rows(tot_I_sum_age_true,
                  tot_I_sum_true) %>% 
          mutate(what = "tot_comm_est",
                 age_cat = age_cat_dict,
                 where = "Full community")
      ),
    # Seroinfections data
    num_exposed_overall_stats %>% 
      mutate(age_cat = age_cat_dict[age_cat],
             what = "infections",
             where = "Full community")
  ) %>% 
    inner_join(pop_dat %>% mutate(age_cat = age_cat_dict[age_cat])) %>% 
    mutate(dt = ifelse(str_detect(what, "_est"), clinical_dt, poi_dt),
           value = str_c(
             formatC(mean, format = "f", digits = 1),
             " (",
             formatC(q5, format = "f", digits = 1),
             "-",
             formatC(q95, format = "f", digits = 1),
             ")"
           )) %>% 
    # Annualize
    mutate(across(c("mean", "q5", "q95"), ~ .x *365.25/dt)) %>% 
    mutate(
      value_annual = str_c(
        formatC(mean, format = "f", digits = 1),
        " (",
        formatC(q5, format = "f", digits = 1),
        "-",
        formatC(q95, format = "f", digits = 1),
        ")"
      ),
      age_cat = factor(age_cat, levels = age_cat_dict),
      annual_fraction = 365.25/dt) %>%
    select(age_cat, what, value, pop, value_annual, where, annual_fraction)
            )

# Save for figures
save(
  list = c("surv_data", "est_data", "incidence_data"),
  file = here(str_glue("generated_data/table_data_bundle_for_figures_{diarrhea_severity}.rdata")) 
)


# G. Compute state transitions -----------------------------------------------

# Fraction of symptomatic infections in community to infections
frac_symp_age <- bind_rows(tot_I_sum_true_draws %>% mutate(age_cat = "overall"),
                           tot_I_sum_true_age_draws) %>% 
  inner_join(num_exposed_overall) %>% 
  mutate(frac_symp = (val*poi_dt)/(num_exposed*clinical_dt)) 

num_symp_age_stats <- frac_symp_age %>% 
  group_by(age_cat) %>% 
  summarise(mean = mean(1/frac_symp),
            q5 = quantile(1/frac_symp, .05),
            q95 = quantile(1/frac_symp, .95))

frac_symp_age_stats <- frac_symp_age %>% 
  group_by(age_cat) %>% 
  summarise(mean = mean(frac_symp),
            q5 = quantile(frac_symp, .05),
            q95 = quantile(frac_symp, .95))

# Fraction of clinical cases to infections
frac_case_age <- bind_rows(tot_I_sum_draws %>% mutate(age_cat = "overall"),
                           tot_I_sum_age_draws) %>% 
  inner_join(num_exposed_overall) %>% 
  mutate(frac_case = (val*poi_dt)/(num_exposed*clinical_dt)) 

num_case_age_stats <- frac_case_age %>% 
  group_by(age_cat) %>% 
  summarise(mean = mean(1/frac_case),
            q5 = quantile(1/frac_case, .05),
            q95 = quantile(1/frac_case, .95))

frac_case_age_stats <- frac_case_age %>% 
  group_by(age_cat) %>% 
  summarise(mean = mean(frac_case),
            q5 = quantile(frac_case, .05),
            q95 = quantile(frac_case, .95))


# Fraction of positive culture to infections
# Using the posterior of culture sens/spec
culture_perf <- res_incid$draws(c("sens_all[4]", "spec_all[4]")) %>% 
  as_draws_df() %>% 
  as_tibble() %>% 
  rename(culture_sens = `sens_all[4]`,
         culture_spec = `spec_all[4]`)

# Compute fractions
frac_culture_age <- bind_rows(tot_I_sum_draws %>% mutate(age_cat = "overall"),
                              tot_I_sum_age_draws) %>% 
  inner_join(culture_perf) %>% 
  mutate(culture_sim = val * (culture_sens + (1-culture_spec))) %>% 
  inner_join(num_exposed_overall) %>% 
  mutate(frac_culture = (culture_sim*poi_dt)/(num_exposed*clinical_dt)) 

num_culture_age_stats <- frac_culture_age %>% 
  group_by(age_cat) %>% 
  summarise(mean = mean(1/frac_culture),
            q5 = quantile(1/frac_culture, .05),
            q95 = quantile(1/frac_culture, .95))

frac_culture_age_stats <- frac_culture_age %>% 
  group_by(age_cat) %>% 
  summarise(mean = mean(frac_culture),
            q5 = quantile(frac_culture, .05),
            q95 = quantile(frac_culture, .95))


# Fraction seeking healthcare by symptomatic infection
frac_seek_age <- healthprob_draws %>% 
  mutate(age_cat = age_categories[age_cat]) %>% 
  inner_join(pop_dat) %>% 
  {
    x <- .
    bind_rows(x, 
              x %>% 
                group_by(.draw) %>% 
                summarise(p = weighted.mean(p, pop),
                          pop = sum(pop)) %>% 
                mutate(age_cat = "overall"))
  } 

num_seek_age_stats <- frac_seek_age %>% 
  group_by(age_cat)  %>% 
  summarise(mean = mean(1/p),
            q5 = quantile(1/p, .05),
            q95 = quantile(1/p, .95))

frac_seek_age_stats <- frac_seek_age %>% 
  group_by(age_cat)  %>% 
  summarise(mean = mean(p),
            q5 = quantile(p, .05),
            q95 = quantile(p, .95))

# Save for figures
save(
  list = ls() %>% str_subset("frac|num") %>% str_subset("exposed", negate = T),
  file = here(str_glue("generated_data/transition_numbers_bundle_for_figures_{diarrhea_severity}.rdata"))
)


# H. Proportion in each case ----------------------------------------------

prob_case <- prob_case %>% 
  group_by(id) %>% 
  mutate(case_prob = exp(mean - log_sum_exp(mean))) 

prob_case %>% 
  group_by(case_name) %>% 
  summarise(sum = sum(case_prob),
            tot = n(),
            frac = sum/tot) %>% 
  mutate(n_infections = case_when(case_name %in% c("<1,1,1>", "<0,1,1>") ~ 2,
                                  case_name %in% c("<0,0,0>", "<1,0,0>") ~ 0,
                                  T ~ 1)) %>% 
  group_by(n_infections) %>% 
  summarise(frac = sum(frac))
