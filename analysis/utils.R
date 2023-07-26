# Filename formatting -----------------------------------------------------

getTimestamp <- function() {
  Sys.time() %>% 
    as.character() %>% 
    str_remove_all("\\:") %>% 
    str_replace(" ", "-")
}

makeStanDataFile <- function(opt) {
  str_glue("generated_data/stan_data_vibriocidal_agegrp{opt$age_grp}_J{opt$num_draws}_ogawa.rds") 
}

makeStanOutputName <- function(opt) {
  str_c(here("generated_data/stan_output_"), opt$output_prefix, ".rds")
}

makeStanOutputStatsName <- function(opt) {
  str_c(here("generated_data/stan_output_stats_"), opt$output_prefix, ".rds")
}

makeStanOutputCSVBasename <- function(opt) {
  str_c(here("generated_data/draws_"), opt$output_prefix, "")
}

makeGenquantOutputName <- function(opt) {
  str_c(here("generated_data/genquant_output_"), opt$output_prefix, ".rds")
}

makeGenquantOutputStatsName <- function(opt) {
  str_c(here("generated_data/genquant_output_stats_"), opt$output_prefix, ".rds")
}

updateOptPrefix <- function(opt) {
  opt$output_prefix <- str_c(opt$output_prefix, 
                             opt$stan_data %>% 
                               str_extract("(?<=vibriocidal)(.)*(?=\\.)"))
  opt
}

getAgeGroupFromFilename <- function(str) {
  str_extract(str, "(?<=agegrp)[1-3]") %>% as.numeric()
}

# Stan data preparation ---------------------------------------------------

computeCounts <- function(x, y) {
  x <- x %>% arrange(desc(rdt_res), desc(pcr_res), desc(cult_res))
  dat_A <- x %>% 
    filter(is.na(cult_res), is.na(pcr_res))
  
  dat_B <- x %>% 
    filter(is.na(cult_res), !is.na(pcr_res))
  
  counts_B <- dat_B %>% 
    count(rdt_res, pcr_res) %>% 
    complete(rdt_res = F, pcr_res = c(F, T)) %>% 
    replace_na(list(n = 0)) %>% 
    arrange(desc(rdt_res), desc(pcr_res))
  
  dat_C <- x %>% 
    filter(!is.na(cult_res), !is.na(pcr_res))
  
  counts_C <- dat_C %>% 
    count(rdt_res, pcr_res, cult_res) %>% 
    complete(rdt_res = c(T), pcr_res = c(F, T), cult_res = c(F, T)) %>% 
    replace_na(list(n = 0)) %>% 
    arrange(desc(rdt_res), desc(pcr_res), desc(cult_res))
  
  tibble(
    y_A = nrow(dat_A),
    y_B = list(counts_B$n),
    y_C = list(counts_C$n)
  )
}


computeIntervals <- function(t_left, 
                             t_right, 
                             n_intervals) {
  
  dateseq <- seq.Date(t_left, t_right, length.out = n_intervals + 1)
  tr <- dateseq[-1]
  tl <- dateseq[-length(dateseq)]
  
  rbind(tl, tr)
}


computeIntervalProbs <- function(intervals,
                                 cases,
                                 case_col = "cases_smooth") {
  
  # Get all dates in interval
  res <- cases %>%
    filter(date >= min(intervals["tl", ]),
           date < max(intervals["tr", ])) %>% 
    rowwise() %>% 
    mutate(interval = which(date >= intervals["tl", ] & date < intervals["tr", ])) %>% 
    group_by(interval) %>%
    summarise(sum_cases = sum(!!rlang::sym(case_col)),
              date_left = min(date),
              date_right = max(date),
              median_date = median(date)) %>%
    ungroup() %>% 
    # Compute probability
    mutate(tot_cases = sum(sum_cases),
           prob = sum_cases/tot_cases)
  
  if (abs(1 - sum(res$prob)) > 1e-4)
    stop("Probabilities do not sum to 1")
  
  return(res$prob)
}

computeSeqMidpoints <- function(t_left, 
                                t_right, 
                                n_intervals) {
  
  # Compute midpoints
  mid <- computeIntervals(t_left = t_left,
                          t_right = t_right,
                          n_intervals = n_intervals) %>% 
    colMeans() %>%
    as.Date(origin = "1970-01-01")
  
  mid
}

computeSeqMidpointsDT <- function(t_left, t_right, n_intervals) {
  
  computeSeqMidpoints(t_left = t_left,
                      t_right = t_right,
                      n_intervals = n_intervals) %>% 
    difftime(t_right, ., units = "days") %>% 
    as.numeric()
}

makeFullDates <- function(t_min = getTmin(),
                          t_max = getTmax()) {
  seq.Date(t_min, t_max, by = "1 days")
}

getTmin <- function() {
  as.Date("2021-03-27") - 180 
}

getTmax <- function() {
  as.Date("2022-02-13")
}

# For RDT reanalysis
getTminRDT <- function() {
  as.Date("2021-01-24") 
}

getTmaxRDT <- function() {
  as.Date("2022-08-28")
}
getStanDraws <- function(opt) {
  dir(here("generated_data"), 
      pattern = opt$output_prefix, 
      full.names = T) %>% 
    str_subset("draws")
}

# From Kirsten
get_age_group <- function(age, 
                          num_breaks = 2,
                          break1 = 5, 
                          break2 = 65) {
  # 3 age categories
  if (num_breaks == 2) {
    grp <- ifelse(age < break1, 1,
                  ifelse(age >= break1 & age < break2, 2, 
                         3)
    )
    # 2 age categories
  } else if (num_breaks == 1) {
    grp <- ifelse(age < break1, 1, 2)
    # error if different number of categories requested
  } else {
    stop('num_breaks must be 1 or 2')
  }
  return(grp)
} 

extract_one_draw <- function(stanfit, chain = 1, iter = 1) {
  x <- get_inits(stanfit, iter = iter)
  x[[chain]]
}


periodDataToMatrix <- function(df,
                               date_col,
                               val_col) {
  
  if (!(val_col %in% colnames(df))) {
    stop("Colname ", val_col, " not in df.")
  }
  
  df %>% 
    arrange(all_of(date_col)) %>% 
    select(all_of(date_col), age_cat, all_of(val_col)) %>% 
    pivot_wider(names_from = "age_cat",
                values_from = val_col) %>% 
    select(-all_of(date_col)) %>% 
    as.matrix()
}


periodDataToArray <- function(df,
                              date_col,
                              val_col) {
  
  if (!(val_col %in% colnames(df))) {
    stop("Colname ", val_col, " not in df.")
  }
  
  group_levels <- unique(df$age_cat)
  
  map(group_levels, function(x) {
    df %>% 
      arrange(all_of(date_col)) %>% 
      filter(age_cat == x) %>% 
      pull(all_of(val_col)) %>% 
      Reduce(f = rbind) %>% 
      as.matrix()
  }) %>% 
    abind::abind(along = 3)
}


getWeeklyOutputPeriods <- function() {
  tibble(
    tl = seq.Date(getTmin(), getTmax() - 7, by = "1 week"),
    tr = seq.Date(getTmin() + 7, getTmax(), by = "1 week")
  )
}

makeGenOutputPeriodData <- function(stan_data_output,
                                    output_periods,
                                    ref_dates) {
  
  stan_data_output$N_output <- nrow(output_periods)
  stan_data_output$foi_tl_output <- map_dbl(output_periods$tl, ~ which(ref_dates == .))
  stan_data_output$foi_tr_output <- map_dbl(output_periods$tr, ~ which(ref_dates == .))
  
  stan_data_output
}


# Result processing -------------------------------------------------------

extractCountStats <- function(genquant, 
                              period_data, 
                              variable,
                              age_categories) {
  
  genquant$summary(variables = variable, .cores = 4) %>% 
    mutate(time = str_extract(variable, "(?<=\\[)[0-9]+(?=,)") %>% as.numeric(),
           age_cat_num = str_extract(variable, "(?<=,)[0-9]+(?=\\])") %>% as.numeric(),
           age_cat = age_categories[age_cat_num],
           case = str_c("counts_", str_extract(variable, "(?<=,)[0-9]+(?=,)")),
           date = period_data$tl[time],
           epiweek = period_data$epiweek[time],
           var = str_extract(variable, "(.)*(?=\\[)"))
}



# Data processing  --------------------------------------------------------

get_age_group <- function(age, 
                          num_breaks = 2,
                          break1 = 5, 
                          break2 = 65) {
  # 3 age categories
  if (num_breaks == 2) {
    grp <- ifelse(age < break1, 1,
                  ifelse(age >= break1 & age < break2, 2, 
                         3)
    )
    # 2 age categories
  } else if (num_breaks == 1) {
    grp <- ifelse(age < break1, 1, 2)
    # error if different number of categories requested
  } else {
    stop('num_breaks must be 1 or 2')
  }
  return(grp)
}


getHealthSeekingData <- function(severity = "diarrhea2") {
  
  # load household care seeking behavior data for round 1
  files <- dir(here::here("serochit_data"), 
               pattern = "baseline_individual_survey",
               full.names = T)
  
  # Define ind and cat columns
  ind_col <- str_glue("ind_{severity}_care")
  site_col <- str_glue("cat_{severity}_site")
  
  healthcare_seek <- map_df(
    files,
    function(x) {
      dat <- read_csv(x)
      dat %>% 
        dplyr::select(round, ID, hh_ID, n_age_years, one_of(c(ind_col, site_col))) %>% 
        mutate(round = as.character(round))
    })
  
  # add variables and filter missing data
  healthcare_seek <- healthcare_seek %>%
    dplyr::filter(!is.na(!!rlang::sym(ind_col))) %>%
    dplyr::mutate(seek_clinic_care = ifelse(!!rlang::sym(ind_col) == 'yes' &
                                              (!!rlang::sym(site_col) %in% c('bitid_hosp',  'sitakunda_uhc')),
                                            1, 0),
                  age_group = get_age_group(n_age_years))
  
  # number that seek care by age group
  df_cs <- healthcare_seek %>%
    group_by(age_group) %>%
    mutate(count = ifelse(!is.na(!!rlang::sym(ind_col)), 1, 0)) %>%
    summarise(x = sum(seek_clinic_care == 1),
              n = sum(count))
  
  df_cs
}

getHealthSeekingDraws <- function(hs_data = getHealthSeekingData(),
                                  hs_stan_model = "analysis/stan/binomial_symp_reporting.stan") {
  
  healthprob_model <- cmdstan_model(hs_stan_model)
  
  healthprob_fit <- healthprob_model$sample(data = list(N = nrow(hs_data),
                                                        y = hs_data$x,
                                                        count = hs_data$n),
                                            iter_warmup = 250,
                                            iter_sampling = 1250,
                                            refresh = 500)
  
  healthprob_draws <- healthprob_fit %>%
    tidybayes::spread_draws(p[i]) %>%
    ungroup() %>% 
    select(age_cat = i, p, .draw)
  
  healthprob_draws
}


addEpiWeek <- function(df, date_col = "date") {
  
  df %>% 
    mutate(week = lubridate::epiweek(!!rlang::sym(date_col)),
           year = lubridate::epiyear(!!rlang::sym(date_col)),
           epiweek = str_c(year, week, sep = "-"))
}

aggregateIncid <- function(I_draws,
                           chain_subset = NULL,
                           draw_subset = NULL,
                           date_subset = NULL,
                           by_age = FALSE,
                           age_categories = NULL,
                           weekly = TRUE,
                           cumulative = FALSE,
                           do_stats = TRUE
) {
  
  filtered_df <- I_draws %>%
    ungroup() %>%
    # Subset draws
    {
      if (!is.null(draw_subset)) {
        filter(., .draw %in% draw_subset)
      } else {
        .
      }
    } %>%
    # Subset Chains
    {
      if (!is.null(chain_subset)) {
        filter(., .chain %in% chain_subset)
      } else {
        .
      }
    } %>% 
    # Subset date
    {
      if (!is.null(date_subset)) {
        filter(., date %in% date_subset)
      } else {
        .
      }
    } 
  
  grouped_df <- filtered_df %>% 
    {
      x <- filtered_df
      if (weekly) {
        y <- addEpiWeek(x) 
        {
          if (by_age) {
            group_by(y, .draw, epiweek, age_cat) 
          } else {
            group_by(y, .draw, epiweek) 
          }
        }
      } else {
        if (by_age) {
          group_by(x, .draw, age_cat)
        } else {
          group_by(x, .draw)
        }
      }
    } %>%
    summarise(
      tl = min(date, na.rm = T),
      tr = max(date, na.rm = T),
      val = sum(value, na.rm = T)
    ) %>%
    ungroup() 
  
  if (do_stats) {
    res_df <- grouped_df %>% 
      {
        if (weekly) {
          # Cumulative incidence
          {
            if (cumulative) {
              {
                if (by_age) {
                  group_by(., .draw, age_cat) %>% 
                    arrange(tl, age_cat) 
                } else {
                  group_by(., .draw)  %>% 
                    arrange(tl) 
                }
              } %>% 
                mutate(val = cumsum(val)) %>% 
                ungroup()
            } else {
              .
            }
          } %>% 
            {
              if (by_age) {
                group_by(., epiweek, tl, tr, age_cat)
              } else {
                group_by(., epiweek, tl, tr)
              }
            }
        } else {
          if (by_age) {
            group_by(., tl, tr, age_cat)
          } else {
            group_by(., tl, tr)
          }
        }
      } %>%
      summarise(mean = mean(val),
                q5 = quantile(val, .05),
                q95 = quantile(val, .95)) %>%
      ungroup() 
    
  } else {
    res_df <- grouped_df
  }
  
  res_df %>% 
    rowwise() %>%
    mutate(tmid = mean.Date(c(tl, tr))) %>%
    ungroup() %>%
    {
      if (by_age) {
        mutate(., age_cat = age_categories[as.numeric(age_cat)])
      } else {
        .
      }
    }
}

getOutputPeriodOverall <- function() {
  tibble(
    tl = c(as.Date("2021-03-27")),
    tr = c(as.Date("2022-02-13"))
  )
}

getWeeklyOutputPeriods <- function() {
  tibble(
    tl = seq.Date(getTmin(), getTmax() - 7, by = "1 week"),
    tr = seq.Date(getTmin() + 7, getTmax(), by = "1 week")
  )
}

getSitakundatData <- function(data_path = "serochit_data/census_data_sex_age_sitakunda.csv",
                              age_categories) {
  read_csv(data_path) %>% 
    janitor::clean_names() %>% 
    select(age, pop = both_sexes_population) %>% 
    mutate(age = case_when(age == "80+" ~ "80-100", 
                           T ~ age)) %>% 
    mutate(al = str_extract(age, "[0-9]+(?=-)") %>% as.numeric(),
           ar = str_extract(age, "(?<=-)[0-9]+") %>% as.numeric(),
           pop_2021 = pop * 1.015^10,
           pop_2021 = case_when(age == "0-4" ~ pop_2021 * .8, 
                                T ~ pop_2021),
           age_cat = get_age_group(al)) %>% 
    group_by(age_cat) %>% 
    summarise(pop = sum(pop_2021)) %>% 
    ungroup() %>% 
    mutate(age_cat = age_categories[age_cat])
}


# Other utils -------------------------------------------------------------

logit <- function(x) {
  log(x/(1-x))
}

print_file <- function(file) {
  cat(paste(readLines(file), "\n", sep=""), sep="")
}

# From the rethinking package
log_sum_exp <- function(x) {
  xmax <- max(x)
  xsum <- sum(exp(x - xmax))
  xmax + log(xsum)
}

# Figures -----------------------------------------------------------------

getAgeCatDict <- function() {
  age_cat_dict <- c("< 5", "5-64", "65+", "overall")
  names(age_cat_dict) <- c("[0,5)", "[5,65)", "[65,Inf)", "overall")
  age_cat_dict
}

getAgeCatDictFigures <- function() {
  age_cat_dict <- c("1-4", "5-64", "65+", "Overall")
  names(age_cat_dict) <- c("[0,5)", "[5,65)", "[65,Inf)", "overall")
  age_cat_dict
}
