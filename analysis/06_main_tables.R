library(tidyverse)
#library(serochit.data.access)
library(sf)
source(here("analysis/utils.R"))


gtsummary::theme_gtsummary_journal(
    journal = c("jama", "lancet", "nejm", "qjecon")[3],
    set_theme = TRUE
)


# clinical data table --------------------------------------------------

clinical <- read_csv("serochit_data/clinical_surveillance_vars.csv") %>%
    dplyr::left_join(.,read_csv("serochit_data/clinical_24jan2021_13feb2022_sitakunda.csv") %>%
                         dplyr::select(record_id,sitakunda,cat_dehydration_status,duration))

clinical %>% 
    dplyr::filter(sitakunda=="Yes") %>%
    mutate(age_group = case_when(n_criteria_age<5 ~ "<5", 
                                 n_criteria_age>=5 & n_criteria_age < 65 ~ "5-64", 
                                 n_criteria_age>=65 ~ "65+") %>% factor(levels = c("<5", "5-64", "65+")),
           cat_dehydration_status = case_when(cat_dehydration_status==0 ~ "None",
                                              cat_dehydration_status==1 ~ "Some",
                                              cat_dehydration_status==2 ~ "Severe") %>% factor(levels = c("None", "Some", "Severe"))) %>%
    dplyr::select(`Site`=site, 
                  `Sex` = cat_sex,
           `Age` = n_criteria_age, 
           `Age group` = age_group, 
           `Patient Type`=cat_care_type, 
           `Dehydration status`=cat_dehydration_status,
           `Culture` = ind_cul_result, 
           `PCR` = ind_pcr_result,
           # `Reported spending most of their time in Sitakunda in the last 7 days` = ind_time_sitakunda,
           rdt_result, 
           `Self-reported antibiotic use 48 hours prior to admission` = antibiotic_use, 
           `No. of days between symptom onset and seeking care` = days_symptom_care,
           `No. of days between healthcare facility admission and discharge` = duration) %>% 
    gtsummary::tbl_summary(by = "rdt_result") %>%
    gtsummary::add_p()


# serological data and survey data table --------------------------------------------------

survey <- bind_rows(
    read_csv("serochit_data/baseline_individual_survey_careseek1a.csv"),
    read_csv("serochit_data/baseline_individual_survey_careseek1b.csv"),
    read_csv("serochit_data/baseline_individual_survey_careseek2.csv") %>% mutate(round = as.character(round)) %>% select(-txt_injury1_trans),
    read_csv("serochit_data/baseline_individual_survey_careseek3.csv") %>% mutate(round = as.character(round)) %>% select(-txt_injury1_trans)
) 
ids <- read_csv("serochit_data/kav_unique_ids.csv")

grid.shp = read_csv("serochit_data/sampling_grid.csv")

# Shapefile of Bangladesh and clusters
bang.map.3 = readRDS("serochit_data/BGD_adm3.rds")
#cluster.shp = readRDS(here::here("data/serosurvey/cluster_shp.rds"))
cluster_pop.shp = readRDS("serochit_data/cluster_shp_pop.rds")

#pop.raster = raster::raster(here::here("data/rasters_100m/BGD_ppp_2020_adj_v2.tif"))
bgd_raster <- bang.map.3 %>% 
    st_transform(., crs = "+proj=longlat +datum=WGS84 +no_defs") %>% 
    as_Spatial() %>%
    malariaAtlas::getRaster(
        surface = "Global friction surface enumerating land-based travel speed with access to motorized transport for a nominal year 2019", 
        shp = .)

bgd_values <- cluster_pop.shp %>% 
    dplyr::mutate(urbanicity2 = exactextractr::exact_extract(bgd_raster, cluster_pop.shp, 'mean'), 
                  area = units::set_units(st_area(geometry), km^2), 
                  pop_density2 = pop_sum/area) %>%
    dplyr::select(id_grid,urbanicity2, pop_density2, pop_sum, area) %>%
    as_tibble()

quant_labels <- function(x){
    temp <- list()
    
    for(i in 2:length(x)){
        temp[[i]] <- paste0(if_else(x[i-1] < 1, signif(x[i-1], 2), round(x[i-1])), 
                            '-', 
                            if_else(x[i] < 1, signif(x[i], 2), round(x[i])))
    }
    
    unlist(temp)
}

quant_levels <- function(x){
    c(min(as.numeric(x)), quantile(as.numeric(x), probs = c(0.25, 0.5, 0.75), na.rm=TRUE), max(as.numeric(x)*1.01))
}

covid <- survey %>% 
    mutate(cat_sex = factor(cat_sex, levels = c("male", "female"), 
                            labels = c("Male", "Female")),
           age_group = case_when(n_age<5 ~ "<5", 
                                 n_age>=5 & n_age <65~ "5-64", 
                                 n_age>=65 ~ "65+") %>% factor(levels = c("<5", "5-64", "65+")),
           # age_group = cut(n_age,c(seq(0,15,by=5), seq(25, 65, by = 10), 100),right = FALSE) %>% ordered, 
           # age_group = factor(age_cat, levels = pop_levels, labels = pop_labels),
           cat_activity_month = if_else(cat_activity_month=="decline", NA_character_, cat_activity_month), 
           cat_activity_month = factor(cat_activity_month, 
                                       levels = c("business", "child", "farmer", "homeworker", "none", "other", "outside", "student"), 
                                       labels = c("Business at home", "Child", "Farmer", "Homemaker", "No work (adult)", "Other","Business outside home", "Student")),
           cat_edu_attainment = if_else(ind_edu_school=="no", "none", cat_edu_attainment), 
           cat_edu_attainment = factor(cat_edu_attainment, 
                                       levels = c("none", "primary", "lower_secondary", "upper_secondary", "bachelors", "post_graduate"), 
                                       labels = c("No schooling", "Primary", "Lower secondary", "Upper secondary", "Bachelors", "Postgraduate")),
           across(c("ind_diarrhea1_care", "ind_diarrhea2_care", "ind_diarrhea3_care"), ~factor(.x, levels = c("yes", "no"), labels = c("Yes", "No")))) %>%
    select(hh_ID, ID, n_age, age_group, cat_sex, cat_activity_month, cat_edu_attainment, ind_diarrhea1_care, ind_diarrhea2_care, ind_diarrhea3_care) %>% 
    left_join(grid.shp %>% select(hh_ID=id_structure, id_grid)) %>%
    mutate(id_grid = if_else(is.na(id_grid) & hh_ID == "H000001", 21L, id_grid)) %>%
    left_join(bgd_values)

covid <- covid %>%
    mutate(pop_density = cut(as.numeric(pop_density2), 
                             breaks = quant_levels(covid$pop_density2[!is.na(covid$pop_density2)]), 
                             labels = quant_labels(quant_levels(covid$pop_density2[!is.na(covid$pop_density2)])),
                             include.lowest = TRUE, right = FALSE),
           urbanicity = cut(urbanicity2, 
                            breaks = quant_levels(covid$urbanicity2[!is.na(covid$urbanicity2)]),
                            labels = quant_labels(quant_levels(covid$urbanicity2[!is.na(covid$urbanicity2)])),
                            include.lowest = TRUE, right = FALSE))

covid <- covid %>% 
    mutate(num_rounds = if_else(ID %in% ids$id, "Serological Data", "Survey"))

covid %>% 
    bind_rows(covid %>% filter(num_rounds=="Serological Data") %>% mutate(num_rounds = "Survey")) %>% 
    dplyr::select(num_rounds, 
                  
                  `Sex` = cat_sex, 
           `Age` = n_age, 
           `Age group` = age_group, 
           `% that would seek medical care if experiencing 3+ loose stools in a day`=ind_diarrhea1_care,
           `% that would seek medical care if experiencing 3+ loose stools in a day and dehydration`=ind_diarrhea2_care,
           `% that would seek medical care if experiencing 3+ loose stools in a day for >3 days`=ind_diarrhea3_care,
           `Main activity in previous month` = cat_activity_month, 
           `Highest educational attainment` = cat_edu_attainment, 
           Friction = urbanicity, 
           `Population density, per 1 km` = pop_density) %>% 
    gtsummary::tbl_summary(by = "num_rounds", 
                           type = list(c("Population density, per 1 km") ~ "categorical")) %>%
    gtsummary::add_p()


# summary data output table --------------------------------------------------

age_cat_dict <- getAgeCatDictFigures() 

load(here("generated_data/table_data_bundle_for_figures.rdata"))

var_dict <- c(
    tot_cases = "Suspected cholera cases",
    tot_rtd_pos = "RDT positive suspected cases",
    tot_clin_est = "True cholera incidence (clinical)",
    tot_comm_est = "True cholera incidence (community + clinics)",
    infections = "Infections"
)

table_data <- bind_rows(
    surv_data %>% dplyr::mutate(value = formatC(value, format = "f", digits = 1),age_cat=recode(age_cat, "< 5"="1-4","overall"="Overall")),
    est_data %>% dplyr::mutate(age_cat=recode(age_cat, "< 5"="1-4","overall"="Overall"))
) %>% 
    mutate(what = var_dict[what],
           age_cat = factor(age_cat, levels = age_cat_dict)) %>% 
    pivot_wider(names_from = "age_cat",
                values_from = "value")

table_data %>% 
    select(where, what, one_of(age_cat_dict)) %>% 
    kbl(caption = "Estimates of annualized incidence rates per 1,000 population") %>% 
    kable_styling() %>% 
    pack_rows(index = table(table_data$where))


