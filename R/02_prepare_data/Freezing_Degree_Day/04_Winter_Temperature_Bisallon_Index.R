



# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: WORKING DIRECTORY
# ══════════════════════════════════════════════════════════════════════════════

# ── Working directory ─────────────────────────────────────────────────────────
# Automatically detects the repository root based on the script location.
# No hardcoded path — works on any machine.
current_file <- rstudioapi::getActiveDocumentContext()$path
path_parts   <- unlist(strsplit(normalizePath(current_file, winslash = "/"), "/"))
target_index <- which(path_parts == "phd-chapter3-Atlantic-Salmon")
base_path    <- paste(path_parts[1:target_index], collapse = "/")
setwd(base_path)



# |========================|
# |     LOAD LIBRARIES     |
# |========================|

library(fst)
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(lubridate)
library(ggpubr)
library(readxl)
library(purrr)


# |=====================|
# |     IMPORT DATA     |
# |=====================|

# Air temperature data
sj_combined <- read_fst("02_PrepareData_Scripts/Intermediate_Cleaned_Data/Winter_Temperature/sj_air_temp_daily_mean_combined.fst")
tr_combined <- read_fst("02_PrepareData_Scripts/Intermediate_Cleaned_Data/Winter_Temperature/tr_air_temp_daily_mean_combined.fst")

# Smolt abundance data
prop_smolt <- as.data.frame(read_excel("01_RawData/MELCCFP_Data/BDFUSION_2024.xlsm", sheet = 1, n_max = 91))


# |===================================|
# |      ADAPT SALMON DATAFRAMES      |
# |===================================|

colnames(prop_smolt) <- iconv(colnames(prop_smolt),       from = "UTF-8",        to = "ASCII//TRANSLIT")
colnames(prop_smolt) <- gsub(pattern = "\\.",     replacement = "",      x = colnames(prop_smolt))
colnames(prop_smolt) <- gsub(pattern = "\\(",     replacement = "",      x = colnames(prop_smolt))
colnames(prop_smolt) <- gsub(pattern = "\\)",     replacement = "",      x = colnames(prop_smolt))
colnames(prop_smolt) <- gsub(pattern = " ",       replacement = "_",     x = colnames(prop_smolt))
colnames(prop_smolt) <- gsub(pattern = "\\%",     replacement = "PC",    x = colnames(prop_smolt))
colnames(prop_smolt) <- gsub(pattern = "\\+",     replacement = "P",     x = colnames(prop_smolt))
colnames(prop_smolt) <- gsub(pattern = "/",       replacement = "_Par_", x = colnames(prop_smolt))
colnames(prop_smolt) <- gsub(pattern = "-",       replacement = "_",     x = colnames(prop_smolt))
colnames(prop_smolt) <- gsub(pattern = "l'annee", replacement = "annee", x = colnames(prop_smolt))
colnames(prop_smolt) <- gsub(pattern = "oufs",    replacement = "oeufs", x = colnames(prop_smolt))
colnames(prop_smolt) <- gsub(pattern = "__",      replacement = "_",     x = colnames(prop_smolt))
prop_smolt$Riviere   <- iconv(prop_smolt$Riviere,         from = "UTF-8",        to = "ASCII//TRANSLIT")
colnames(prop_smolt) <- tolower(colnames(prop_smolt))
prop_smolt$riviere   <- tolower(prop_smolt$riviere)
prop_smolt$riviere <- ifelse(prop_smolt$riviere == "saint-jean", "stj", prop_smolt$riviere)
prop_smolt$riviere <- ifelse(prop_smolt$riviere == "trinite", "tri", prop_smolt$riviere)

# Subset columns
prop_smolt_sub <- prop_smolt %>%
  dplyr::select("riviere", "annee", "smolt_nb_en_devalaison_2_p", "smolt_nb_en_devalaison_3_p", "smolt_nb_en_devalaison_4_p", "smolt_nb_en_devalaison_5_p") %>%
  dplyr::rename("river"="riviere", "year"="annee", "smolt_nb_2yr"="smolt_nb_en_devalaison_2_p", "smolt_nb_3yr"="smolt_nb_en_devalaison_3_p",
         "smolt_nb_4yr"="smolt_nb_en_devalaison_4_p", "smolt_nb_5yr"="smolt_nb_en_devalaison_5_p") 

# Make sure all years are there
prop_smolt_sub <- prop_smolt_sub %>%
  dplyr::group_by(river) %>%
  # complete years between min et max
  complete(year = seq(min(year), max(year)), fill = list(
    smolt_nb_2yr = NA,
    smolt_nb_3yr = NA,
    smolt_nb_4yr = NA,
    smolt_nb_5yr = NA
  )) %>%
  ungroup()

# Create a dataset for each river, with a lag to know from what cohort each smolt come from, and then calculate the proportion of individuals from this cohort smoltified at each age (2-5yr)
prop_smolt_sub_sj <- prop_smolt_sub %>%
  dplyr::filter(river == "stj") %>%
  dplyr::arrange(year) %>%
  dplyr::mutate(
    cohort = year,                 
    smolt_prop_2yr = lead(smolt_nb_2yr, n = 3), # 2yr smolts correspond to the egg-laying cohort of t-3 (t being the year of smoltification).
    smolt_prop_3yr = lead(smolt_nb_3yr, n = 4), # 3yr smolts correspond to the egg-laying cohort of t-4 (t being the year of smoltification).
    smolt_prop_4yr = lead(smolt_nb_4yr, n = 5), # 4yr smolts correspond to the egg-laying cohort of t-5 (t being the year of smoltification).
    smolt_prop_5yr = lead(smolt_nb_5yr, n = 6)  # 5yr smolts correspond to the egg-laying cohort of t-6 (t being the year of smoltification).
  ) %>%
  dplyr::select(cohort, smolt_prop_2yr, smolt_prop_3yr, smolt_prop_4yr, smolt_prop_5yr) %>%
  rowwise() %>%
  dplyr::mutate(
    # Total with na.rm=TRUE only for 2yr and 5yr
    total = sum(c(smolt_prop_2yr, smolt_prop_5yr), na.rm = TRUE) +
      smolt_prop_3yr + smolt_prop_4yr,  # si NA ici, total devient NA
    smolt_prop_2yr = smolt_prop_2yr / total,
    smolt_prop_3yr = smolt_prop_3yr / total,
    smolt_prop_4yr = smolt_prop_4yr / total,
    smolt_prop_5yr = smolt_prop_5yr / total
  ) %>%
  ungroup() %>%
  dplyr::select(-total) %>%
  as.data.frame()

prop_smolt_sub_tr <- prop_smolt_sub %>%
  dplyr::filter(river == "tri") %>%
  dplyr::arrange(year) %>%
  dplyr::mutate(
    cohort = year,                 
    smolt_prop_2yr = lead(smolt_nb_2yr, n = 3), # 2yr smolts correspond to the egg-laying cohort of t-3 (t being the year of smoltification).
    smolt_prop_3yr = lead(smolt_nb_3yr, n = 4), # 3yr smolts correspond to the egg-laying cohort of t-4 (t being the year of smoltification).
    smolt_prop_4yr = lead(smolt_nb_4yr, n = 5), # 4yr smolts correspond to the egg-laying cohort of t-5 (t being the year of smoltification).
    smolt_prop_5yr = lead(smolt_nb_5yr, n = 6)  # 5yr smolts correspond to the egg-laying cohort of t-6 (t being the year of smoltification).
  ) %>%
  dplyr::select(cohort, smolt_prop_2yr, smolt_prop_3yr, smolt_prop_4yr, smolt_prop_5yr) %>%
  rowwise() %>%
  dplyr::mutate(
    # Total with na.rm=TRUE only for 2yr and 5yr
    total = sum(c(smolt_prop_2yr, smolt_prop_5yr), na.rm = TRUE) +
      smolt_prop_3yr + smolt_prop_4yr,  # si NA ici, total devient NA
    smolt_prop_2yr = smolt_prop_2yr / total,
    smolt_prop_3yr = smolt_prop_3yr / total,
    smolt_prop_4yr = smolt_prop_4yr / total,
    smolt_prop_5yr = smolt_prop_5yr / total
  ) %>%
  ungroup() %>%
  dplyr::select(-total) %>%
  as.data.frame()


# |============================================================================================================================================|
# |      CUMULATIVE FREEZING DEGREE DAY INDEX BETWEEN OCTOBER 1ST AND END OF EACH MONTH (NOVEMBER TO APRIL), BASED ON BISAILLON (2005)         |
# |============================================================================================================================================|

# According to Bisallon (2005)'s master theses
### Egg-fry mortality: Significant negative relationship between egg-fry mortality and FDDo-j (cumulated freezing degree days from October 1st to 
### January 1st) -> mortality decreases during cold winters.

# Then, the freezing dd is summed over each cohort and weighted by proportion of indivual leaving the river

# smolt_data = prop_smolt_sub_tr
# temperature_data = tr_combined
# temp_col = "Mean_T"

calc_dd <- function(smolt_data, temperature_data, temp_col){
  
  # Years vector, to manage in the loop later
  years <- unique(smolt_data$cohort)
  
  # Adding columns for periods between October and later months in winter
  smolt_data$co_fdd_per_oc_no <- NA
  smolt_data$co_fdd_per_oc_de <- NA
  smolt_data$co_fdd_per_oc_ja <- NA
  smolt_data$co_fdd_per_oc_fe <- NA
  smolt_data$co_fdd_per_oc_ma <- NA
  smolt_data$co_fdd_per_oc_ap <- NA
  
  # Adding monthly column
  smolt_data$co_fdd_mon_oc <- NA
  smolt_data$co_fdd_mon_no <- NA
  smolt_data$co_fdd_mon_de <- NA
  smolt_data$co_fdd_mon_ja <- NA
  smolt_data$co_fdd_mon_fe <- NA
  smolt_data$co_fdd_mon_ma <- NA
  smolt_data$co_fdd_mon_ap <- NA
  
  # Adding bimonthly column
  smolt_data$co_fdd_bimon_no_de <- NA
  smolt_data$co_fdd_bimon_ja_fe <- NA
  smolt_data$co_fdd_bimon_ma_ap <- NA
  
  # Adding bimonthly column
  smolt_data$co_d_minus_10_no <- NA
  smolt_data$co_d_minus_10_de <- NA
  
  # Preparation of the temperature dataset
  temperature_data$year <- year(temperature_data$date)
  temperature_data$month <- month(temperature_data$date)
  temperature_data$day <- day(temperature_data$date)
  temperature_data$doy <- yday(temperature_data$date)
  
  for(i in seq_along(years)){
    # Extract the proportions, replacing NA values with 0 to consider that this age doesn't contribute to the population for this cohort
    # Note that cohorts with NA (converted to 0) at 3+ and/or 4+ smolting year will be removed latter, since those ages are considered to contribute significently to the cohort and thus replacing NA by 0 at these ages would cause important biases
    smolt_data_i <- smolt_data %>% filter(cohort == years[i])
    p2 <- ifelse(is.na(smolt_data_i$smolt_prop_2yr), 0, smolt_data_i$smolt_prop_2yr)
    p3 <- ifelse(is.na(smolt_data_i$smolt_prop_3yr), 0, smolt_data_i$smolt_prop_3yr)
    p4 <- ifelse(is.na(smolt_data_i$smolt_prop_4yr), 0, smolt_data_i$smolt_prop_4yr)
    p5 <- ifelse(is.na(smolt_data_i$smolt_prop_5yr), 0, smolt_data_i$smolt_prop_5yr)
    
    # Convert temperatures to negative values (for those that are negative), and 0 if positive
    temperature_data$dd <- ifelse(temperature_data[,temp_col] < 0, temperature_data[,temp_col], 0)
    
    # Definition of periods between October and later months in winter
    periods <- c(
      "oc_no" = "11",
      "oc_de" = "12",
      "oc_ja" = "01",
      "oc_fe" = "02",
      "oc_ma" = "03",
      "oc_ap" = "04"
    )
    
    # Definition of monthly periods
    monthly_periods <- c(
      "oc" = "10",
      "no" = "11",
      "de" = "12",
      "ja" = "01",
      "fe" = "02",
      "ma" = "03",
      "ap" = "04"
    )
    
    # Definition of bimonthly periods
    bi_monthly_periods <- list(
      "no_de" = c("11", "12"),
      "ja_fe" = c("01", "02"),
      "ma_ap" = c("03", "04")
    )
    
    # Definition of frazil and anchor periods
    minus10_periods <- c(
      "no" = "11",
      "de" = "12"
      )
    
    # Function to extract a period between October and later months in winter
    extract_temp_period <- function(start_year, end_month, data) {
      start_date <- as.Date(paste0(start_year, "-10-01"))
      end_year <- ifelse(as.numeric(end_month) < 10, start_year + 1, start_year)
      end_date <- as.Date(paste0(end_year, "-", end_month, "-", days_in_month(as.Date(paste0(end_year, "-", end_month, "-01")))))
      data %>%
        dplyr::filter(date >= start_date & date <= end_date) %>%
        dplyr::pull(dd)
    }
    
    # Function to extract monthly periods
    extract_month <- function(start_year, month, data){
      end_year <- ifelse(as.numeric(month) < 10, start_year + 1, start_year)
      data %>%
        dplyr::filter(year == end_year & month == as.numeric(month)) %>%
        dplyr::pull(dd)
    }
    
    # Function to extract bimonthly periods (2 consecutive months)
    extract_two_months <- function(start_year, months, data){
      end_years <- ifelse(as.numeric(months) < 10, start_year + 1, start_year)
      data %>%
        dplyr::filter((year == end_years[1] & month == as.numeric(months[1])) |
                 (year == end_years[2] & month == as.numeric(months[2]))) %>%
        dplyr::pull(dd)
    }
    
    # Application of periods between October and later months in winter
    t_0 <- map(periods, ~ extract_temp_period(start_year = years[i]+1, end_month = .x, data = temperature_data))
    t_1 <- map(periods, ~ extract_temp_period(start_year = years[i]+2, end_month = .x, data = temperature_data))
    t_2 <- map(periods, ~ extract_temp_period(start_year = years[i]+3, end_month = .x, data = temperature_data))
    t_3 <- map(periods, ~ extract_temp_period(start_year = years[i]+4, end_month = .x, data = temperature_data))
    t_4 <- map(periods, ~ extract_temp_period(start_year = years[i]+5, end_month = .x, data = temperature_data))
    t_5 <- map(periods, ~ extract_temp_period(start_year = years[i]+6, end_month = .x, data = temperature_data))
    
    # Application of monthly periods
    m_0 <- map(monthly_periods, ~ extract_month(start_year = years[i]+1, month = .x, data = temperature_data))
    m_1 <- map(monthly_periods, ~ extract_month(start_year = years[i]+2, month = .x, data = temperature_data))
    m_2 <- map(monthly_periods, ~ extract_month(start_year = years[i]+3, month = .x, data = temperature_data))
    m_3 <- map(monthly_periods, ~ extract_month(start_year = years[i]+4, month = .x, data = temperature_data))
    m_4 <- map(monthly_periods, ~ extract_month(start_year = years[i]+5, month = .x, data = temperature_data))
    m_5 <- map(monthly_periods, ~ extract_month(start_year = years[i]+6, month = .x, data = temperature_data))
    
    # Application of bimonthly periods
    c_0 <- map(bi_monthly_periods, ~ extract_two_months(start_year = years[i]+1, months = .x, data = temperature_data))
    c_1 <- map(bi_monthly_periods, ~ extract_two_months(start_year = years[i]+2, months = .x, data = temperature_data))
    c_2 <- map(bi_monthly_periods, ~ extract_two_months(start_year = years[i]+3, months = .x, data = temperature_data))
    c_3 <- map(bi_monthly_periods, ~ extract_two_months(start_year = years[i]+4, months = .x, data = temperature_data))
    c_4 <- map(bi_monthly_periods, ~ extract_two_months(start_year = years[i]+5, months = .x, data = temperature_data))
    c_5 <- map(bi_monthly_periods, ~ extract_two_months(start_year = years[i]+6, months = .x, data = temperature_data))
    
    # Application of frazil and anchor ice (-10C) periods
    minus_0 <- map(minus10_periods, ~ extract_month(start_year = years[i]+1, month = .x, data = temperature_data))
    minus_1 <- map(minus10_periods, ~ extract_month(start_year = years[i]+2, month = .x, data = temperature_data))
    minus_2 <- map(minus10_periods, ~ extract_month(start_year = years[i]+3, month = .x, data = temperature_data))
    minus_3 <- map(minus10_periods, ~ extract_month(start_year = years[i]+4, month = .x, data = temperature_data))
    minus_4 <- map(minus10_periods, ~ extract_month(start_year = years[i]+5, month = .x, data = temperature_data))
    minus_5 <- map(minus10_periods, ~ extract_month(start_year = years[i]+6, month = .x, data = temperature_data))
    
    # Fonction to perform the calculation (to avoid repetition)
    calc_FDD <- function(x0, x1, x2, x3, x4, x5, threshold = 0){
      sum(
        c(
          1      * abs(x0[x0 < threshold]),
          1      * abs(x1[x1 < threshold]),
          max(0, (1-p2)) * abs(x2[x2 < threshold]),
          max(0, (1-p2-p3)) * abs(x3[x3 < threshold]),
          max(0, (1-p2-p3-p4)) * abs(x4[x4 < threshold]),
          max(0, (1-p2-p3-p4-p5)) * abs(x5[x5 < threshold])
        )
      )
    }
    
    # Application of periods between October and later months in winter
    for(nm in names(periods)){
      smolt_data[i,paste0("co_fdd_per_", nm)] <- calc_FDD(
        t_0[[nm]], t_1[[nm]], t_2[[nm]], t_3[[nm]], t_4[[nm]], t_5[[nm]], threshold = 0
      )
    }
    
    # Application of monthly periods
    for(nm in names(monthly_periods)){
      smolt_data[i,paste0("co_fdd_mon_", nm)] <- calc_FDD(
        m_0[[nm]], m_1[[nm]], m_2[[nm]], m_3[[nm]], m_4[[nm]], m_5[[nm]], threshold = 0
      )
    }
    
    # Application of bimonthly periods
    for(nm in names(bi_monthly_periods)){
      smolt_data[i,paste0("co_fdd_bimon_", nm)] <- calc_FDD(
        c_0[[nm]], c_1[[nm]], c_2[[nm]], c_3[[nm]], c_4[[nm]], c_5[[nm]], threshold = 0
      )
    }
    
    # Application of frazil and anchor ice (-10C) periods
    for(nm in names(minus10_periods)){
      smolt_data[i,paste0("co_d_minus_10_", nm)] <- calc_FDD(
        minus_0[[nm]], minus_1[[nm]], minus_2[[nm]], minus_3[[nm]], minus_4[[nm]], minus_5[[nm]], threshold = 10
      )
    }
  }
  
  # Remove years with just NA in smolt proportion, since I replace NA by 0 in the calculation, these years are super biased. (Changing NA for 0 is considered ok only for 2yrs and 5yrs because they are less individuals in these categories)
  smolt_data <- smolt_data %>%
    dplyr::filter(!(is.na(smolt_prop_2yr) &
               is.na(smolt_prop_3yr) &
               is.na(smolt_prop_4yr) &
               is.na(smolt_prop_5yr))) %>%
    dplyr::select(cohort, starts_with("co_fdd")) %>%
    as.data.frame()
  
  return(smolt_data)
}

tr_dd <- calc_dd(smolt_data = prop_smolt_sub_tr, temperature_data = tr_combined, temp_col = "Mean_T")
sj_dd <- calc_dd(smolt_data = prop_smolt_sub_sj, temperature_data = sj_combined, temp_col = "Mean_T")

# Plot Trinite data FDD
tr_dd_long <- tr_dd %>% pivot_longer(cols = -cohort, names_to = "Month", values_to = "FDD")
tr_p <- ggplot(tr_dd_long, aes(x = cohort, y = FDD, color = Month)) +
  geom_line(linewidth = 1) +         
  labs(y = "FDD", x = "Cohort", title = "Trinite") + 
  theme_minimal() +             
  theme(
    legend.title = element_blank(), 
    text = element_text(size = 12)  
  )
tr_p

# Plot St-Jean data FDD
sj_dd_long <- sj_dd %>% pivot_longer(cols = -cohort, names_to = "Month", values_to = "FDD")
sj_p <- ggplot(sj_dd_long, aes(x = cohort, y = FDD, color = Month)) +
  geom_line(linewidth = 1) +         
  labs(y = "FDD", x = "Cohort", title = "St-Jean") + 
  theme_minimal() +            
  theme(
    legend.title = element_blank(), 
    text = element_text(size = 12)  
  )
sj_p

ggarrange(tr_p, sj_p, ncol = 2, common.legend = TRUE)


# |===============================================|
# |      EXPORT COMBINED TEMPERATURE RESULTS      |
# |===============================================|

# Combine results
tr_dd$river <- "tri"
sj_dd$river <- "stj"

data_all <- bind_rows(tr_dd, sj_dd) %>%
  dplyr::select(river, cohort, everything())

write_fst(data_all, path = "03_Clean_Data/River/River_Winter_Temperatures.fst")



###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################






