




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
# |     Load libraries     |
# |========================|

library(readxl)
library(lubridate)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(broom)
library(fst)


# |=========================|
# |     Import raw data     |
# |=========================|

# Trinite
df_t <- as.data.frame(read_excel(
  "R/00_raw_data/salmon/Rapport St-Jean_Trinite_2024.xlsx", 
  sheet = "T2024FIG5", range = "O4:BF57"))

# St-Jean
df_sj <- as.data.frame(read_excel(
  "R/00_raw_data/salmon/Rapport St-Jean_Trinite_2024.xlsx", 
  sheet = "S2024FIG4", range = "P4:BB58"))


# |=================================================================|
# |     Create a function to rearange dataframe from both river     |
# |=================================================================|

arrange_bd <- function(data){

  # Make sure there are no spaces and accents in column names
  colnames(data) <- tolower(colnames(data))
  colnames(data) <- iconv(colnames(data), from = "UTF-8", to = "ASCII//TRANSLIT")
  colnames(data) <- gsub(pattern = " ", replacement = "_", x = colnames(data))
  
  # Since the year in column date is not the right one, create a column for month and a column for day
  data <- data %>%
    dplyr::mutate(date = as.Date(date),
           month = month(date),
           day = day(date)) %>%
    # Keep only columns "Month", "Day" and all year columns
    dplyr::select(month, day, matches("^[0-9]+$")) %>%
    pivot_longer(
      cols = -c(month, day),
      names_to = "cohort",
      values_to = "captures"
    ) %>%
    dplyr::mutate(across(everything(), as.numeric)) %>%
    arrange(cohort, month, day) %>%
    dplyr::mutate(date = as.Date(paste(cohort, month, day, sep = "-")),
           doy = yday(date)) %>%
    dplyr::select(date, cohort, month, day, doy, captures) %>%
    as.data.frame()
  
  # Calculate median date
  median_dates <- data %>%
    dplyr::group_by(cohort) %>%
    dplyr::filter(!is.na(captures)) %>%
    dplyr::summarise(
      median_date_doy = median(rep(doy, captures)),
      sd_date_doy    = sd(rep(doy, captures)),
      quantile_5_date = as.numeric(quantile(rep(doy, captures), 0.05)),
      quantile_25_date = as.numeric(quantile(rep(doy, captures))[2]),
      quantile_75_date = as.numeric(quantile(rep(doy, captures))[4]),
      quantile_95_date = as.numeric(quantile(rep(doy, captures), 0.95)),
      .groups = "drop"
    ) %>%
    complete(cohort = full_seq(range(cohort), 1)) %>%
    dplyr::mutate(
      median_date = as.Date(median_date_doy, origin = paste0(cohort, "-01-01")),
      quantile_5_date = as.Date(quantile_5_date, origin = paste0(cohort, "-01-01")),
      quantile_25_date = as.Date(quantile_25_date, origin = paste0(cohort, "-01-01")),
      quantile_75_date = as.Date(quantile_75_date, origin = paste0(cohort, "-01-01")),
      quantile_95_date = as.Date(quantile_95_date, origin = paste0(cohort, "-01-01"))
    ) %>%
    dplyr::select(cohort, median_date_doy, sd_date_doy, median_date, quantile_5_date, quantile_25_date, quantile_75_date, quantile_95_date) %>%
    as.data.frame()
  
  return(list(data = data, median_dates = median_dates))
}



df_t_arranged <- arrange_bd(data = df_t)
df_sj_arranged <- arrange_bd(data = df_sj)

# Add river name
df_t_arranged$median_dates <- df_t_arranged$median_dates %>%
  dplyr::mutate(river = "tri") %>%
  dplyr::select(river, everything())
df_sj_arranged$median_dates <- df_sj_arranged$median_dates %>%
  dplyr::mutate(river = "stj") %>%
  dplyr::select(river, everything())


# |=================================================================|
# |      IMPUTATION OF MISSING YEARS (±5-year sliding window)       |
# |=================================================================|

impute_missing_phenology <- function(arranged_obj) {
  
  median_dates <- arranged_obj$median_dates
  raw_data     <- arranged_obj$data
  
  # Identify missing cohorts (those with NA in median_date_doy after complete())
  missing_cohorts <- median_dates %>%
    filter(is.na(median_date_doy)) %>%
    pull(cohort)
  
  for (yr in missing_cohorts) {
    
    # Define the ±5-year window, excluding the missing year itself
    window_years <- setdiff((yr - 5):(yr + 5), yr)
    
    # Extract raw capture data for those years only
    window_data <- raw_data %>%
      filter(cohort %in% window_years, !is.na(captures)) %>%
      # Rebuild a synthetic DOY vector (repeat doy by captures)
      group_by(cohort) %>%
      reframe(doy_rep = rep(doy, captures))
    
    if (nrow(window_data) == 0) {
      message(paste("No data available in ±5-year window for cohort", yr, "— skipped."))
      next
    }
    
    # Pool all DOY values across window years
    pooled_doy <- window_data$doy_rep
    
    # Compute summary stats on pooled window
    imputed_stats <- tibble(
      cohort           = yr,
      median_date_doy  = median(pooled_doy),
      sd_date_doy      = sd(pooled_doy),
      median_date      = as.Date(median(pooled_doy), origin = paste0(yr, "-01-01")),
      quantile_5_date  = as.Date(quantile(pooled_doy, 0.05),
                                 origin = paste0(yr, "-01-01")),
      quantile_25_date = as.Date(quantile(pooled_doy, 0.25),
                                 origin = paste0(yr, "-01-01")),
      quantile_75_date = as.Date(quantile(pooled_doy, 0.75),
                                 origin = paste0(yr, "-01-01")),
      quantile_95_date = as.Date(quantile(pooled_doy, 0.95),
                                 origin = paste0(yr, "-01-01"))
    )
    
    # Inject into median_dates, preserving river column if present
    if ("river" %in% colnames(median_dates)) {
      imputed_stats <- imputed_stats %>%
        mutate(river = median_dates$river[1]) %>%
        select(river, everything())
    }
    
    median_dates <- median_dates %>%
      rows_update(imputed_stats, by = "cohort")
    
    message(paste("Cohort", yr, "imputed from window years:",
                  paste(intersect(window_years, unique(raw_data$cohort)), collapse = ", ")))
  }
  
  # Return updated object
  arranged_obj$median_dates <- median_dates
  return(arranged_obj)
}

# Apply to both rivers
df_sj_arranged <- impute_missing_phenology(df_sj_arranged)
df_t_arranged  <- impute_missing_phenology(df_t_arranged)


# |=======================|
# |     Export result     |
# |=======================|

write_fst(df_t_arranged$median_dates, path = "R/01_derived_data/Year_Median_Date_Trinite.fst")
write_fst(df_sj_arranged$median_dates, path = "R/01_derived_data/Year_Median_Date_StJean.fst")



