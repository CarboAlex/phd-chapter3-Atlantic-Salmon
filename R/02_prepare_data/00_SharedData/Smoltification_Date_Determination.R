

# ══════════════════════════════════════════════════════════════════════════════
# HEADER
# ══════════════════════════════════════════════════════════════════════════════
#
# Smolt migration dates per cohort
#
# Description:
#   Computes, for each cohort, the median, standard deviation and quantiles
#   (5, 25, 75, 95 %) of the smolt migration day of year, from the daily
#   captures at the downstream traps of the St. Jean and Trinite rivers.
#   Cohorts without captures are imputed by pooling the captures of the
#   years within a ±5-year window.
#
# Inputs:
#   R/00_raw_data/salmon/Smolt_Daily_Captures_StJean.csv
#   R/00_raw_data/salmon/Smolt_Daily_Captures_Trinite.csv
#
# Outputs:
#   R/01_derived_data/Year_Median_Date_StJean.fst
#   R/01_derived_data/Year_Median_Date_Trinite.fst


# ══════════════════════════════════════════════════════════════════════════════
# LIBRARIES
# ══════════════════════════════════════════════════════════════════════════════

library(tidyr)
library(dplyr)
library(fst)


# ══════════════════════════════════════════════════════════════════════════════
# WORKING DIRECTORY
# ══════════════════════════════════════════════════════════════════════════════

# Sets the working directory to the repository root (the folder containing the
# .Rproj file), whatever the folder is named. 
find_repo_root <- function() {
  start <- tryCatch(dirname(rstudioapi::getSourceEditorContext()$path),
                    error = function(e) "")
  if (!nzchar(start)) {
    f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    start <- if (length(f)) dirname(normalizePath(f)) else getwd()
  }
  d <- normalizePath(start, winslash = "/")
  while (length(list.files(d, pattern = "\\.Rproj$")) == 0) {
    if (dirname(d) == d) stop("Repository root not found (no .Rproj file above this script).")
    d <- dirname(d)
  }
  d
}
setwd(find_repo_root())

# ══════════════════════════════════════════════════════════════════════════════
# DATA IMPORT
# ══════════════════════════════════════════════════════════════════════════════

col_types <- c(date = "Date", cohort = "numeric", doy = "numeric", captures = "numeric")

# Trinite
df_t <- read.csv("R/00_raw_data/salmon/Smolt_Daily_Captures_Trinite.csv",
                 colClasses = col_types)

# St-Jean
df_sj <- read.csv("R/00_raw_data/salmon/Smolt_Daily_Captures_StJean.csv",
                  colClasses = col_types)


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

# Function computing the migration dates of each cohort
arrange_bd <- function(data){
  
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

# Imputation of missing years (±5-year sliding window)
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


# ══════════════════════════════════════════════════════════════════════════════
# MIGRATION DATES
# ══════════════════════════════════════════════════════════════════════════════

df_t_arranged <- arrange_bd(data = df_t)
df_sj_arranged <- arrange_bd(data = df_sj)

# Add river name
df_t_arranged$median_dates <- df_t_arranged$median_dates %>%
  dplyr::mutate(river = "tri") %>%
  dplyr::select(river, everything())
df_sj_arranged$median_dates <- df_sj_arranged$median_dates %>%
  dplyr::mutate(river = "stj") %>%
  dplyr::select(river, everything())

# Apply imputation to both rivers
df_sj_arranged <- impute_missing_phenology(df_sj_arranged)
df_t_arranged  <- impute_missing_phenology(df_t_arranged)


# ══════════════════════════════════════════════════════════════════════════════
# EXPORT
# ══════════════════════════════════════════════════════════════════════════════

write_fst(df_t_arranged$median_dates, path = "R/01_derived_data/Year_Median_Date_Trinite.fst")
write_fst(df_sj_arranged$median_dates, path = "R/01_derived_data/Year_Median_Date_StJean.fst")

