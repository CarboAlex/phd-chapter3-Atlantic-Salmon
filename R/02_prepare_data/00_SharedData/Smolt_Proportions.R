

# ══════════════════════════════════════════════════════════════════════════════
# HEADER
# ══════════════════════════════════════════════════════════════════════════════

# Computes cohort-level smolt age-class proportions for two rivers (stj, tri).
# Missing age classes are imputed with the river mean when their mean
# contribution is below a user-defined threshold; otherwise the cohort
# is excluded. Final proportions are renormalized to sum to 1.
#
# Inputs:
#   - R/00_raw_data/salmon/Salmon_Monitoring_StJean_Trinite_2024.csv
#       Smolt counts by age class, river and year (MELCCFP 2024 report, Table 5)
#
# Outputs:
#   - R/01_derived_data/Smolt_Proportions.fst  [fst data frame]
#     Columns: river, cohort, smolt_prop_2yr, smolt_prop_3yr,
#              smolt_prop_4yr, smolt_prop_5yr


# ══════════════════════════════════════════════════════════════════════════════
# LIBRARIES
# ══════════════════════════════════════════════════════════════════════════════

library(dplyr)
library(tidyr)
library(fst)


# ══════════════════════════════════════════════════════════════════════════════
# WORKING DIRECTORY
# ══════════════════════════════════════════════════════════════════════════════

# ── Working directory ─────────────────────────────────────────────────────────
# Automatically detects the repository root based on the script location.
# No hardcoded path — works on any machine.
current_file <- rstudioapi::getActiveDocumentContext()$path
path_parts   <- unlist(strsplit(normalizePath(current_file, winslash = "/"), "/"))
target_index <- which(path_parts == "phd-chapter3-Atlantic-Salmon")
base_path    <- paste(path_parts[1:target_index], collapse = "/")
setwd(base_path)


# ══════════════════════════════════════════════════════════════════════════════
# DATA IMPORT
# ══════════════════════════════════════════════════════════════════════════════

# ── Smolt counts by smoltification year ───────────────────────────────────────
# SOURCE: MELCCFP annual report 2024 (Table 5)
prop_smolt <- read.csv(
  file.path(base_path,
            "R/00_raw_data/salmon/Salmon_Monitoring_StJean_Trinite_2024.csv")
)


# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

# Age classes whose mean proportional contribution is below this threshold
# are eligible for mean imputation when their smoltification-year count is
# missing. Cohorts with any ineligible NA are excluded entirely.
NA_THRESHOLD <- 0.10   # 10 % — adjust here for sensitivity analysis


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

# ── build_cohort_props ────────────────────────────────────────────────────────
# Purpose: Aligns smolt raw counts to the egg-deposition cohort year by
#          applying lead() offsets matching each age class's smoltification
#          schedule (2yr -> +3 yrs, 3yr -> +4 yrs, 4yr -> +5 yrs,
#          5yr -> +6 yrs).
# Parameters:
#   river_data  data frame for a single river containing columns year,
#               smolt_nb_2yr, smolt_nb_3yr, smolt_nb_4yr, smolt_nb_5yr.
build_cohort_props <- function(river_data) {
  river_data %>%
    dplyr::arrange(year) %>%
    dplyr::mutate(
      cohort         = year,
      smolt_nb_2yr_c = dplyr::lead(smolt_nb_2yr, n = 3),
      smolt_nb_3yr_c = dplyr::lead(smolt_nb_3yr, n = 4),
      smolt_nb_4yr_c = dplyr::lead(smolt_nb_4yr, n = 5),
      smolt_nb_5yr_c = dplyr::lead(smolt_nb_5yr, n = 6)
    ) %>%
    dplyr::select(cohort,
                  smolt_nb_2yr_c, smolt_nb_3yr_c,
                  smolt_nb_4yr_c, smolt_nb_5yr_c)
}


# ── impute_and_normalize ──────────────────────────────────────────────────────
# Purpose: Applies the imputation and exclusion rules to cohort-level raw
#          counts, then computes and renormalizes smolt age-class proportions.
# Parameters:
#   cohort_raw    data frame produced by build_cohort_props() for one river.
#   river_id      character string identifying the river ("stj" or "tri").
#   ref_stats     data frame of per-river summary statistics from
#                 prop_by_devalaison_year (used to assess eligibility).
#   ref_means_raw data frame of per-river mean raw counts over complete years
#                 (used as imputation values).
#   threshold     numeric; NA_THRESHOLD defined in Section 5.
# Returns: named list with two elements:
#   $proportions    data frame: cohort, smolt_prop_2yr–5yr.
#   $imputation_log data frame recording every imputed cell.
impute_and_normalize <- function(cohort_raw, river_id,
                                 ref_stats, ref_means_raw, threshold) {
  
  stats_r   <- ref_stats     %>% dplyr::filter(river == river_id)
  means_raw <- ref_means_raw %>% dplyr::filter(river == river_id)
  
  # Determine eligibility: age classes whose mean proportion is below threshold
  # are imputed with the river mean when their count is missing.
  eligible <- list(
    `2yr` = stats_r$mean_2yr < threshold,
    `3yr` = stats_r$mean_3yr < threshold,
    `4yr` = stats_r$mean_4yr < threshold,
    `5yr` = stats_r$mean_5yr < threshold
  )
  
  # Initialise the imputation log
  imputation_log <- data.frame(
    river  = character(),
    cohort = integer(),
    age    = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
  
  result <- cohort_raw %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      
      # Flag cohorts where every age class is missing (no smoltification data)
      all_na = all(is.na(c(smolt_nb_2yr_c, smolt_nb_3yr_c,
                           smolt_nb_4yr_c, smolt_nb_5yr_c))),
      
      # Impute raw counts where eligible; leave as NA otherwise so the cohort
      # is subsequently excluded (see `exclude` below)
      nb2 = dplyr::case_when(
        all_na                  ~ NA_real_,
        !is.na(smolt_nb_2yr_c)  ~ smolt_nb_2yr_c,
        eligible[["2yr"]]       ~ means_raw$mean_raw_2yr,
        TRUE                    ~ NA_real_
      ),
      nb3 = dplyr::case_when(
        all_na                  ~ NA_real_,
        !is.na(smolt_nb_3yr_c)  ~ smolt_nb_3yr_c,
        eligible[["3yr"]]       ~ means_raw$mean_raw_3yr,
        TRUE                    ~ NA_real_
      ),
      nb4 = dplyr::case_when(
        all_na                  ~ NA_real_,
        !is.na(smolt_nb_4yr_c)  ~ smolt_nb_4yr_c,
        eligible[["4yr"]]       ~ means_raw$mean_raw_4yr,
        TRUE                    ~ NA_real_
      ),
      nb5 = dplyr::case_when(
        all_na                  ~ NA_real_,
        !is.na(smolt_nb_5yr_c)  ~ smolt_nb_5yr_c,
        eligible[["5yr"]]       ~ means_raw$mean_raw_5yr,
        TRUE                    ~ NA_real_
      ),
      
      # Exclude cohorts where at least one ineligible age class remains NA
      exclude = !all_na & any(is.na(c(nb2, nb3, nb4, nb5))),
      
      # Compute proportions; set to NA for excluded or fully missing cohorts
      total          = ifelse(all_na | exclude, NA_real_,
                              nb2 + nb3 + nb4 + nb5),
      smolt_prop_2yr = ifelse(all_na | exclude, NA_real_, nb2 / total),
      smolt_prop_3yr = ifelse(all_na | exclude, NA_real_, nb3 / total),
      smolt_prop_4yr = ifelse(all_na | exclude, NA_real_, nb4 / total),
      smolt_prop_5yr = ifelse(all_na | exclude, NA_real_, nb5 / total)
    ) %>%
    dplyr::ungroup()
  
  # Build the imputation log by checking which eligible NAs were replaced
  for (age in c("2yr", "3yr", "4yr", "5yr")) {
    raw_col       <- paste0("smolt_nb_", age, "_c")
    eligible_flag <- eligible[[age]]
    if (eligible_flag) {
      imputed_rows <- result %>%
        dplyr::filter(is.na(.data[[raw_col]]) & !all_na & !exclude)
      if (nrow(imputed_rows) > 0) {
        imputation_log <- rbind(
          imputation_log,
          data.frame(
            river  = river_id,
            cohort = imputed_rows$cohort,
            age    = age,
            reason = paste0(
              "NA replaced by river mean (mean = ",
              round(eligible_flag * 100, 1),
              "% < ", threshold * 100, "% threshold)"
            ),
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
  
  list(
    proportions    = result %>%
      dplyr::select(cohort, smolt_prop_2yr, smolt_prop_3yr,
                    smolt_prop_4yr, smolt_prop_5yr) %>%
      as.data.frame(),
    imputation_log = imputation_log,
    eligible       = eligible
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# MAIN SCRIPT
# ══════════════════════════════════════════════════════════════════════════════

# ── Select and rename smolt count columns ────────────────────────────────────
prop_smolt_sub <- prop_smolt %>%
  dplyr::select(river, year,
                smolt_n_age2, smolt_n_age3, smolt_n_age4, smolt_n_age5) %>%
  dplyr::rename(
    smolt_nb_2yr = smolt_n_age2,
    smolt_nb_3yr = smolt_n_age3,
    smolt_nb_4yr = smolt_n_age4,
    smolt_nb_5yr = smolt_n_age5
  )

# Ensure every year between the observed minimum and maximum is represented;
# gaps are filled with NA so the lead() offsets in build_cohort_props() are
# applied to the correct calendar positions.
prop_smolt_sub <- prop_smolt_sub %>%
  dplyr::group_by(river) %>%
  tidyr::complete(year = seq(min(year), max(year)), fill = list(
    smolt_nb_2yr = NA,
    smolt_nb_3yr = NA,
    smolt_nb_4yr = NA,
    smolt_nb_5yr = NA
  )) %>%
  dplyr::ungroup()


# ── Compute raw proportions by smoltification year ───────────────────────────
# Proportions here are per smoltification year (not per cohort). They serve
# only to derive reference means and ranges used in the methods text, and to
# assess eligibility for imputation. Only years where all four age classes
# were sampled (complete years) are included.
prop_by_devalaison_year <- prop_smolt_sub %>%
  dplyr::filter(
    !is.na(smolt_nb_2yr) & !is.na(smolt_nb_3yr) &
      !is.na(smolt_nb_4yr) & !is.na(smolt_nb_5yr)
  ) %>%
  dplyr::mutate(
    total    = smolt_nb_2yr + smolt_nb_3yr + smolt_nb_4yr + smolt_nb_5yr,
    prop_2yr = smolt_nb_2yr / total,
    prop_3yr = smolt_nb_3yr / total,
    prop_4yr = smolt_nb_4yr / total,
    prop_5yr = smolt_nb_5yr / total
  )

# Summary statistics by river — used to decide eligibility for imputation
# (mean < NA_THRESHOLD) and to report variability in the methods text.
ref_stats <- prop_by_devalaison_year %>%
  dplyr::group_by(river) %>%
  dplyr::summarise(
    n_complete_years = n(),
    mean_2yr = mean(prop_2yr), sd_2yr = sd(prop_2yr),
    cv_2yr   = sd(prop_2yr) / mean(prop_2yr),
    min_2yr  = min(prop_2yr),  max_2yr = max(prop_2yr),
    mean_3yr = mean(prop_3yr), sd_3yr = sd(prop_3yr),
    cv_3yr   = sd(prop_3yr) / mean(prop_3yr),
    min_3yr  = min(prop_3yr),  max_3yr = max(prop_3yr),
    mean_4yr = mean(prop_4yr), sd_4yr = sd(prop_4yr),
    cv_4yr   = sd(prop_4yr) / mean(prop_4yr),
    min_4yr  = min(prop_4yr),  max_4yr = max(prop_4yr),
    mean_5yr = mean(prop_5yr), sd_5yr = sd(prop_5yr),
    cv_5yr   = sd(prop_5yr) / mean(prop_5yr),
    min_5yr  = min(prop_5yr),  max_5yr = max(prop_5yr),
    .groups = "drop"
  )


# ── Align counts to egg-deposition cohorts ───────────────────────────────────
# A missing smoltification year creates a single NA in exactly one age class
# for up to five different cohorts, never multiple NAs in the same cohort
# (unless two consecutive smoltification years are missing).
cohort_raw_sj <- build_cohort_props(
  prop_smolt_sub %>% dplyr::filter(river == "stj")
)
cohort_raw_tr <- build_cohort_props(
  prop_smolt_sub %>% dplyr::filter(river == "tri")
)


# ── Compute per-river mean raw counts for imputation ─────────────────────────
# Imputation is performed on the raw count scale (before dividing by total)
# to avoid distorting the denominator.
ref_means_raw <- prop_by_devalaison_year %>%
  dplyr::group_by(river) %>%
  dplyr::summarise(
    mean_raw_2yr = mean(smolt_nb_2yr),
    mean_raw_3yr = mean(smolt_nb_3yr),
    mean_raw_4yr = mean(smolt_nb_4yr),
    mean_raw_5yr = mean(smolt_nb_5yr),
    .groups = "drop"
  )


# ── Apply imputation and compute cohort proportions ──────────────────────────
result_sj <- impute_and_normalize(
  cohort_raw_sj, "stj", ref_stats, ref_means_raw, NA_THRESHOLD
)
result_tr <- impute_and_normalize(
  cohort_raw_tr, "tri", ref_stats, ref_means_raw, NA_THRESHOLD
)

prop_smolt_sub_sj <- result_sj$proportions
prop_smolt_sub_tr <- result_tr$proportions


# ── Combine rivers and export ─────────────────────────────────────────────────
prop_smolt_sub_sj$river <- "stj"
prop_smolt_sub_tr$river <- "tri"

smolt_props_all <- dplyr::bind_rows(prop_smolt_sub_sj, prop_smolt_sub_tr) %>%
  dplyr::select(river, cohort, smolt_prop_2yr, smolt_prop_3yr,
                smolt_prop_4yr, smolt_prop_5yr)

write_fst(
  smolt_props_all,
  file.path(base_path, "R/01_derived_data/Smolt_Proportions.fst")
)








