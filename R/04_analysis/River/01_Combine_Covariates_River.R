

# ══════════════════════════════════════════════════════════════════════════════
# HEADER
# ══════════════════════════════════════════════════════════════════════════════
#
# Merges all cleaned predictor and response datasets by cohort year for each
# river, producing a single wide-format table per river ready for variable
# selection.
#
# Variable name conventions used across all merged datasets:
#   co         = cohort-weighted metric
#   gr         = growing season period
#   em         = fry emergence period
#   dd_u       = degree days below threshold
#   dd_o       = degree days above threshold
#   temp       = water temperature
#   disc       = discharge
#   fdd_per    = freezing degree days Oct 1 to end of each month (air temp)
#   fdd_mon    = monthly freezing degree days (air temp)
#   fdd_bimon  = bimonthly freezing degree days (air temp)
#   oc/no/de/ja/fe/ma/ap = Oct/Nov/Dec/Jan/Feb/Mar/Apr
#   d_minus_10 = number of days below -10 degrees C (air temp)
#
# Inputs:
#   - R/00_raw_data/salmon/Salmon_Monitoring_StJean_Trinite_2024.csv (CSV)
#   - R/03_clean_data/River/River_Temperatures.fst (FST)
#   - R/03_clean_data/River/River_Discharge.fst (FST)
#   - R/03_clean_data/River/ForestActivity_Cohort.fst (FST)
#   - R/03_clean_data/River/River_Winter_Temperatures.fst (FST)
#
# Outputs:
#   - R/04_analysis/River/Final_Data/PreVariablesSelection_River_Trinite.fst (FST)
#   - R/04_analysis/River/Final_Data/PreVariablesSelection_River_StJean.fst (FST)
#
# ══════════════════════════════════════════════════════════════════════════════


# ══════════════════════════════════════════════════════════════════════════════
# LIBRARIES
# ══════════════════════════════════════════════════════════════════════════════

library(fst)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(purrr)


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

# ── Egg-to-smolt survival rates ───────────────────────────────────────────────
# SOURCE: MELCCFP 2024 annual report (Tables 15-16)
# year = egg-deposition year of the cohort
surv <- read.csv(
  file.path(base_path,
            "R/00_raw_data/salmon/Salmon_Monitoring_StJean_Trinite_2024.csv")) %>%
  dplyr::select(river, year, surv_egg_smolt = surv_egg_smolt_pct)

# ── Cohort-weighted river temperature metrics ─────────────────────────────────
# SOURCE: CEQEAU hydrological model outputs
riv_temp <- read_fst(
  file.path(base_path, "R/03_clean_data/River/River_Temperatures.fst"))


# ── Cohort-weighted discharge metrics ────────────────────────────────────────
# SOURCE: CEQEAU hydrological model outputs
disch <- read_fst(
  file.path(base_path, "R/03_clean_data/River/River_Discharge.fst"))


# ── Cohort-weighted forest disturbance area ───────────────────────────────────
# SOURCE: Quebec ecoforestry map (MFFP)
forest <- read_fst(
  file.path(base_path, "R/03_clean_data/River/ForestActivity_Cohort.fst"))


# ── Cohort-weighted winter freezing degree day indices ────────────────────────
# SOURCE: Combined climate station and ERA5 raster air temperature
riv_winter_cond <- read_fst(
  file.path(base_path,
            "R/03_clean_data/River/River_Winter_Temperatures.fst"))


# ══════════════════════════════════════════════════════════════════════════════
# MAIN SCRIPT
# ══════════════════════════════════════════════════════════════════════════════

# ── Exploratory plot of egg-to-smolt survival ────────────────────────────────

# Rename year to cohort before plotting so the axis label is consistent
surv <- surv %>% dplyr::rename(cohort = year)

ggplot(data = surv, aes(x = cohort, y = surv_egg_smolt, color = river)) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.8) +
  labs(x = "Year", y = "Survival rate (%)", color = "River") +
  theme_minimal(base_size = 13)


# ── Standardize column names across all datasets ─────────────────────────────

# Ensure all column names are lowercase for consistent joining
colnames(surv)            <- tolower(colnames(surv))
colnames(riv_temp)        <- tolower(colnames(riv_temp))
colnames(disch)           <- tolower(colnames(disch))
colnames(riv_winter_cond) <- tolower(colnames(riv_winter_cond))
colnames(forest)          <- tolower(colnames(forest))

# Rename any remaining "year" column to "cohort" to align join keys
if ("year" %in% colnames(surv)) {
  colnames(surv)[which(colnames(surv) == "year")] <- "cohort"
}
if ("year" %in% colnames(riv_temp)) {
  colnames(riv_temp)[which(colnames(riv_temp) == "year")] <- "cohort"
}
if ("year" %in% colnames(disch)) {
  colnames(disch)[which(colnames(disch) == "year")] <- "cohort"
}
if ("year" %in% colnames(riv_winter_cond)) {
  colnames(riv_winter_cond)[
    which(colnames(riv_winter_cond) == "year")] <- "cohort"
}


# ── Filter by river and merge all datasets ────────────────────────────────────

# Build per-river table lists; river labels may appear as full names or codes
# depending on the source script, so both variants are included in the filter
river_tables_stj <- list(
  surv            %>% filter(river %in% c("Saint-Jean", "stj")),
  riv_temp        %>% filter(river %in% c("Saint-Jean", "stj")),
  disch           %>% filter(river %in% c("Saint-Jean", "stj")),
  riv_winter_cond %>% filter(river %in% c("Saint-Jean", "stj")),
  forest          %>% filter(river %in% c("Saint-Jean", "stj"))
)

river_tables_tri <- list(
  surv            %>% filter(river %in% c("Trinite", "tri")),
  riv_temp        %>% filter(river %in% c("Trinite", "tri")),
  disch           %>% filter(river %in% c("Trinite", "tri")),
  riv_winter_cond %>% filter(river %in% c("Trinite", "tri")),
  forest          %>% filter(river %in% c("Trinite", "tri"))
)

# Full join on cohort to retain all years present in any dataset,
# then drop the redundant river columns before adding a clean river label
river_merged_stj <- reduce(river_tables_stj, full_join, by = "cohort") %>%
  dplyr::select(!dplyr::contains("river")) %>%
  dplyr::mutate(river = "stj") %>%
  dplyr::select(river, everything())

river_merged_tri <- reduce(river_tables_tri, full_join, by = "cohort") %>%
  dplyr::select(!dplyr::contains("river")) %>%
  dplyr::mutate(river = "tri") %>%
  dplyr::select(river, everything())


# ── Export ────────────────────────────────────────────────────────────────────

write_fst(river_merged_tri,
          file.path(base_path,
                    "R/04_analysis/River/Final_Data",
                    "PreVariablesSelection_River_Trinite.fst"))

write_fst(river_merged_stj,
          file.path(base_path,
                    "R/04_analysis/River/Final_Data",
                    "PreVariablesSelection_River_StJean.fst"))

