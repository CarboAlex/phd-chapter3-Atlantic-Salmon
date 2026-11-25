


# ==============================================================================
# Processes raw CPR (Continuous Plankton Recorder) data for the NW Atlantic:
# renames taxa columns, converts presence-only flags to NA, converts abundance
# to biomass using allometric formulas, and exports the tidy long-format table.
#
# Inputs:
#   - R/00_raw_data/Plankton_LabradorSea/NWAtl_CPR_Data.xlsx
#       Raw CPR dataset (NW Atlantic, 1958-2021). Excel workbook.
#
# Outputs:
#   - R/01_derived_data/Food_Availability/cpr_data_raw.fst
#       Long-format CPR table with per-sample biomass estimates. FST format.
# ==============================================================================


# ══════════════════════════════════════════════════════════════════════════════
# LIBRARIES
# ══════════════════════════════════════════════════════════════════════════════

library(readxl)
library(fst)
library(data.table)
library(dplyr)
library(tidyr)
library(lubridate)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(viridis)
library(units)


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

# ── CPR raw data ──────────────────────────────────────────────────────────────
# SOURCE: Continuous Plankton Recorder (CPR) Survey, NW Atlantic.
#   David Johns (Marine Biological Association) (2023): Selected NW Atlantic
#   CPR taxa. The Archive for Marine Species and Habitats Data (DASSH).
#   https://doi.org/10.17031/6513f104926b6
cpr_raw <- read_excel(
  path = file.path(base_path,
                   "R/00_raw_data/Plankton_LabradorSea/NWAtl_CPR_Data.xlsx")
)
cpr_raw <- as.data.table(cpr_raw)


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

# ── Biomass conversion: Richardson et al. (2006) ──────────────────────────────
# Purpose: Computes dry mass (g) per individual using the allometric formula
#   mass = 0.08 * length^2.1, then multiplies by abundance to get biomass.
#   Applied to all taxa except Oithona spp.
# Parameters:
#   abundance  - numeric; number of individuals per CPR sample
#   length_mm  - numeric; body length in millimetres
richardson_et_al_2006 <- function(abundance, length_mm) {
  mass    <- 0.08 * length_mm^2.1
  biomass <- abundance * mass
  return(biomass)
}

# ── Biomass conversion: Krylov (1968) ─────────────────────────────────────────
# Purpose: Computes dry mass (g) per individual using the allometric formula
#   mass = 0.008 * length^3, then multiplies by abundance to get biomass.
#   Applied exclusively to Oithona spp.
# Parameters:
#   abundance  - numeric; number of individuals per CPR sample
#   length_mm  - numeric; body length in millimetres
krylov_1968 <- function(abundance, length_mm) {
  mass    <- 0.008 * length_mm^3
  biomass <- abundance * mass
  return(biomass)
}


# ══════════════════════════════════════════════════════════════════════════════
# MAIN SCRIPT
# ══════════════════════════════════════════════════════════════════════════════

# ── Rename taxa columns ────────────────────────────────────────────────────────
# Original column names are numeric taxon codes; replace with descriptive names.
setnames(cpr_raw, old = "1",     new = "Calanus_I_IV")
setnames(cpr_raw, old = "3",     new = "Para_Pseudocalanus_spp")
setnames(cpr_raw, old = "10",    new = "Oithona")
setnames(cpr_raw, old = "40",    new = "Calanus_finmarchicus")
setnames(cpr_raw, old = "41",    new = "Calanus_helgolandicus")
setnames(cpr_raw, old = "42",    new = "Calanus_glacialis")
setnames(cpr_raw, old = "44",    new = "Calanus_hyperboreus")
setnames(cpr_raw, old = "55",    new = "Metridia_lucens")
setnames(cpr_raw, old = "56",    new = "Metridia_longa")
setnames(cpr_raw, old = "82",    new = "Hyperiidea")
setnames(cpr_raw, old = "88",    new = "Euphausiacea")
setnames(cpr_raw, old = "314",   new = "Metridia_I_IV")
setnames(cpr_raw, old = "10742", new = "Pseudocalanus_spp_ad_tot")

# ── Recode presence-only flags to NA ──────────────────────────────────────────
# CPR convention: 0.001 indicates presence with no quantitative abundance data.
# Replace with NA so downstream biomass calculations are not biased.
cpr_raw[, (7:ncol(cpr_raw)) := lapply(.SD, function(x) {
  fifelse(x == 0.001, NA_real_, x)
}), .SDcols = 7:ncol(cpr_raw)]

# ── Visualise sampling locations ───────────────────────────────────────────────
# Retain one row per unique sample for spatial display.
cpr_points <- unique(
  cpr_raw[, .(Sample_Id, Latitude, Longitude, Midpoint_Date_Local)]
)

# Convert to sf for spatial operations and plotting.
cpr_sf  <- st_as_sf(cpr_points, coords = c("Longitude", "Latitude"), crs = 4326)
na_map  <- ne_countries(scale = "medium", continent = "North America",
                        returnclass = "sf")

ggplot() +
  geom_sf(data = na_map, fill = "grey90", color = "grey50") +
  geom_sf(data = cpr_sf, color = "red", size = 0.5, alpha = 0.4) +
  coord_sf(xlim = c(-66, -15), ylim = c(44, 70), expand = FALSE) +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude")

# ── Reshape to long format ─────────────────────────────────────────────────────
# One row per sample-taxon combination to facilitate per-taxon biomass
# computation and subsequent aggregation.
taxa_cols <- setdiff(
  names(cpr_raw),
  c("Sample_Id", "Latitude", "Longitude", "Midpoint_Date_Local",
    "Year", "Month", "Chlorophyll_Index")
)

cpr_long <- melt(
  cpr_raw,
  id.vars     = c("Sample_Id", "Latitude", "Longitude", "Midpoint_Date_Local",
                  "Year", "Month", "Chlorophyll_Index"),
  measure.vars  = taxa_cols,
  variable.name = "taxa",
  value.name    = "count"
)

# ── Convert abundance to biomass ──────────────────────────────────────────────
# Body lengths follow Imlay et al. (2024) and the sources cited therein:
#   Calanus I-IV        : 1.65 mm  (Richardson et al., 2006)
#   C. finmarchicus     : 2.70 mm  (Richardson et al., 2006)
#   C. glacialis        : 4.60 mm  (Richardson et al., 2006)
#   C. hyperboreus      : 6.95 mm  (Richardson et al., 2006)
#   C. helgolandicus    : 2.68 mm  (Richardson et al., 2006)
#   Metridia I-IV       : 0.93 mm  (Richardson et al., 2006)
#   M. longa            : 4.10 mm  (Richardson et al., 2006)
#   M. lucens           : 2.27 mm  (Richardson et al., 2006)
#   Para-Pseudocalanus  : 0.70 mm  (Richardson et al., 2006)
#   Pseudocalanus ad.   : 1.20 mm  (Richardson et al., 2006)
#   Euphausiacea        : 8.67 mm  (Lindley, 1978)
#   Hyperiidea          : 3.00 mm  (Williams & Robins, 1981; Imlay et al., 2024)
#   Oithona spp.        : 0.68 mm  (Krylov formula)
zooLength <- data.frame(
  taxa = c("Calanus_I_IV", "Calanus_finmarchicus", "Calanus_glacialis",
           "Calanus_hyperboreus", "Calanus_helgolandicus",
           "Metridia_I_IV", "Metridia_longa", "Metridia_lucens",
           "Para_Pseudocalanus_spp", "Pseudocalanus_spp_ad_tot",
           "Euphausiacea", "Hyperiidea", "Oithona"),
  length_mm = c(1.65, 2.70, 4.60,
                6.95, 2.68,
                0.93, 4.10, 2.27,
                0.7,  1.20,
                8.67, 3.00, 0.68)
)

# Join lengths and apply the appropriate formula per taxon.
cpr_long <- cpr_long %>%
  left_join(zooLength, by = "taxa") %>%
  mutate(biomass_g = if_else(
    taxa == "Oithona",
    krylov_1968(count, length_mm),
    richardson_et_al_2006(count, length_mm)
  ))

# ── Aggregate biomass per sample ──────────────────────────────────────────────
# Sum abundance and biomass across all taxa within each CPR sample.
cpr_sample <- cpr_long %>%
  group_by(Sample_Id, Latitude, Longitude, Year, Month) %>%
  summarise(
    countSample    = sum(count,      na.rm = TRUE),
    biomass_gSample = sum(biomass_g, na.rm = TRUE),
    .groups = "drop"
  )

# ── Export ─────────────────────────────────────────────────────────────────────
fst::write_fst(
  cpr_long,
  path     = file.path(base_path,
                       "R/01_derived_data/Food_Availability/cpr_data_raw.fst"),
  compress = 100
)





