

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: HEADER
# ══════════════════════════════════════════════════════════════════════════════
#
# Imports and cleans salmon demographic rates from the MELCCFP BDFUSION
# database, produces exploratory plots, and exports egg-to-smolt survival and
# smolt-to-adult return rates for use in downstream analyses.
#
# Note: the Trinite River dam was built in 1977 and modified in 2010.
#
# Inputs:
#   - R/00_raw_data/salmon/BDFUSION_2024.xlsm
#       MELCCFP salmon monitoring database (Excel)
#
# Outputs:
#   - R/03_clean_data/River/SalmonSurvival_River.fst
#       Annual egg-to-smolt survival rate by river and cohort year (FST)
#   - R/03_clean_data/Sea/SalmonSurvival_Sea.fst
#       Annual smolt-to-adult return rates (1SW and 2SW) by river and cohort
#       year (FST)
#
# ══════════════════════════════════════════════════════════════════════════════


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: LIBRARIES
# ══════════════════════════════════════════════════════════════════════════════

library(readxl)
library(dplyr)
library(ggplot2)
library(fst)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: WORKING DIRECTORY
# ══════════════════════════════════════════════════════════════════════════════

# ── Working directory ──────────────────────────────────────────────────────
# Automatically detects the repository root based on the script location.
# No hardcoded path — works on any machine.
current_file <- rstudioapi::getActiveDocumentContext()$path
path_parts   <- unlist(strsplit(normalizePath(current_file, winslash = "/"), "/"))
target_index <- which(path_parts == "phd-chapter3-Atlantic-Salmon")
base_path    <- paste(path_parts[1:target_index], collapse = "/")
setwd(base_path)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: DATA IMPORT
# ══════════════════════════════════════════════════════════════════════════════

# ── MELCCFP salmon monitoring database ────────────────────────────────────────
# SOURCE: Quebec Ministry of Environment (MELCCFP), BDFUSION database
# STATUS: Confidential — not included in this repository.
#         An example file with identical structure and simulated data is
#         available at the same path with the suffix _EXAMPLE.xlsm.
#         To request access, contact the MELCCFP once the thesis is published.
df <- as.data.frame(
  read_excel(file.path(base_path, "R/00_raw_data/salmon/BDFUSION_2024.xlsm"),
             sheet = 1, n_max = 91)
  )


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5: CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

dam_year <- 2010  # Year of the Trinite River dam modification


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 7: MAIN SCRIPT
# ══════════════════════════════════════════════════════════════════════════════

# ── Clean column names ────────────────────────────────────────────────────────

# Strip accents, punctuation, and special characters from column headers
# so that all names are ASCII-safe and usable as R identifiers
colnames(df) <- iconv(colnames(df), from = "UTF-8", to = "ASCII//TRANSLIT")
colnames(df) <- gsub("\\.",     "",       colnames(df))
colnames(df) <- gsub("\\(",     "",       colnames(df))
colnames(df) <- gsub("\\)",     "",       colnames(df))
colnames(df) <- gsub(" ",       "_",      colnames(df))
colnames(df) <- gsub("\\%",     "PC",     colnames(df))
colnames(df) <- gsub("\\+",     "P",      colnames(df))
colnames(df) <- gsub("/",       "_Par_",  colnames(df))
colnames(df) <- gsub("-",       "_",      colnames(df))
colnames(df) <- gsub("l'annee", "annee",  colnames(df))
colnames(df) <- gsub("oufs",    "oeufs",  colnames(df))
colnames(df) <- gsub("__",      "_",      colnames(df))

df$Riviere <- iconv(df$Riviere, from = "UTF-8", to = "ASCII//TRANSLIT")


# ── Subset and rescale rate columns ──────────────────────────────────────────

# Keep only the columns needed for survival and return rate analyses
df <- subset(df, select = c(
  "Riviere", "Annee",
  "Taux_survie_PC_Oeuf_Smolt",
  "Taux_de_retour_PC_Smolt_Adulte_Madeleinaux",
  "Taux_de_retour_PC_Smolt_Adulte_Dibermarins",
  "Taux_de_retour_PC_Smolt_Adulte_Tribermarins",
  "Taux_de_retour_PC_Smolt_Adulte_Tous"
))

# Convert proportions to percentages for readability in plots
df$Taux_survie_PC_Oeuf_Smolt                  <-
  df$Taux_survie_PC_Oeuf_Smolt                  * 100
df$Taux_de_retour_PC_Smolt_Adulte_Tous        <-
  df$Taux_de_retour_PC_Smolt_Adulte_Tous        * 100
df$Taux_de_retour_PC_Smolt_Adulte_Madeleinaux <-
  df$Taux_de_retour_PC_Smolt_Adulte_Madeleinaux * 100
df$Taux_de_retour_PC_Smolt_Adulte_Dibermarins <-
  df$Taux_de_retour_PC_Smolt_Adulte_Dibermarins * 100


# ── Exploratory plots ─────────────────────────────────────────────────────────

# Shared theme applied to all plots
theme_surv <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggplot(data = df,
       aes(x = Annee, y = Taux_survie_PC_Oeuf_Smolt,
           color = factor(Riviere))) +
  geom_point(size = 1.5) +
  geom_line() +
  labs(color = "River", x = "Cohort",
       y = "Survival rate egg to smolt (%)") +
  theme_surv

ggplot(data = df,
       aes(x = Annee, y = Taux_de_retour_PC_Smolt_Adulte_Tous,
           color = factor(Riviere))) +
  geom_point(size = 1.5) +
  geom_line() +
  geom_vline(xintercept = dam_year, color = "grey30",
             linetype = "dashed", linewidth = 0.8) +
  geom_text(aes(x = dam_year, y = 5,
                label = " Modification of the\n Trinite R. dam"),
            color = "grey30", hjust = 0, inherit.aes = FALSE) +
  labs(color = "River", x = "Cohort",
       y = "Return rate smolt to adult (%)") +
  theme_surv

ggplot(data = df,
       aes(x = Annee, y = Taux_de_retour_PC_Smolt_Adulte_Madeleinaux,
           color = factor(Riviere))) +
  geom_point(size = 1.5) +
  geom_line() +
  geom_vline(xintercept = dam_year, color = "grey30",
             linetype = "dashed", linewidth = 0.8) +
  geom_text(aes(x = dam_year, y = 3.5,
                label = " Modification of the\n Trinite R. dam"),
            color = "grey30", hjust = 0, inherit.aes = FALSE) +
  labs(color = "River", x = "Cohort",
       y = "Return rate smolt to 1SW (%)") +
  theme_surv

ggplot(data = df,
       aes(x = Annee, y = Taux_de_retour_PC_Smolt_Adulte_Dibermarins,
           color = factor(Riviere))) +
  geom_point(size = 1.5) +
  geom_line() +
  geom_vline(xintercept = dam_year, color = "grey30",
             linetype = "dashed", linewidth = 0.8) +
  geom_text(aes(x = dam_year, y = 3,
                label = " Modification of the\n Trinite R. dam"),
            color = "grey30", hjust = 0, inherit.aes = FALSE) +
  labs(color = "River", x = "Cohort",
       y = "Return rate smolt to 2SW (%)") +
  theme_surv


# ── Rename columns and split by life stage ────────────────────────────────────

colnames(df) <- tolower(colnames(df))

df <- df %>%
  dplyr::rename(
    river          = riviere,
    year           = annee,
    surv_egg_smolt = taux_survie_pc_oeuf_smolt,
    return_1sw     = taux_de_retour_pc_smolt_adulte_madeleinaux,
    return_2sw     = taux_de_retour_pc_smolt_adulte_dibermarins
  ) %>%
  dplyr::select(river, year, surv_egg_smolt, return_1sw, return_2sw)

# Separate freshwater survival from marine return rates for export
df_river <- df %>% dplyr::select(river, year, surv_egg_smolt)
df_sea   <- df %>% dplyr::select(river, year, return_1sw, return_2sw)


# ── Export ────────────────────────────────────────────────────────────────────

write_fst(df_river,
          path = file.path(base_path, "R/03_clean_data/River/SalmonSurvival_River.fst"))

write_fst(df_sea,
          path = file.path(base_path, "R/03_clean_data/Sea/SalmonSurvival_Sea.fst"))

