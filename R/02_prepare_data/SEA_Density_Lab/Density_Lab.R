


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


# |========================|
# |     LOAD LIBRARIES     |
# |========================|

library(tidyverse)
library(readxl)
library(ggplot2)

# ---- Output path ----
out_path   <- file.path(base_path, "R/03_clean_data/Sea")


# ---- Data import ----

# NAC (North America) — already at January 1st
df_NAC <- as.data.frame(
  read_excel(
    file.path(base_path, "R/00_raw_data/Exploitation_Sea/WGNAS_2025_Table4_3_6_1_EstimatesPFA_Jan1st_NAC.xlsx")
  )
) %>%
  dplyr::select(
    YEAR      = Year,
    NACcohort = Est_PFA1SWcohort,  # 1st year at sea (mat + nmat)
    NACnonmat = Est_PFA1SWnmat     # 2nd year at sea (nmat only)
  )

# SNEAC (Southern Europe) — already at January 1st
df_SNEAC <- as.data.frame(
  read_excel(
    file.path(base_path, "R/00_raw_data/Exploitation_Sea/WGNAS_2025_Table3_3_5_3_EstimatesPFA_Jan1st_SNEAC.xlsx")
  )
) %>%
  dplyr::select(YEAR = Year, SNEACnonmat = SNEAC_median)

# ---- Processing — Tara Imlay's method ----
# pfa_1sw : 1st year abundance = NAC smolt cohort (mat + nmat)
#            -> only North American salmon present in the NW Atlantic
# pfa_2sw : 2nd year abundance = NAC non-matures + SNEAC non-matures
#            -> both stocks co-occur (Labrador, Greenland, Iceland, Faroes)

Density_Lab <- left_join(df_NAC, df_SNEAC, by = "YEAR") %>%
  mutate(
    pfa_1sw = NACcohort,
    pfa_2sw = NACnonmat + SNEACnonmat
  ) %>%
  dplyr::select(YEAR, pfa_1sw, pfa_2sw)

# ---- Export ----
write_fst(Density_Lab, path = file.path(out_path, "Density_Lab.fst"), compress = 100)

# ---- Figure ----
Density_Lab_long <- Density_Lab %>%
  pivot_longer(
    cols      = c(pfa_1sw, pfa_2sw),
    names_to  = "SW",
    values_to = "PFA"
  ) %>%
  mutate(SW = dplyr::recode(SW,
                            "pfa_1sw" = "1st year at sea (1SW)",
                            "pfa_2sw" = "2nd year at sea (2SW)"
  ))

ggplot(Density_Lab_long, aes(x = YEAR, y = PFA, color = SW)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("1st year at sea (1SW)" = "#2196F3",
                                "2nd year at sea (2SW)" = "#E53935")) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title    = "Pre-fishery abundance (PFA) of Atlantic salmon",
    subtitle = "Separated by year at sea — Tara's method",
    x        = "Year",
    y        = "PFA (thousands of individuals)",
    color    = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position   = c(0.75, 0.9),
    legend.background = element_rect(fill = "white", color = NA),
    plot.title        = element_text(face = "bold"),
    plot.subtitle     = element_text(color = "grey40")
  )



