

# NE PAS OUBLIER DE MODIFIER LES DATES AUXQUELLES JE PRENDS LES DONNEES DE DEBIT
# EN FONCTION DES PERIODES D'EMERGENCE ET DE GROWTH SI POSSIBLE. 
# SI IMPOSSIBLE, ALLER AVEC LES LIMITES SUGGEREES PAR ILIAS. 
# DANS TOUS LES CAS, LES CHANGEMENTS DE DATES DOIVENT SE FAIRE DIRECT DANS LES DATAFRAMES "PERIOD..."



# |======================================================|
# |    SET THE WORKING DIRECTORY SO IT'S "Chapitre_3"    |
# |======================================================|

current_file <- rstudioapi::getActiveDocumentContext()$path # Get the location of this file (not the working directory of the open R session)
path <- normalizePath(current_file, winslash = "/", mustWork = TRUE) # Get the complete name (normalize)
path_parts <- unlist(strsplit(path, "/")) # Separate the path into its components
target_index <- which(path_parts == "Chapitre_3") # Search for "Chapitre_3" index
target_path <- paste(path_parts[1:target_index], collapse = "/") # Rebuild path to "Chapitre_3"
setwd(target_path) # Set the working directory


# |========================|
# |     LOAD LIBRARIES     |
# |========================|

library(fst)
library(dplyr)
library(lubridate)
library(ggplot2)
library(readxl)
library(stringr)
library(ggpubr)
library(tidyverse)
library(tidyr)
library(modi)
library(zoo)



# |===========================|
# |      IMPORT RAW DATA      |
# |===========================|

disch_sj <- read.csv("01_RawData/Temperature_Discharge_River/CEQEAU_Temperature_Discharge_1979_2024_SJ.csv")
disch_tr <- read.csv("01_RawData/Temperature_Discharge_River/CEQEAU_Temperature_Discharge_1979_2024_TR.csv")

prop_smolt <- as.data.frame(read_excel("~/Alex/Ecole/Doctorat/Chapitre_3/01_RawData/MELCCFP_Data/BDFUSION_2024.xlsm", sheet = 1, n_max = 91))

# Periods (growth)
gr_period_stj <- read_fst("02_PrepareData_Scripts/Intermediate_Cleaned_Data/Growth_Periods_StJean.fst")
gr_period_tri <- read_fst("02_PrepareData_Scripts/Intermediate_Cleaned_Data/Growth_Periods_Trinite.fst")

# Periods (hatch)
hatch_period <- data.frame(
  year = sort(unique(gr_period_stj$year)),
  start_doy = 98,
  end_doy = 136
)

# Periods (emergence)
em_period <- data.frame(
  year = sort(unique(gr_period_stj$year)),
  start_doy = 137,
  end_doy = 154
)


# |=========================================|
# |      ADAPT TEMPERATURES DATAFRAMES      |
# |=========================================|

# Change column names
colnames(disch_sj)[which(colnames(disch_sj) == "Observed.Discharge..m3.s.")] <- "Observed_Discharge_m3s"
colnames(disch_sj)[which(colnames(disch_sj) == "Simulated.Discharge..m3.s.")] <- "Simulated_Discharge_m3s"
colnames(disch_sj)[which(colnames(disch_sj) == "Tw..degree.C.")] <- "Temperature_C"

colnames(disch_tr)[which(colnames(disch_tr) == "Observed.Discharge..RatingCurve...m3.s.")] <- "Observed_Discharge_m3s"
colnames(disch_tr)[which(colnames(disch_tr) == "Simulated.Discharge..m3.s.")] <- "Simulated_Discharge_m3s"
colnames(disch_tr)[which(colnames(disch_tr) == "Simulated.Tw..degree.C.")] <- "Temperature_C"

colnames(disch_sj) <- tolower(colnames(disch_sj))
colnames(disch_tr) <- tolower(colnames(disch_tr))

# Change date column in class date 
old_locale <- Sys.getlocale("LC_TIME") # Saving the current locale
Sys.setlocale("LC_TIME", "C") # Set locale to C to use English months
disch_sj$date <- as.Date(disch_sj$date, format = "%d-%b-%Y") # Convert character to date
disch_tr$date <- as.Date(disch_tr$date, format = "%d-%b-%Y") # Convert character to date
Sys.setlocale("LC_TIME", old_locale) # Restores the original locale

# Keep only date and temperature data
disch_sj <- subset(disch_sj, select = c("date", "simulated_discharge_m3s"))
disch_tr <- subset(disch_tr, select = c("date", "simulated_discharge_m3s"))

# Add a column for year, month, day and day of year
disch_sj$year <- year(disch_sj$date)
disch_sj$month <- month(disch_sj$date)
disch_sj$day <- day(disch_sj$date)
disch_sj$doy <- yday(disch_sj$date)

disch_tr$year <- year(disch_tr$date)
disch_tr$month <- month(disch_tr$date)
disch_tr$day <- day(disch_tr$date)
disch_tr$doy <- yday(disch_tr$date)

annee <- 1994
ggplot(data = disch_sj %>% filter((month >= 5 & month <= 10) & year == annee), aes(x = date, y = simulated_discharge_m3s)) +
  geom_line()
ggplot(data = disch_tr %>% filter((month >= 6 & month <= 9) & year == annee), aes(x = date, y = simulated_discharge_m3s)) +
  geom_line()


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
  rename("river"="riviere", "year"="annee", "smolt_nb_2yr"="smolt_nb_en_devalaison_2_p", "smolt_nb_3yr"="smolt_nb_en_devalaison_3_p",
         "smolt_nb_4yr"="smolt_nb_en_devalaison_4_p", "smolt_nb_5yr"="smolt_nb_en_devalaison_5_p") 

# Make sure all years are there
prop_smolt_sub <- prop_smolt_sub %>%
  group_by(river) %>%
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
  filter(river == "stj") %>%
  arrange(year) %>%
  mutate(
    cohort = year,                 
    smolt_prop_2yr = lead(smolt_nb_2yr, n = 3), # 2yr smolts correspond to the egg-laying cohort of t-3 (t being the year of smoltification).
    smolt_prop_3yr = lead(smolt_nb_3yr, n = 4), # 3yr smolts correspond to the egg-laying cohort of t-4 (t being the year of smoltification).
    smolt_prop_4yr = lead(smolt_nb_4yr, n = 5), # 4yr smolts correspond to the egg-laying cohort of t-5 (t being the year of smoltification).
    smolt_prop_5yr = lead(smolt_nb_5yr, n = 6)  # 5yr smolts correspond to the egg-laying cohort of t-6 (t being the year of smoltification).
  ) %>%
  dplyr::select(cohort, smolt_prop_2yr, smolt_prop_3yr, smolt_prop_4yr, smolt_prop_5yr) %>%
  rowwise() %>%
  mutate(
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
  filter(river == "tri") %>%
  arrange(year) %>%
  mutate(
    cohort = year,                 
    smolt_prop_2yr = lead(smolt_nb_2yr, n = 3), # 2yr smolts correspond to the egg-laying cohort of t-3 (t being the year of smoltification).
    smolt_prop_3yr = lead(smolt_nb_3yr, n = 4), # 3yr smolts correspond to the egg-laying cohort of t-4 (t being the year of smoltification).
    smolt_prop_4yr = lead(smolt_nb_4yr, n = 5), # 4yr smolts correspond to the egg-laying cohort of t-5 (t being the year of smoltification).
    smolt_prop_5yr = lead(smolt_nb_5yr, n = 6)  # 5yr smolts correspond to the egg-laying cohort of t-6 (t being the year of smoltification).
  ) %>%
  dplyr::select(cohort, smolt_prop_2yr, smolt_prop_3yr, smolt_prop_4yr, smolt_prop_5yr) %>%
  rowwise() %>%
  mutate(
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


# |============================================================================|
# |      COHORT MEAN, MEDIAN AND VARIANCE DISCHARGE DURING GROWING SEASON      |
# |============================================================================|
# Average, median ans variance river water temperature during the growth period (between April and October) for the cohort, weighted by the percentage of individuals smoltifying each year.

# Weighted median function 
weighted_median <- function(x, w) {
  # Remove NA
  keep <- !is.na(x) & !is.na(w)
  x <- x[keep]
  w <- w[keep]
  
  # Order
  o <- order(x)
  x <- x[o]
  w <- w[o]
  
  # Calculation of the weighted median
  cum_w <- cumsum(w) / sum(w)
  x[which(cum_w >= 0.5)[1]]
}

weighted_co_disch <- function(smolt_data, discharge_data, periods){
  # Calculate average/median discharge during the desired period, as well as variance
  years <- unique(smolt_data$cohort)
  smolt_data$co_mean_disch <- NA
  smolt_data$co_median_disch <- NA
  smolt_data$co_var_disch <- NA
  for(i in seq_along(years)){
    smolt_data_i <- smolt_data %>%
      filter(cohort == years[i])
    
    # Note that i is the year of eggs deposition (not hatching), so growth season t+0 is the year t+1 in discharge dataset
    period_i <- periods %>% filter(year == smolt_data_i$cohort + 1) # Retrieve the specific period for year i
    if(nrow(period_i) == 0 || is.na(period_i$start_doy) || is.na(period_i$end_doy)){ # Verify that the period exists for this cohort.
      t_0 <- NA_real_
    }else{
      if(!years[i] + 1 %in% discharge_data$year){ # Check that the year exists in the discharge data
        t_0 <- NA_real_
      }else{
        t_0 <- discharge_data %>%
          dplyr::filter(year == years[i] + 1, doy >= period_i$start_doy, doy <= period_i$end_doy) %>%
          dplyr::pull(simulated_discharge_m3s)
        if (length(t_0) == 0) t_0 <- NA_real_ # If no data found in the interval
      }
    }
    # Following the same idea, growth season t+1 is the year t+2 in discharge dataset
    period_i <- periods %>% filter(year == smolt_data_i$cohort + 2) # Retrieve the specific period for year i
    if(nrow(period_i) == 0 || is.na(period_i$start_doy) || is.na(period_i$end_doy)){ # Verify that the period exists for this cohort.
      t_1 <- NA_real_
    }else{
      if(!years[i] + 2 %in% discharge_data$year){ # Check that the year exists in the discharge data
        t_1 <- NA_real_
      }else{
        t_1 <- discharge_data %>%
          dplyr::filter(year == years[i] + 2, doy >= period_i$start_doy, doy <= period_i$end_doy) %>%
          dplyr::pull(simulated_discharge_m3s)
        if (length(t_1) == 0) t_1 <- NA_real_ # If no data found in the interval
      }
    }
    # Following the same idea, growth season t+2 is the year t+3 in discharge dataset
    period_i <- periods %>% filter(year == smolt_data_i$cohort + 3) # Retrieve the specific period for year i
    if(nrow(period_i) == 0 || is.na(period_i$start_doy) || is.na(period_i$end_doy)){ # Verify that the period exists for this cohort.
      t_2 <- NA_real_
    }else{
      if(!years[i] + 3 %in% discharge_data$year){ # Check that the year exists in the discharge data
        t_2 <- NA_real_
      }else{
        t_2 <- discharge_data %>%
          dplyr::filter(year == years[i] + 3, doy >= period_i$start_doy, doy <= period_i$end_doy) %>%
          dplyr::pull(simulated_discharge_m3s)
        if (length(t_2) == 0) t_2 <- NA_real_ # If no data found in the interval
      }
    }
    # Following the same idea, growth season t+3 is the year t+4 in discharge dataset
    period_i <- periods %>% filter(year == smolt_data_i$cohort + 4) # Retrieve the specific period for year i
    if(nrow(period_i) == 0 || is.na(period_i$start_doy) || is.na(period_i$end_doy)){ # Verify that the period exists for this cohort.
      t_3 <- NA_real_
    }else{
      if(!years[i] + 4 %in% discharge_data$year){ # Check that the year exists in the discharge data
        t_3 <- NA_real_
      }else{
        t_3 <- discharge_data %>%
          dplyr::filter(year == years[i] + 4, doy >= period_i$start_doy, doy <= period_i$end_doy) %>%
          dplyr::pull(simulated_discharge_m3s)
        if (length(t_3) == 0) t_3 <- NA_real_ # If no data found in the interval
      }
    }
    # Following the same idea, growth season t+4 is the year t+5 in discharge dataset
    period_i <- periods %>% filter(year == smolt_data_i$cohort + 5) # Retrieve the specific period for year i
    if(nrow(period_i) == 0 || is.na(period_i$start_doy) || is.na(period_i$end_doy)){ # Verify that the period exists for this cohort.
      t_4 <- NA_real_
    }else{
      if(!years[i] + 5 %in% discharge_data$year){ # Check that the year exists in the discharge data
        t_4 <- NA_real_
      }else{
        t_4 <- discharge_data %>%
          dplyr::filter(year == years[i] + 5, doy >= period_i$start_doy, doy <= period_i$end_doy) %>%
          dplyr::pull(simulated_discharge_m3s)
        if (length(t_4) == 0) t_4 <- NA_real_ # If no data found in the interval
      }
    }
    # Following the same idea, growth season t+5 is the year t+6 in discharge dataset
    period_i <- periods %>% filter(year == smolt_data_i$cohort + 6) # Retrieve the specific period for year i
    if(nrow(period_i) == 0 || is.na(period_i$start_doy) || is.na(period_i$end_doy)){ # Verify that the period exists for this cohort.
      t_5 <- NA_real_
    }else{
      if(!years[i] + 6 %in% discharge_data$year){ # Check that the year exists in the discharge data
        t_5 <- NA_real_
      }else{
        t_5 <- discharge_data %>%
          dplyr::filter(year == years[i] + 6, doy >= period_i$start_doy, doy <= period_i$end_doy) %>%
          dplyr::pull(simulated_discharge_m3s)
        if (length(t_5) == 0) t_5 <- NA_real_ # If no data found in the interval
      }
    }
    
    # Extract the proportions, replacing NA values with 0
    p2 <- ifelse(is.na(smolt_data_i$smolt_prop_2yr), 0, smolt_data_i$smolt_prop_2yr)
    p3 <- ifelse(is.na(smolt_data_i$smolt_prop_3yr), 0, smolt_data_i$smolt_prop_3yr)
    p4 <- ifelse(is.na(smolt_data_i$smolt_prop_4yr), 0, smolt_data_i$smolt_prop_4yr)
    p5 <- ifelse(is.na(smolt_data_i$smolt_prop_5yr), 0, smolt_data_i$smolt_prop_5yr)
    
    # Calculate weighted mean
    smolt_data$co_mean_disch[i] <- weighted.mean(
      x = c(t_0, t_1, t_2, t_3, t_4, t_5),
      w = c(rep(1, length(t_0)), # weight = 1 because all individuals from the cohort i were exposed to year t
            rep(1, length(t_1)), # weight = 1 because all individuals from the cohort i were exposed to year t
            rep(max(c(0, 1-p2)),          length(t_2)), 
            rep(max(c(0, 1-p2-p3)),       length(t_3)),
            rep(max(c(0, 1-p2-p3-p4)),    length(t_4)),
            rep(max(c(0, 1-p2-p3-p4-p5)), length(t_5))
      ),
      na.rm = FALSE
    )
    
    # Calculate weighted median
    smolt_data$co_median_disch[i] <- weighted_median(
      x = c(t_0, t_1, t_2, t_3, t_4, t_5),
      w = c(rep(1, length(t_0)), # weight = 1 because all individuals from the cohort i were exposed to year t
            rep(1, length(t_1)), # weight = 1 because all individuals from the cohort i were exposed to year t
            rep(max(c(0, 1-p2)),          length(t_2)), 
            rep(max(c(0, 1-p2-p3)),       length(t_3)),
            rep(max(c(0, 1-p2-p3-p4)),    length(t_4)),
            rep(max(c(0, 1-p2-p3-p4-p5)), length(t_5))
      )
    )
    
    # Calculate weighted variation
    x_all <- c(t_0, t_1, t_2, t_3, t_4, t_5)
    w_all <- c(
      rep(1, length(t_0)),
      rep(1, length(t_1)),
      rep(max(c(0, 1-p2)),          length(t_2)), 
      rep(max(c(0, 1-p2-p3)),       length(t_3)),
      rep(max(c(0, 1-p2-p3-p4)),    length(t_4)),
      rep(max(c(0, 1-p2-p3-p4-p5)), length(t_5))
    )
    if(length(x_all) > 0 && length(w_all) > 0){
      smolt_data$co_var_disch[i] <- weighted.var(
        x = x_all,
        w = w_all,
        na.rm = FALSE
      )
    } else {
      smolt_data$co_var_disch[i] <- NA
    }
  }
  
  # Remove years with just NA
  smolt_data <- smolt_data %>%
    dplyr::filter(!(is.na(smolt_data$smolt_prop_2yr) & 
                      is.na(smolt_data$smolt_prop_3yr) &
                      is.na(smolt_data$smolt_prop_4yr) & 
                      is.na(smolt_data$smolt_prop_5yr))) %>%
    dplyr::select(cohort, co_mean_disch, co_median_disch, co_var_disch) %>%
    as.data.frame()
  
  return(smolt_data)
}

weighted_co_season_sj <- weighted_co_disch(
  smolt_data = prop_smolt_sub_sj,
  discharge_data = disch_sj,
  periods = gr_period_stj
) %>%
  rename(co_gr_mean_disc = co_mean_disch, co_gr_median_disc = co_median_disch, co_gr_var_disc = co_var_disch)

weighted_co_season_tr <- weighted_co_disch(
  smolt_data = prop_smolt_sub_tr,
  discharge_data = disch_tr,
  periods = gr_period_tri
) %>%
  rename(co_gr_mean_disc = co_mean_disch, co_gr_median_disc = co_median_disch, co_gr_var_disc = co_var_disch)

ggarrange(
  ggarrange(
    ggplot(data = weighted_co_season_sj, aes(x = cohort, y = co_gr_mean_disc)) +
      geom_point() +
      geom_line() +
      labs(title = "St-Jean River"),
    ggplot(data = weighted_co_season_tr, aes(x = cohort, y = co_gr_mean_disc)) +
      geom_point() +
      geom_line() +
      labs(title = "Trinite River"),
    ncol = 2),
  
  ggarrange(
    ggplot(data = weighted_co_season_sj, aes(x = cohort, y = co_gr_median_disc)) +
      geom_point() +
      geom_line() +
      labs(title = "St-Jean River"),
    ggplot(data = weighted_co_season_tr, aes(x = cohort, y = co_gr_median_disc)) +
      geom_point() +
      geom_line() +
      labs(title = "Trinite River"),
    ncol = 2),
  ncol = 1)

ggarrange(
  ggplot(data = weighted_co_season_sj, aes(x = cohort, y = co_gr_var_disc)) +
    geom_point() +
    geom_line() +
    labs(title = "St-Jean River"),
  ggplot(data = weighted_co_season_tr, aes(x = cohort, y = co_gr_var_disc)) +
    geom_point() +
    geom_line() +
    labs(title = "Trinite River"),
  ncol = 2)


# |===========================================================|
# |      COHORT MEDIAN DISCHARGE DURING GRAVEL EMERGENCE      |
# |===========================================================|
# Average/median during key periods:
### Gravel emergence period -> one value per cohort

emergence_sj <- disch_sj %>%
  dplyr::filter(doy >= unique(em_period$start_doy), doy <= unique(em_period$end_doy)) %>%
  dplyr::rename(cohort = year) %>%
  dplyr::group_by(cohort) %>%
  dplyr::summarise(em_mean_disc = mean(simulated_discharge_m3s),
                   em_median_disc = median(simulated_discharge_m3s),
                   em_var_disc = var(simulated_discharge_m3s))

emergence_tr <- disch_tr %>%
  dplyr::filter(doy >= unique(em_period$start_doy), doy <= unique(em_period$end_doy)) %>%
  dplyr::rename(cohort = year) %>%
  dplyr::group_by(cohort) %>%
  dplyr::summarise(em_mean_disc = mean(simulated_discharge_m3s),
                   em_median_disc = median(simulated_discharge_m3s),
                   em_var_disc = var(simulated_discharge_m3s))

ggarrange(
  ggplot(data = emergence_sj, aes(x = cohort, y = em_mean_disc)) +
    geom_point() +
    geom_line() +
    labs(title = "St-Jean River"),
  ggplot(data = emergence_tr, aes(x = cohort, y = em_mean_disc)) +
    geom_point() +
    geom_line() +
    labs(title = "Trinite River"),
  ncol = 2)

ggarrange(
  ggplot(data = emergence_sj, aes(x = cohort, y = em_median_disc)) +
    geom_point() +
    geom_line() +
    labs(title = "St-Jean River"),
  ggplot(data = emergence_tr, aes(x = cohort, y = em_median_disc)) +
    geom_point() +
    geom_line() +
    labs(title = "Trinite River"),
  ncol = 2)


# |==========================================================================================================|
# |      YEARLY 75TH PERCENTILE AVERAGE DAILY FLOW, BASED ON OBEDZINSKI ET AL. (2018): DURING EMERGENCE      |
# |==========================================================================================================|

qt75_em_sj <- disch_sj %>%
  dplyr::filter(doy >= unique(em_period$start_doy), doy <= unique(em_period$end_doy)) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(em_quant_75 = as.numeric(quantile(simulated_discharge_m3s, 0.75))) %>%
  dplyr::rename(cohort = year)

qt75_em_tr <- disch_tr %>%
  dplyr::filter(doy >= unique(em_period$start_doy), doy <= unique(em_period$end_doy)) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(em_quant_75 = as.numeric(quantile(simulated_discharge_m3s, 0.75))) %>%
  dplyr::rename(cohort = year)

ggarrange(
  ggplot(data = qt75_em_sj, aes(x = cohort, y = em_quant_75)) +
    geom_line() +
    labs(title = "St. Jean"),
  ggplot(data = qt75_em_tr, aes(x = cohort, y = em_quant_75)) +
    geom_line() +
    labs(title = "Trinite"),
  ncol = 1)


# |=======================================================================================================|
# |      YEARLY 75TH PERCENTILE AVERAGE DAILY FLOW, BASED ON OBEDZINSKI ET AL. (2018): DURING GROWTH      |
# |=======================================================================================================|

weighted_co_disch_75qt <- function(smolt_data, discharge_data, periods){
  # Calculate average/median discharge during the desired period, as well as variance
  years <- unique(smolt_data$cohort)
  smolt_data$co_mean_disch <- NA
  for(i in seq_along(years)){
    smolt_data_i <- smolt_data %>%
      filter(cohort == years[i])
    
    # Note that i is the year of eggs deposition (not hatching), so growth season t+0 is the year t+1 in temperatures dataset
    period_i <- periods %>% filter(year == smolt_data_i$cohort + 1) # Retrieve the specific period for year i
    if(nrow(period_i) == 0 || is.na(period_i$start_doy) || is.na(period_i$end_doy)){ # Verify that the period exists for this cohort.
      t_0 <- NA_real_
    }else{
      if(!years[i] + 1 %in% discharge_data$year){ # Check that the year exists in the temperature data
        t_0 <- NA_real_
      }else{
        t_0 <- discharge_data %>%
          dplyr::filter(year == years[i] + 1, doy >= period_i$start_doy, doy <= period_i$end_doy) %>%
          dplyr::pull(simulated_discharge_m3s)
        t_0 <- as.numeric(quantile(t_0, 0.75))
        if (length(t_0) == 0) t_0 <- NA_real_ # If no data found in the interval
      }
    }
    # Following the same idea, growth season t+1 is the year t+2 in temperatures dataset
    period_i <- periods %>% filter(year == smolt_data_i$cohort + 2) # Retrieve the specific period for year i
    if(nrow(period_i) == 0 || is.na(period_i$start_doy) || is.na(period_i$end_doy)){ # Verify that the period exists for this cohort.
      t_1 <- NA_real_
    }else{
      if(!years[i] + 2 %in% discharge_data$year){ # Check that the year exists in the temperature data
        t_1 <- NA_real_
      }else{
        t_1 <- discharge_data %>%
          dplyr::filter(year == years[i] + 2, doy >= period_i$start_doy, doy <= period_i$end_doy) %>%
          dplyr::pull(simulated_discharge_m3s)
        t_1 <- as.numeric(quantile(t_1, 0.75))
        if (length(t_1) == 0) t_1 <- NA_real_ # If no data found in the interval
      }
    }
    # Following the same idea, growth season t+2 is the year t+3 in temperatures dataset
    period_i <- periods %>% filter(year == smolt_data_i$cohort + 3) # Retrieve the specific period for year i
    if(nrow(period_i) == 0 || is.na(period_i$start_doy) || is.na(period_i$end_doy)){ # Verify that the period exists for this cohort.
      t_2 <- NA_real_
    }else{
      if(!years[i] + 3 %in% discharge_data$year){ # Check that the year exists in the temperature data
        t_2 <- NA_real_
      }else{
        t_2 <- discharge_data %>%
          dplyr::filter(year == years[i] + 3, doy >= period_i$start_doy, doy <= period_i$end_doy) %>%
          dplyr::pull(simulated_discharge_m3s)
        t_2 <- as.numeric(quantile(t_2, 0.75))
        if (length(t_2) == 0) t_2 <- NA_real_ # If no data found in the interval
      }
    }
    # Following the same idea, growth season t+3 is the year t+4 in temperatures dataset
    period_i <- periods %>% filter(year == smolt_data_i$cohort + 4) # Retrieve the specific period for year i
    if(nrow(period_i) == 0 || is.na(period_i$start_doy) || is.na(period_i$end_doy)){ # Verify that the period exists for this cohort.
      t_3 <- NA_real_
    }else{
      if(!years[i] + 4 %in% discharge_data$year){ # Check that the year exists in the temperature data
        t_3 <- NA_real_
      }else{
        t_3 <- discharge_data %>%
          dplyr::filter(year == years[i] + 4, doy >= period_i$start_doy, doy <= period_i$end_doy) %>%
          dplyr::pull(simulated_discharge_m3s)
        t_3 <- as.numeric(quantile(t_3, 0.75))
        if (length(t_3) == 0) t_3 <- NA_real_ # If no data found in the interval
      }
    }
    # Following the same idea, growth season t+4 is the year t+5 in temperatures dataset
    period_i <- periods %>% filter(year == smolt_data_i$cohort + 5) # Retrieve the specific period for year i
    if(nrow(period_i) == 0 || is.na(period_i$start_doy) || is.na(period_i$end_doy)){ # Verify that the period exists for this cohort.
      t_4 <- NA_real_
    }else{
      if(!years[i] + 5 %in% discharge_data$year){ # Check that the year exists in the temperature data
        t_4 <- NA_real_
      }else{
        t_4 <- discharge_data %>%
          dplyr::filter(year == years[i] + 5, doy >= period_i$start_doy, doy <= period_i$end_doy) %>%
          dplyr::pull(simulated_discharge_m3s)
        t_4 <- as.numeric(quantile(t_4, 0.75))
        if (length(t_4) == 0) t_4 <- NA_real_ # If no data found in the interval
      }
    }
    # Following the same idea, growth season t+5 is the year t+6 in temperatures dataset
    period_i <- periods %>% filter(year == smolt_data_i$cohort + 6) # Retrieve the specific period for year i
    if(nrow(period_i) == 0 || is.na(period_i$start_doy) || is.na(period_i$end_doy)){ # Verify that the period exists for this cohort.
      t_5 <- NA_real_
    }else{
      if(!years[i] + 6 %in% discharge_data$year){ # Check that the year exists in the temperature data
        t_5 <- NA_real_
      }else{
        t_5 <- discharge_data %>%
          dplyr::filter(year == years[i] + 6, doy >= period_i$start_doy, doy <= period_i$end_doy) %>%
          dplyr::pull(simulated_discharge_m3s)
        t_5 <- as.numeric(quantile(t_5, 0.75))
        if(length(t_5) == 0) t_5 <- NA_real_ # If no data found in the interval
      }
    }
    
    # Extract the proportions, replacing NA values with 0
    p2 <- ifelse(is.na(smolt_data_i$smolt_prop_2yr), 0, smolt_data_i$smolt_prop_2yr)
    p3 <- ifelse(is.na(smolt_data_i$smolt_prop_3yr), 0, smolt_data_i$smolt_prop_3yr)
    p4 <- ifelse(is.na(smolt_data_i$smolt_prop_4yr), 0, smolt_data_i$smolt_prop_4yr)
    p5 <- ifelse(is.na(smolt_data_i$smolt_prop_5yr), 0, smolt_data_i$smolt_prop_5yr)
    
    # Calculate weighted mean
    smolt_data$co_mean_disch[i] <- weighted.mean(
      x = c(t_0, t_1, t_2, t_3, t_4, t_5),
      w = c(rep(1, length(t_0)), # weight = 1 because all individuals from the cohort i were exposed to year t
            rep(1, length(t_1)), # weight = 1 because all individuals from the cohort i were exposed to year t
            rep(max(c(0, 1-p2)),          length(t_2)), 
            rep(max(c(0, 1-p2-p3)),       length(t_3)),
            rep(max(c(0, 1-p2-p3-p4)),    length(t_4)),
            rep(max(c(0, 1-p2-p3-p4-p5)), length(t_5))
      ),
      na.rm = FALSE
    )
  }
  
  # Remove years with just NA
  smolt_data <- smolt_data %>%
    dplyr::filter(!(is.na(smolt_data$smolt_prop_2yr) & 
                      is.na(smolt_data$smolt_prop_3yr) &
                      is.na(smolt_data$smolt_prop_4yr) & 
                      is.na(smolt_data$smolt_prop_5yr))) %>%
    dplyr::select(cohort, co_mean_disch) %>%
    as.data.frame()
  
  return(smolt_data)
}

qt75_gr_co_sj <- weighted_co_disch_75qt(
  smolt_data = prop_smolt_sub_sj,
  discharge_data = disch_sj,
  periods = gr_period_stj
) %>%
  rename(co_gr_mean_75qt_disc = co_mean_disch)

qt75_gr_co_tr <- weighted_co_disch_75qt(
  smolt_data = prop_smolt_sub_tr,
  discharge_data = disch_tr,
  periods = gr_period_tri
) %>%
  rename(co_gr_mean_75qt_disc = co_mean_disch)


ggarrange(
  ggplot(data = qt75_gr_co_sj, aes(x = cohort, y = co_gr_mean_75qt_disc)) +
    geom_line() +
    labs(title = "St. Jean"),
  ggplot(data = qt75_gr_co_tr, aes(x = cohort, y = co_gr_mean_75qt_disc)) +
    geom_line() +
    labs(title = "Trinite"),
  ncol = 1)



# |=======================================================================================|
# |      MINIMUM 30-D AVERAGE FLOW, BASED ON OBEDZINSKI ET AL. (2018): DURING GROWTH      |
# |=======================================================================================|

# Minimum 30-d average flow:
### 1- Create as many 30-day windows as possible within the period of interest
### 2- Calculate the AVERAGE discharge over these 30 days in each window, which gives as many values as there are windows. So for one window: mean(30 values).
### 3- Take the lowest value of the windows, which gives a value per year (interest period in the year).

weighted_co_disch <- function(smolt_data, discharge_data, periods){
  # Calculate average/median discharge during the desired period, as well as variance
  years <- unique(smolt_data$cohort)
  smolt_data$co_30d_mean_disch <- NA
  for(i in seq_along(years)){
    smolt_data_i <- smolt_data %>%
      filter(cohort == years[i])
    
    # Note that i is the year of eggs deposition (not hatching), so growth season t+0 is the year t+1 in temperatures dataset
    period_i <- periods %>% filter(year == smolt_data_i$cohort + 1) # Retrieve the specific period for year i
    if(nrow(period_i) == 0 || is.na(period_i$start_doy) || is.na(period_i$end_doy)){ # Verify that the period exists for this cohort.
      t_0 <- NA_real_
    }else{
      if(!years[i] + 1 %in% discharge_data$year){ # Check that the year exists in the temperature data
        t_0 <- NA_real_
      }else{
        t_0 <- discharge_data %>%
          dplyr::filter(year == years[i] + 1, doy >= period_i$start_doy, doy <= period_i$end_doy) %>%
          dplyr::pull(simulated_discharge_m3s)
        t_0 <- rollapply(t_0, width = 30, FUN = mean, fill = NA,
                         align = "right"   # the value is associated with the last day of the window
        )
        t_0 <- as.numeric(na.omit(t_0)) # Remove the first NA, because they represent the first windows that are under 30 days
        t_0 <- min(t_0)
        if (length(t_0) == 0) t_0 <- NA_real_ # If no data found in the interval
      }
    }
    # Following the same idea, growth season t+1 is the year t+2 in temperatures dataset
    period_i <- periods %>% filter(year == smolt_data_i$cohort + 2) # Retrieve the specific period for year i
    if(nrow(period_i) == 0 || is.na(period_i$start_doy) || is.na(period_i$end_doy)){ # Verify that the period exists for this cohort.
      t_1 <- NA_real_
    }else{
      if(!years[i] + 2 %in% discharge_data$year){ # Check that the year exists in the temperature data
        t_1 <- NA_real_
      }else{
        t_1 <- discharge_data %>%
          dplyr::filter(year == years[i] + 2, doy >= period_i$start_doy, doy <= period_i$end_doy) %>%
          dplyr::pull(simulated_discharge_m3s)
        t_1 <- rollapply(t_1, width = 30, FUN = mean, fill = NA,
                         align = "right"   # the value is associated with the last day of the window
        )
        t_1 <- as.numeric(na.omit(t_1)) # Remove the first NA, because they represent the first windows that are under 30 days
        t_1 <- min(t_1)
        if (length(t_1) == 0) t_1 <- NA_real_ # If no data found in the interval
      }
    }
    # Following the same idea, growth season t+2 is the year t+3 in temperatures dataset
    period_i <- periods %>% filter(year == smolt_data_i$cohort + 3) # Retrieve the specific period for year i
    if(nrow(period_i) == 0 || is.na(period_i$start_doy) || is.na(period_i$end_doy)){ # Verify that the period exists for this cohort.
      t_2 <- NA_real_
    }else{
      if(!years[i] + 3 %in% discharge_data$year){ # Check that the year exists in the temperature data
        t_2 <- NA_real_
      }else{
        t_2 <- discharge_data %>%
          dplyr::filter(year == years[i] + 3, doy >= period_i$start_doy, doy <= period_i$end_doy) %>%
          dplyr::pull(simulated_discharge_m3s)
        t_2 <- rollapply(t_2, width = 30, FUN = mean, fill = NA,
                         align = "right"   # the value is associated with the last day of the window
        )
        t_2 <- as.numeric(na.omit(t_2)) # Remove the first NA, because they represent the first windows that are under 30 days
        t_2 <- min(t_2)
        if (length(t_2) == 0) t_2 <- NA_real_ # If no data found in the interval
      }
    }
    # Following the same idea, growth season t+3 is the year t+4 in temperatures dataset
    period_i <- periods %>% filter(year == smolt_data_i$cohort + 4) # Retrieve the specific period for year i
    if(nrow(period_i) == 0 || is.na(period_i$start_doy) || is.na(period_i$end_doy)){ # Verify that the period exists for this cohort.
      t_3 <- NA_real_
    }else{
      if(!years[i] + 4 %in% discharge_data$year){ # Check that the year exists in the temperature data
        t_3 <- NA_real_
      }else{
        t_3 <- discharge_data %>%
          dplyr::filter(year == years[i] + 4, doy >= period_i$start_doy, doy <= period_i$end_doy) %>%
          dplyr::pull(simulated_discharge_m3s)
        t_3 <- rollapply(t_3, width = 30, FUN = mean, fill = NA,
                         align = "right"   # the value is associated with the last day of the window
        )
        t_3 <- as.numeric(na.omit(t_3)) # Remove the first NA, because they represent the first windows that are under 30 days
        t_3 <- min(t_3)
        if (length(t_3) == 0) t_3 <- NA_real_ # If no data found in the interval
      }
    }
    # Following the same idea, growth season t+4 is the year t+5 in temperatures dataset
    period_i <- periods %>% filter(year == smolt_data_i$cohort + 5) # Retrieve the specific period for year i
    if(nrow(period_i) == 0 || is.na(period_i$start_doy) || is.na(period_i$end_doy)){ # Verify that the period exists for this cohort.
      t_4 <- NA_real_
    }else{
      if(!years[i] + 5 %in% discharge_data$year){ # Check that the year exists in the temperature data
        t_4 <- NA_real_
      }else{
        t_4 <- discharge_data %>%
          dplyr::filter(year == years[i] + 5, doy >= period_i$start_doy, doy <= period_i$end_doy) %>%
          dplyr::pull(simulated_discharge_m3s)
        t_4 <- rollapply(t_4, width = 30, FUN = mean, fill = NA,
                         align = "right"   # the value is associated with the last day of the window
        )
        t_4 <- as.numeric(na.omit(t_4)) # Remove the first NA, because they represent the first windows that are under 30 days
        t_4 <- min(t_4)
        if (length(t_4) == 0) t_4 <- NA_real_ # If no data found in the interval
      }
    }
    # Following the same idea, growth season t+5 is the year t+6 in temperatures dataset
    period_i <- periods %>% filter(year == smolt_data_i$cohort + 6) # Retrieve the specific period for year i
    if(nrow(period_i) == 0 || is.na(period_i$start_doy) || is.na(period_i$end_doy)){ # Verify that the period exists for this cohort.
      t_5 <- NA_real_
    }else{
      if(!years[i] + 6 %in% discharge_data$year){ # Check that the year exists in the temperature data
        t_5 <- NA_real_
      }else{
        t_5 <- discharge_data %>%
          dplyr::filter(year == years[i] + 6, doy >= period_i$start_doy, doy <= period_i$end_doy) %>%
          dplyr::pull(simulated_discharge_m3s)
        t_5 <- rollapply(t_5, width = 30, FUN = mean, fill = NA,
                         align = "right"   # the value is associated with the last day of the window
        )
        t_5 <- as.numeric(na.omit(t_5)) # Remove the first NA, because they represent the first windows that are under 30 days
        t_5 <- min(t_5)
        if (length(t_5) == 0) t_5 <- NA_real_ # If no data found in the interval
      }
    }
    
    # Extract the proportions, replacing NA values with 0
    p2 <- ifelse(is.na(smolt_data_i$smolt_prop_2yr), 0, smolt_data_i$smolt_prop_2yr)
    p3 <- ifelse(is.na(smolt_data_i$smolt_prop_3yr), 0, smolt_data_i$smolt_prop_3yr)
    p4 <- ifelse(is.na(smolt_data_i$smolt_prop_4yr), 0, smolt_data_i$smolt_prop_4yr)
    p5 <- ifelse(is.na(smolt_data_i$smolt_prop_5yr), 0, smolt_data_i$smolt_prop_5yr)
    
    # Calculate weighted mean
    smolt_data$co_30d_mean_disch[i] <- weighted.mean(
      x = c(t_0, t_1, t_2, t_3, t_4, t_5),
      w = c(rep(1, length(t_0)), # weight = 1 because all individuals from the cohort i were exposed to year t
            rep(1, length(t_1)), # weight = 1 because all individuals from the cohort i were exposed to year t
            rep(max(c(0, 1-p2)),          length(t_2)), 
            rep(max(c(0, 1-p2-p3)),       length(t_3)),
            rep(max(c(0, 1-p2-p3-p4)),    length(t_4)),
            rep(max(c(0, 1-p2-p3-p4-p5)), length(t_5))
      ),
      na.rm = FALSE
    )
  }
  
  # Remove years with just NA
  smolt_data <- smolt_data %>%
    dplyr::filter(!(is.na(smolt_data$smolt_prop_2yr) & 
                      is.na(smolt_data$smolt_prop_3yr) &
                      is.na(smolt_data$smolt_prop_4yr) & 
                      is.na(smolt_data$smolt_prop_5yr))) %>%
    dplyr::select(cohort, co_30d_mean_disch) %>%
    as.data.frame()
  
  return(smolt_data)
}

co_30d_mean_disch_sj <- weighted_co_disch(
  smolt_data = prop_smolt_sub_sj,
  discharge_data = disch_sj,
  periods = gr_period_stj
)

co_30d_mean_disch_tr <- weighted_co_disch(
  smolt_data = prop_smolt_sub_tr,
  discharge_data = disch_tr,
  periods = gr_period_tri
)

ggarrange(
  ggplot(data = co_30d_mean_disch_sj, aes(x = cohort, y = co_30d_mean_disch)) +
    geom_line() +
    labs(title = "St. Jean"),
  ggplot(data = co_30d_mean_disch_tr, aes(x = cohort, y = co_30d_mean_disch)) +
    geom_line() +
    labs(title = "Trinite"),
  ncol = 1)


# |=============================================|
# |      EXPORT COMBINED DISCHARGE RESULTS      |
# |=============================================|

# Combine all results
sj_all <- Reduce(function(x, y) merge(x, y, by = "cohort", all = TRUE), 
                 list(weighted_co_season_sj, emergence_sj, qt75_em_sj, qt75_gr_co_sj, co_30d_mean_disch_sj))
tr_all <- Reduce(function(x, y) merge(x, y, by = "cohort", all = TRUE), 
                 list(weighted_co_season_tr, emergence_tr, qt75_em_tr, qt75_gr_co_tr, co_30d_mean_disch_tr))

sj_all$river <- "stj"
tr_all$river <- "tri"

data_all <- bind_rows(sj_all, tr_all) %>%
  dplyr::select(river, cohort, everything())

write_fst(data_all, path = "03_Clean_Data/River/River_Discharge.fst")






#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################
#############################################################################################################################










