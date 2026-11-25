# Determinants of interannual survival variation in two Atlantic salmon populations monitored over three decades


## Table of contents
- [Repository structure](#repository-structure)
- [1. Raw data](#1-raw-data)
  - [1.1 Salmon data (MELCCFP)](#11-salmon-data-melccfp)
  - [1.2 Discharge and water temperature (CEQUEAU model)](#12-discharge-and-water-temperature-cequeau-model)
  - [1.3 Air temperature](#13-air-temperature)
  - [1.4 Forestry activity](#14-forestry-activity)
  - [1.5 Food availability at sea](#15-food-availability-at-sea)
  - [1.6 Acoustic telemetry](#16-acoustic-telemetry)
  - [1.7 Sea surface temperature](#17-sea-surface-temperature)
  - [1.8 Fisheries and abundance at sea (ICES WGNAS)](#18-fisheries-and-abundance-at-sea-ices-wgnas)
- [2. Data preparation and variable construction](#2-data-preparation-and-variable-construction)
  - [2.1 Shared data](#21-shared-data)
  - [2.2 Freshwater variables](#22-freshwater-variables)
    - [2.2.1 Discharge](#221-discharge)
    - [2.2.2 Forestry activity](#222-forestry-activity)
    - [2.2.3 Temperature](#223-temperature)
    - [2.2.4 Winter temperature (freezing degree-days)](#224-winter-temperature-freezing-degree-days)
  - [2.3 Marine variables](#23-marine-variables)
    - [2.3.1 Thermal conditions](#231-thermal-conditions)
    - [2.3.2 Food availability](#232-food-availability)
    - [2.3.3 Temperature difference between sea and river at sea entry](#233-temperature-difference-between-sea-and-river-at-sea-entry)
    - [2.3.4 Exploitation rates at sea](#234-exploitation-rates-at-sea)
    - [2.3.5 Pre-Fishery abundance (density at sea)](#235-pre-fishery-abundance-density-at-sea)
- [3. Modelling](#3-modelling)
  - [3.1 Freshwater survival](#31-freshwater-survival)
  - [3.2 Return from sea](#32-return-from-sea)
- [4. Software environment](#4-software-environment)


## Repository structure

```
phd-chapter3-Atlantic-Salmon/
└── R/
    ├── 00_raw_data/       Raw data (section 1)
    ├── 01_derived_data/   Intermediate data produced by the scripts
    ├── 02_prepare_data/   Data cleaning and variable construction scripts (section 2)
    ├── 03_clean_data/     Cleaned variables, by theme
    └── 04_analysis/       Dataset assembly, modelling and figures (section 3)
```

Note: intermediate data are not described in a separate section. Each one is listed in the outputs of the script that produces it, with its publication status.

---

## 1. Raw data

Folder: `R/00_raw_data/`

Raw data files from which the variables are built. Each dataset is described with the same template: file(s), description, columns used (if applicable), source and publication.


### 1.1 Salmon data (MELCCFP)

#### Smolt monitoring, survival and return rates — St. Jean and Trinité

- **File(s)**: `R/00_raw_data/salmon/Salmon_Monitoring_StJean_Trinite_2024.csv`
- **Description**: annual monitoring data for the salmon of the St. Jean and Trinité rivers, extracted from the tables of the MELCCFP 2024 annual report: number of downstream-migrating smolts by age class, smolt condition factor, egg-to-smolt survival and smolt-to-adult return rates.
- **Columns**:

  | Column | Description | Report table |
  |---|---|---|
  | `river` | River (`stj` = St. Jean, `tri` = Trinité) | — |
  | `year` | Year, whose meaning depends on the column: year of sea entry (smolt migration) for `smolt_n_age2` to `smolt_n_age5`, `smolt_k`, `return_1sw_pct`, `return_2sw_pct`, `return_3sw_pct` and `return_all_pct`; spawning year of the cohort for `surv_egg_smolt_pct` | — |
  | `smolt_n_age2` to `smolt_n_age5` | Number of downstream-migrating smolts of age 2+, 3+, 4+ and 5+ | Table 5 |
  | `smolt_k` | Mean Fulton's condition factor of the collected smolts | Table 5 |
  | `surv_egg_smolt_pct` | Egg-to-smolt survival (%) | Tables 15 (St. Jean) and 16 (Trinité) |
  | `return_1sw_pct` | Smolt-to-adult return rate, one-sea-winter (1SW) (%) | Tables 17 (St. Jean) and 18 (Trinité) |
  | `return_2sw_pct` | Smolt-to-adult return rate, two-sea-winter (2SW) (%) | Tables 17 (St. Jean) and 18 (Trinité) |
  | `return_3sw_pct` | Smolt-to-adult return rate, three-sea-winter (3SW) (%) | Tables 17 (St. Jean) and 18 (Trinité) |
  | `return_all_pct` | Smolt-to-adult return rate, all sea ages combined (%) | Tables 17 (St. Jean) and 18 (Trinité) |

Caution: depending on the column, `year` is either the year of sea entry or the spawning year (see the column table).

- **Source**: Ministère de l'Environnement, de la Lutte contre les changements climatiques, de la Faune et des Parcs (MELCCFP). Cauchon, V., April, J., 2025. Suivi des populations témoins de saumon atlantique au Québec: Rapport scientifique 2024. Ministère des Forêts, de la Faune et des Parcs, Secteur de la faune et des parcs.
- **Publication**: Published

#### Daily smolt captures at the traps — St. Jean and Trinité

- **File(s)**:
  - St. Jean: `R/00_raw_data/salmon/Smolt_Daily_Captures_StJean.csv`
  - Trinité: `R/00_raw_data/salmon/Smolt_Daily_Captures_Trinite.csv`
- **Description**: number of smolts captured per day at the downstream traps, for each year, extracted from the MELCCFP 2024 annual report (Figures 4 and 5). One row per day of trap operation: days without operation are absent and days without captures are 0. Years when the trap was not operated (e.g., St. Jean 2015) are absent. Used to compute the median migration dates and the date quantiles (migration windows).
- **Columns**:

  | Column | Description |
  |---|---|
  | `date` | Capture date (YYYY-MM-DD) |
  | `cohort` | Year of downstream migration |
  | `doy` | Day of year |
  | `captures` | Number of smolts captured that day |

- **Source**: Ministère de l'Environnement, de la Lutte contre les changements climatiques, de la Faune et des Parcs (MELCCFP). Cauchon, V., April, J., 2025. Suivi des populations témoins de saumon atlantique au Québec: Rapport scientifique 2024. Ministère des Forêts, de la Faune et des Parcs, Secteur de la faune et des parcs.
- **Publication**: Published.

### 1.2 Discharge and water temperature (CEQUEAU model)

#### Simulated discharge and water temperature — St. Jean

- **File(s)**: `R/00_raw_data/Temperature_Discharge_River/CEQEAU_Temperature_Discharge_1979_2024_SJ.csv`
- **Description**: daily discharge and water temperature simulated by the CEQUEAU model, 1979–2024.
- **Columns used**:

  | Column | Description |
  |---|---|
  | `Date` | Date |
  | `Observed Discharge (m3/s)` | Observed discharge (m³/s), when available, at the river station (data from RivTemp) |
  | `Simulated Discharge (m3/s)` | Discharge simulated by CEQUEAU (m³/s) |
  | `Tw (degree C)` | Water temperature simulated by CEQUEAU (°C) |

- **Publication**:

#### Simulated discharge and water temperature — Trinité

- **File(s)**: `R/00_raw_data/Temperature_Discharge_River/CEQEAU_Temperature_Discharge_1979_2024_TR.csv`
- **Description**: daily discharge and water temperature simulated by the CEQUEAU model, 1979–2024.
- **Columns used**:

  | Column | Description |
  |---|---|
  | `Date` | Date |
  | `Observed Discharge (RatingCurve) (m3/s)` | Observed discharge (m³/s), when available, computed from a rating curve established with the river level and the discharges sampled in summer 2025 |
  | `Simulated Discharge (m3/s)` | Discharge simulated by CEQUEAU (m³/s) |
  | `Simulated Tw (degree C)` | Water temperature simulated by CEQUEAU (°C) |

- **Publication**:

### 1.3 Air temperature

#### Weather stations — St. Jean

- **File(s)**:
  - `R/00_raw_data/Air_Temperature_Stations/StJean_Air_Temperature_Stations_Records_1_to_10000.csv`
  - `R/00_raw_data/Air_Temperature_Stations/StJean_Air_Temperature_Stations_Records_10001_to_20000.csv`
  - `R/00_raw_data/Air_Temperature_Stations/StJean_Air_Temperature_Stations_Records_20001_to_30000.csv`
  - `R/00_raw_data/Air_Temperature_Stations/StJean_Air_Temperature_Stations_Records_30001_to_37060.csv`
- **Description**: daily air temperatures from several climate stations near the St. Jean River.
- **Source**: ClimateData.ca (https://climatedata.ca/)
- **Publication**: not published (downloadable from the source)

#### Weather stations — Trinité

- **File(s)**:
  - `R/00_raw_data/Air_Temperature_Stations/Trinite_Air_Temperature_Stations_Records_1_to_10000.csv`
  - `R/00_raw_data/Air_Temperature_Stations/Trinite_Air_Temperature_Stations_Records_10001_to_20000.csv`
  - `R/00_raw_data/Air_Temperature_Stations/Trinite_Air_Temperature_Stations_Records_20001_to_30000.csv`
  - `R/00_raw_data/Air_Temperature_Stations/Trinite_Air_Temperature_Stations_Records_30001_to_31306.csv`
- **Description**: daily air temperatures from several climate stations near the Trinité River.
- **Source**: ClimateData.ca (https://climatedata.ca/)
- **Publication**: not published (downloadable from the source).


#### Air temperature rasters (ERA5) — St. Jean

- **File(s)**: `R/00_raw_data/Air_Temperature_Raster/sj_1978_2024.zip`
- **Description**: hourly air temperature rasters (ERA5 reanalysis, GRIB format compressed as ZIP) in and around the St. Jean River watershed (0.25° grid, 2 x 6 cells; extent: -65.73 to -64.22 °E, 48.52 to 49.02 °N), from 1978-01-01 to 2024-12-31.
- **Source**: Copernicus Climate Data Store, ERA5 hourly data on single levels (https://cds.climate.copernicus.eu/ ; DOI: 10.24381/cds.adbb2d47)
- **Publication**: not published (downloadable from the source).

#### Air temperature rasters (ERA5) — Trinité

- **File(s)**:
  - `R/00_raw_data/Air_Temperature_Raster/tr_1975_2008.zip`
  - `R/00_raw_data/Air_Temperature_Raster/tr_2009_2025.zip`
- **Description**: hourly air temperature rasters (ERA5 reanalysis, GRIB format compressed as ZIP) in and around the Trinité River watershed (0.25° grid, 3 x 4 cells; extent: -68.19 to -67.19 °E, 49.20 to 49.95 °N), from 1975-01-01 to 2025-04-30.
- **Source**: Copernicus Climate Data Store, ERA5 hourly data on single levels (https://cds.climate.copernicus.eu/ ; DOI: 10.24381/cds.adbb2d47)
- **Publication**: not published (downloadable from the source).

### 1.4 Forestry activity

#### Ecoforest map with disturbances, map sheet 22A — St. Jean

- **File(s)**: `R/00_raw_data/Forest/CARTE_ECO_MAJ_22A_GPKG/CARTE_ECO_MAJ_22A.gpkg` (layer `pee_maj_22a`)
- **Description**: ecoforest map with disturbances, used to extract the area of disturbances in the St. Jean River watershed.
- **Source**: download the gpkg.zip file of map sheet 22A at https://diffusion.mffp.gouv.qc.ca/Diffusion/DonneeGratuite/Foret/DONNEES_FOR_ECO_SUD/Cartes_ecoforestieres_perturbations/02-Donnees/Decoupage250K/
- **Publication**: Downloadable from the source (too large for the repository).

#### Ecoforest map with disturbances, map sheet 22G — Trinité

- **File(s)**: `R/00_raw_data/Forest/CARTE_ECO_MAJ_22G_GPKG/CARTE_ECO_MAJ_22G.gpkg` (layer `pee_maj_22g`)
- **Description**: ecoforest map with disturbances, used to extract the area of disturbances in the Trinité River watershed.
- **Source**: download the gpkg.zip file of map sheet 22G at https://diffusion.mffp.gouv.qc.ca/Diffusion/DonneeGratuite/Foret/DONNEES_FOR_ECO_SUD/Cartes_ecoforestieres_perturbations/02-Donnees/Decoupage250K/
- **Publication**: Downloadable from the source (too large for the repository).

### 1.5 Food availability at sea

#### Continuous Plankton Recorder (CPR) data — Northwest Atlantic

- **File(s)**: `R/00_raw_data/Plankton_LabradorSea/NWAtl_CPR_Data.xlsx`
- **Description**: abundance of zooplankton taxa per CPR sample in the Northwest Atlantic, 1958–2021.
- **Columns used**:

  | Column | Description |
  |---|---|
  | `Sample_Id` | Sample identifier |
  | `Latitude`, `Longitude` | Sample position |
  | `Midpoint_Date_Local` | Sample date |
  | `Year`, `Month`, `Chlorophyll_Index` | Year, month and chlorophyll index (kept, not used in the calculations) |
  | `1`, `3`, `10`, `40`, `41`, `42`, `44`, `55`, `56`, `82`, `88`, `314`, `10742` | Abundance per taxon (CPR codes): *Calanus* I–IV, *Para-Pseudocalanus* spp., *Oithona* spp., *C. finmarchicus*, *C. helgolandicus*, *C. glacialis*, *C. hyperboreus*, *Metridia lucens*, *M. longa*, Hyperiidea, Euphausiacea, *Metridia* I–IV, *Pseudocalanus* spp. adults |

- **Source**: Johns, D. (Marine Biological Association) (2023). Selected NW Atlantic CPR taxa. The Archive for Marine Species and Habitats Data (DASSH). https://doi.org/10.17031/6513f104926b6
- **Publication**: Not published. To reproduce the analysis, download the dataset from the DOI above and save it as `R/00_raw_data/Plankton_LabradorSea/NWAtl_CPR_Data.xlsx`.


### 1.6 Acoustic telemetry

#### Post-smolt detections at acoustic receivers

- **File(s)**: `R/00_raw_data/telemetry/master_dataset_smolt_manuscript_dec_10_2025 QC.csv`
- **Description**: detections of tagged post-smolts at the acoustic receivers of the Gulf of St. Lawrence, the Strait of Belle Isle, St. Anthony, southern Labrador and Newfoundland. The fish come from the St. Jean, Restigouche and Cascapédia rivers (grouped for the St. Jean) and from the Trinité and Saguenay rivers (grouped for the Trinité). Semicolon-separated CSV file, with the comma as decimal separator.
- **Columns used**: `serialno_tag`, `Date_Hit`, `Location`, `river`, `deploy_year`, `Station_Name`, `SerialNo_Rec`, `Latitude`, `Longitude`, `fl_mm`
- **Source**: Lilly, J., Burke, C., Dufresne, C., Lavoie, D., April, J., Carr, J., ... Robertson, M. J. (2026). Atlantic Salmon post-smolt migration from Canadian rivers to the Labrador Sea: routes, timing, and influence of sea surface temperature. Canadian Journal of Fisheries and Aquatic Sciences, 83, 1-31. https://doi.org/10.1139/cjfas-2025-0197
- **Publication**: Not published — confidential data, excluded via `.gitignore`.
- **Example file**: `R/00_raw_data/telemetry/master_dataset_smolt_manuscript_EXAMPLE.csv` — a few rows in the same format as the original file (same columns, `;` separator, decimal comma), to show the structure expected by `01_MigrationDates_Determination.R`. Published.

### 1.7 Sea surface temperature

#### Sea surface temperature and sea ice fraction (ESA CCI SST, level 4)

- **File(s)**: `R/00_raw_data/SST_ORAS5_Labrador/ESACCI_SST_L4_<YYYY>_<MM>.nc` (one file per month)
- **Description**: daily sea surface temperature and sea ice fraction (ESA CCI SST, level 4, version 3.0, combined product, 0.05° grid) for the Labrador Sea and Gulf of St. Lawrence region (43 to 65 °N, -68 to -38 °E). Each monthly file contains the daily NetCDF files of the month.
- **Source**: Copernicus Climate Data Store, "satellite-sea-surface-temperature" dataset (https://cds.climate.copernicus.eu/)
- **Publication**: Downloadable from the source, with the script `05_Download_ORAS5_LabradorSea.R` (too large for the repository).

### 1.8 Fisheries and abundance at sea (ICES WGNAS)

#### Catches of North American salmon in West Greenland

- **File(s)**: `R/00_raw_data/Exploitation_Sea/WGNAS_2025_Table5_2_2_3_EstimatedPercentageNumbersSalmonCaughtGreenland.xlsx`
- **Description**: estimated number of salmon of North American origin caught each year in West Greenland (Table 5.2.2.3 of the WGNAS 2025 report, transcribed into Excel).
- **Columns used**: `Year`, `N_NA`
- **Source**: ICES (2025). Working Group on North Atlantic Salmon (WGNAS). ICES Scientific Reports.
- **Publication**: Downloadable from the source (Table 5.2.2.3 of the report).

#### Pre-Fishery abundance (PFA) on 1 January — North American stock complex (NAC)

- **File(s)**: `R/00_raw_data/Exploitation_Sea/WGNAS_2025_Table4_3_6_1_EstimatesPFA_Jan1st_NAC.xlsx`
- **Description**: Pre-Fishery abundance (PFA) on 1 January of the North American stock complex, in thousands of fish (Table 4.3.6.1 of the WGNAS 2025 report, transcribed into Excel).
- **Columns used**: `Year`, `Est_PFA1SWcohort` (1SW cohort, mature and non-mature), `Est_PFA1SWnmat` (non-mature 1SW)
- **Source**: ICES (2025). Working Group on North Atlantic Salmon (WGNAS). ICES Scientific Reports.
- **Publication**: Downloadable from the source (Table 4.3.6.1 of the report).

#### Pre-Fishery abundance (PFA) on 1 January — Southern European stock complex (SNEAC)

- **File(s)**: `R/00_raw_data/Exploitation_Sea/WGNAS_2025_Table3_3_5_3_EstimatesPFA_Jan1st_SNEAC.xlsx`
- **Description**: Pre-Fishery abundance (PFA) on 1 January of the non-mature fish of the Southern European stock complex, in thousands of fish (Table 3.3.5.3 of the WGNAS 2025 report, transcribed into Excel).
- **Columns used**: `Year`, `SNEAC_median`
- **Source**: ICES (2025). Working Group on North Atlantic Salmon (WGNAS). ICES Scientific Reports.
- **Publication**: Downloadable from the source (Table 3.3.5.3 of the report).

#### Exploitation rates of North American salmon (Imlay et al., 2024)

- **File(s)**: `R/00_raw_data/Exploitation_Sea/FisheryExploitationRates.rds`
- **Description**: annual exploitation rates of North American one-sea-winter (1SW) and two-sea-winter (2SW) salmon.
- **Columns used**: `YEAR`, `er1NAC`, `er2NAC`
- **Source**: data associated with Imlay, T. L., et al. (2024). Body length changes for Atlantic salmon (*Salmo salar*) over five decades exhibit weak spatial synchrony over a broad latitudinal gradient. *Ecology and Evolution*, 14(6), e11538.
- **Publication**: Downloadable from the source (file `FisheryExploitationRates.rds` of the OSF project).

---

## 2. Data preparation and variable construction

Folder: `R/02_prepare_data/`

Data cleaning and variable construction scripts. Each script is described with the same template: description, inputs, outputs (with their publication status) and, if applicable, variables created. Intermediate data are listed in the outputs of the script that produces them.


### 2.1 Shared data

Folder: `R/02_prepare_data/00_SharedData/`

Values or data used by several scripts to build the freshwater or marine variables.

#### `Smolt_Proportions.R`

- **Description**: computes the proportions of smolts by age class, at the cohort level, for the St. Jean and Trinité rivers. Missing age classes are filled with the river mean when their mean contribution is below 10%; otherwise, the cohort is excluded. The final proportions are renormalized so that they sum to 1.
- **Inputs**: `R/00_raw_data/salmon/Salmon_Monitoring_StJean_Trinite_2024.csv`
- **Outputs**: `R/01_derived_data/Smolt_Proportions.fst` — Not published (regenerated by this script)

#### `Growth_Period_Determination.R`

- **Description**: determines, for each river and each year, the growth period, i.e., the period during which the water temperature predicted by a GAM exceeds 6 °C. The cyclic GAM (`bs = "cc"`, `k = 30`, REML) is fitted year by year; the period runs from the first upward crossing to the last downward crossing of the 6 °C threshold by the predicted curve.
- **Inputs**:
  - St. Jean: `R/00_raw_data/Temperature_Discharge_River/CEQEAU_Temperature_Discharge_1979_2024_SJ.csv`
  - Trinité: `R/00_raw_data/Temperature_Discharge_River/CEQEAU_Temperature_Discharge_1979_2024_TR.csv`
- **Outputs**:
  - St. Jean: `R/01_derived_data/Growth_Periods_StJean.fst` — Not published (regenerated by this script)
  - Trinité: `R/01_derived_data/Growth_Periods_Trinite.fst` — Not published (regenerated by this script)

#### `Smoltification_Date_Determination.R`

- **Description**: for each year, computes the median, the standard deviation and the 5, 25, 75 and 95% quantiles of the smolt capture day (weighted by the number of captures). Years without data are imputed by pooling the captures of the years within a ±5-year window.
- **Inputs**:
  - St. Jean: `R/00_raw_data/salmon/Smolt_Daily_Captures_StJean.csv`
  - Trinité: `R/00_raw_data/salmon/Smolt_Daily_Captures_Trinite.csv`
- **Outputs**:
  - St. Jean: `R/01_derived_data/Year_Median_Date_StJean.fst` — Published
  - Trinité: `R/01_derived_data/Year_Median_Date_Trinite.fst` — Published

### 2.2 Freshwater variables

#### 2.2.1 Discharge

Folder: `R/02_prepare_data/RIVER_Discharge/`

##### `River_Discharge.R`

- **Description**: creates the discharge variables of the St. Jean and Trinité rivers during two periods of interest: emergence (Julian days 137 to 154) and growth (annual periods produced by `Growth_Period_Determination.R`). The growth-period variables are weighted at the cohort level with the proportions of smolts by age class (`Smolt_Proportions.fst`).
- **Inputs**:
  - `R/00_raw_data/Temperature_Discharge_River/CEQEAU_Temperature_Discharge_1979_2024_SJ.csv`
  - `R/00_raw_data/Temperature_Discharge_River/CEQEAU_Temperature_Discharge_1979_2024_TR.csv`
  - `R/01_derived_data/Growth_Periods_StJean.fst`
  - `R/01_derived_data/Growth_Periods_Trinite.fst`
  - `R/01_derived_data/Smolt_Proportions.fst`
- **Outputs**:
  - `R/03_clean_data/River/River_Discharge.fst` — Published.
- **Variables created**:

  | Variable | Period | Description |
  |---|---|---|
  | `co_gr_mean_disc` | Growth (cohort-weighted) | Mean daily discharge |
  | `co_gr_median_disc` | Growth (cohort-weighted) | Median daily discharge |
  | `co_gr_var_disc` | Growth (cohort-weighted) | Variance of daily discharge |
  | `co_gr_mean_75qt_disc` | Growth (cohort-weighted) | 75th percentile of daily discharge |
  | `co_30d_mean_disch` | Growth (cohort-weighted) | Minimum of the 30-day moving average of daily discharge |
  | `em_mean_disc` | Emergence | Mean daily discharge |
  | `em_median_disc` | Emergence | Median daily discharge |
  | `em_var_disc` | Emergence | Variance of daily discharge |
  | `em_quant_75qt_disc` | Emergence | 75th percentile of daily discharge |


#### 2.2.2 Forestry activity

Folder: `R/02_prepare_data/RIVER_Forest_Activities/`

##### `01_Watershed_Delim_STJ.R`

- **Description**: downloads a digital elevation model (DEM) and the river course, delineates the St. Jean River watershed with WhiteboxTools, then exports it as a polygon.
- **Inputs**:
  - No file. The DEM (`elevatr`, `z = 12`) and the river course (OpenStreetMap, via `osmdata`) are downloaded by the script, within a bounding box defined in the script.
- **Outputs**:
  - `R/01_derived_data/ForestActivities/poly_watershed_STJ.shp` — Not published (regenerated by this script)
  - WhiteboxTools working files, not used by the other scripts: `dem_utm_STJ`, `dem_filled_STJ`, `d8_pointer_STJ`, `flow_accum_STJ` and `raster_watershed_STJ` (.tif); `outlet_STJ` and `snapped_outlet_STJ` (.shp) — Not published (regenerated by this script)
- **Requirements**: WhiteboxTools must be installed once with `whitebox::install_whitebox()`.

##### `02_Watershed_Delim_TRI.R`

- **Description**: downloads a DEM and the river course, delineates the Trinité River watershed with WhiteboxTools, then exports it as a polygon. Unlike the St. Jean script, the points of the river course at an elevation of 5 m or less (`elev_min`) are excluded before the outlet is chosen, to avoid the tidal zone.
- **Inputs**:
  - No file. The DEM (`elevatr`, `z = 12`) and the river course (OpenStreetMap, via `osmdata`) are downloaded by the script, within a bounding box defined in the script.
- **Outputs**:
  - `R/01_derived_data/ForestActivities/poly_watershed_TRI.shp` — Not published (regenerated by this script)
  - WhiteboxTools working files, not used by the other scripts: `dem_utm_TRI`, `dem_filled_TRI`, `d8_pointer_TRI`, `flow_accum_TRI` and `raster_watershed_TRI` (.tif); `outlet_TRI` and `snapped_outlet_TRI` (.shp) — Not published (regenerated by this script)
- **Requirements**: WhiteboxTools must be installed once with `whitebox::install_whitebox()`.

##### `03_Extract_Forest_Activities_STJ.R`

- **Description**: extracts the forest disturbances of the St. Jean River watershed. The script keeps stand-replacing and partial disturbances that occurred since 1980, excludes codes P, PRR, REA, RIA and RPS, clips the polygons to the watershed boundaries and computes the area of each one.
- **Inputs**:
  - `R/00_raw_data/Forest/CARTE_ECO_MAJ_22A_GPKG/CARTE_ECO_MAJ_22A.gpkg` (layer `pee_maj_22a`)
  - `R/01_derived_data/ForestActivities/poly_watershed_STJ.shp`
- **Outputs**:
  - `R/01_derived_data/ForestActivities/ForestActivity_STJ.fst` — Not published (regenerated by this script). One row per disturbed polygon, with the disturbance type, the disturbance year and its area.

##### `04_Extract_Forest_Activities_TRI.R`

- **Description**: same processing as `03_Extract_Forest_Activities_STJ.R`, applied to the Trinité River watershed.
- **Inputs**:
  - `R/00_raw_data/Forest/CARTE_ECO_MAJ_22G_GPKG/CARTE_ECO_MAJ_22G.gpkg` (layer `pee_maj_22g`)
  - `R/01_derived_data/ForestActivities/poly_watershed_TRI.shp`
- **Outputs**:
  - `R/01_derived_data/ForestActivities/ForestActivity_TRI.fst` — Not published (regenerated by this script). One row per disturbed polygon, with the disturbance type, the disturbance year and its area.

##### `05_Weighted_Cohort_ForestActivities.R`

- **Description**: computes, for each salmon cohort of the St. Jean and Trinité rivers, the mean annual area of forest disturbances, weighted by the proportions of smolts by age class.
- **Inputs**:
  - `R/01_derived_data/ForestActivities/ForestActivity_STJ.fst`
  - `R/01_derived_data/ForestActivities/ForestActivity_TRI.fst`
  - `R/01_derived_data/Smolt_Proportions.fst`
- **Outputs**:
  - `R/03_clean_data/River/ForestActivity_Cohort.fst` — Published.
- **Variables created**:

  | Variable | Period | Description |
  |---|---|---|
  | `co_mean_forest` | Years in the river of the cohort (cohort-weighted) | Mean annual area of forest disturbances (km²), weighted by the proportions of smolts by age class |

#### 2.2.3 Temperature

Folder: `R/02_prepare_data/RIVER_Temperature/`

##### `River_Temperature.R`

- **Description**: creates the water temperature variables of the St. Jean and Trinité rivers during three periods: emergence (Julian days 137 to 154), growth (annual periods produced by `Growth_Period_Determination.R`) and the full year. The growth and full-year variables are weighted at the cohort level with the proportions of smolts by age class (`Smolt_Proportions.fst`); the emergence variables refer to the first year in the river and are not weighted. The degree-day thresholds (16 and 20 °C) come from Rooke et al. (2019); the potential growth thermal index (PGTI) follows Ouellet-Proulx et al. (2023), with T<sub>min</sub> = 6 °C, T<sub>opt</sub> = 18 °C and T<sub>max</sub> = 27 °C.
- **Inputs**:
  - `R/00_raw_data/Temperature_Discharge_River/CEQEAU_Temperature_Discharge_1979_2024_SJ.csv`
  - `R/00_raw_data/Temperature_Discharge_River/CEQEAU_Temperature_Discharge_1979_2024_TR.csv`
  - `R/01_derived_data/Growth_Periods_StJean.fst`
  - `R/01_derived_data/Growth_Periods_Trinite.fst`
  - `R/01_derived_data/Smolt_Proportions.fst`
- **Outputs**:
  - `R/03_clean_data/River/River_Temperatures.fst` — Published.
- **Variables created**:

  | Variable | Period | Description |
  |---|---|---|
  | `co_gr_mean_temp` | Growth (cohort-weighted) | Mean daily temperature |
  | `co_gr_median_temp` | Growth (cohort-weighted) | Median daily temperature |
  | `co_gr_var_temp` | Growth (cohort-weighted) | Variance of daily temperature |
  | `co_gr_dd_u_16_temp` | Growth (cohort-weighted) | Sum of daily temperatures ≤ 16 °C |
  | `co_gr_dd_o_20_temp` | Growth (cohort-weighted) | Sum of daily temperatures ≥ 20 °C |
  | `co_year_mean_sum_ptgi` | Full year (cohort-weighted) | Annual sum of the PGTI |
  | `co_year_mean_sum_inv_ptgi` | Full year (cohort-weighted) | Annual sum of 1 − PGTI |
  | `em_mean_temp` | Emergence | Mean daily temperature |
  | `em_median_temp` | Emergence | Median daily temperature |
  | `em_var_temp` | Emergence | Variance of daily temperature |
  | `em_dd_u_16_temp` | Emergence | Sum of daily temperatures ≤ 16 °C |
  | `em_dd_o_20_temp` | Emergence | Sum of daily temperatures ≥ 20 °C |
  | `em_sum_ptgi` | Emergence | Sum of the PGTI |
  | `em_sum_inv_pgti` | Emergence | Sum of 1 − PGTI |


#### 2.2.4 Winter temperature (freezing degree-days)

Folder: `R/02_prepare_data/RIVER_Freezing_Degree_Day/`

##### `01_Winter_Temperature_Extract_Stations_Temperature.R`

- **Description**: imports, merges and cleans the daily air temperatures of several climate stations near the St. Jean and Trinité rivers, checks the correlation between stations, then computes for each river the daily mean of the available stations.
- **Inputs**:
  - St. Jean: `R/00_raw_data/Air_Temperature_Stations/StJean_Air_Temperature_Stations_Records_*.csv` (4 files)
  - Trinité: `R/00_raw_data/Air_Temperature_Stations/Trinite_Air_Temperature_Stations_Records_*.csv` (4 files)
- **Outputs**:
  - St. Jean: `R/01_derived_data/Winter_Temperature/sj_air_temp_daily_mean_stations.fst` — Not published (regenerated by this script)
  - Trinité: `R/01_derived_data/Winter_Temperature/tr_air_temp_daily_mean_stations.fst` — Not published (regenerated by this script)
  - Columns: `LOCAL_DATE`, `Station_Mean_T`, `Station_Min_T`, `Station_Max_T` (station means, °C)

##### `02_Winter_Temperature_Extract_Raster_Temperature.R`

- **Description**: imports the ERA5 hourly air temperature rasters and computes, for each river, the daily mean temperature, defined as the mean of the minimum and maximum temperatures of the day (same convention as Canadian climate stations).
- **Inputs**:
  - St. Jean: `R/00_raw_data/Air_Temperature_Raster/sj_1978_2024.zip`
  - Trinité: `R/00_raw_data/Air_Temperature_Raster/tr_1975_2008.zip`, `R/00_raw_data/Air_Temperature_Raster/tr_2009_2025.zip`
- **Outputs**:
  - St. Jean: `R/01_derived_data/Winter_Temperature/sj_air_temp_daily_mean_rasters.fst` — Not published (regenerated by this script)
  - Trinité: `R/01_derived_data/Winter_Temperature/tr_air_temp_daily_mean_rasters.fst` — Not published (regenerated by this script)
  - Columns: `date`, `Raster_Mean_T` (°C)

##### `03_Winter_Temperature_Merge_Station_Raster.R`

- **Description**: merges the daily temperatures from the stations and the ERA5 rasters, compares the two sources (correlation, absolute difference, Wilcoxon test), then produces a combined series from October to April. The temperature of a given day is the mean of the two sources when both are available, otherwise the value of the only available source.
- **Inputs**:
  - `R/01_derived_data/Winter_Temperature/sj_air_temp_daily_mean_stations.fst`
  - `R/01_derived_data/Winter_Temperature/tr_air_temp_daily_mean_stations.fst`
  - `R/01_derived_data/Winter_Temperature/sj_air_temp_daily_mean_rasters.fst`
  - `R/01_derived_data/Winter_Temperature/tr_air_temp_daily_mean_rasters.fst`
- **Outputs**:
  - St. Jean: `R/01_derived_data/Winter_Temperature/sj_air_temp_daily_mean_combined.fst` — Not published (regenerated by this script)
  - Trinité: `R/01_derived_data/Winter_Temperature/tr_air_temp_daily_mean_combined.fst` — Not published (regenerated by this script)
  - Columns: `date`, `Mean_T` (°C)

##### `04_Winter_Temperature_Bisallon_Index.R`

- **Description**: computes freezing degree-day indices (FDD, sum of the absolute values of negative daily mean temperatures) over several winter periods, following the approach of Bisaillon (2005). Each index is weighted at the cohort level: the winters of each year in the river are weighted by the proportion of the cohort still present in the river at that age (`Smolt_Proportions.fst`).
- **Inputs**:
  - `R/01_derived_data/Winter_Temperature/sj_air_temp_daily_mean_combined.fst`
  - `R/01_derived_data/Winter_Temperature/tr_air_temp_daily_mean_combined.fst`
  - `R/01_derived_data/Smolt_Proportions.fst`
- **Outputs**:
  - `R/03_clean_data/River/River_Winter_Temperatures.fst` — Published.

- **Variables created** (all cohort-weighted):

  | Variable | Period | Description |
  |---|---|---|
  | `co_fdd_per_oc_no` | 1 October to 30 November | Cumulative freezing degree-days |
  | `co_fdd_per_oc_de` | 1 October to 31 December | Cumulative freezing degree-days |
  | `co_fdd_per_oc_ja` | 1 October to 31 January | Cumulative freezing degree-days |
  | `co_fdd_per_oc_fe` | 1 October to the end of February | Cumulative freezing degree-days |
  | `co_fdd_per_oc_ma` | 1 October to 31 March | Cumulative freezing degree-days |
  | `co_fdd_per_oc_ap` | 1 October to 30 April | Cumulative freezing degree-days |
  | `co_fdd_mon_oc` | October | Freezing degree-days |
  | `co_fdd_mon_no` | November | Freezing degree-days |
  | `co_fdd_mon_de` | December | Freezing degree-days |
  | `co_fdd_mon_ja` | January | Freezing degree-days |
  | `co_fdd_mon_fe` | February | Freezing degree-days |
  | `co_fdd_mon_ma` | March | Freezing degree-days |
  | `co_fdd_mon_ap` | April | Freezing degree-days |
  | `co_fdd_bimon_no_de` | November and December | Freezing degree-days |
  | `co_fdd_bimon_ja_fe` | January and February | Freezing degree-days |
  | `co_fdd_bimon_ma_ap` | March and April | Freezing degree-days |



### 2.3 Marine variables

#### 2.3.1 Thermal conditions

Folder: `R/02_prepare_data/SEA_SST/`

The seven scripts in this folder are run in order. Scripts 01 to 04 define the migration windows and the polygons of the migration route; scripts 05 and 06 download the sea surface temperatures and derive a daily thermal habitat index per polygon from them; script 07 sums this index over the migration windows of each cohort.

##### `01_MigrationDates_Determination.R`

- **Description**: estimates the dates of migration at sea from the acoustic telemetry data. The date of entry into the Strait of Belle Isle is the median of the day of the first detection of each fish at the receivers of the strait (or, failing that, of the gulf). The date of exit from the strait, which is also the date of entry into the Labrador Sea, is the median of the day of the last detection at the receivers of southern Labrador, Newfoundland or St. Anthony (two river groups pooled). The mean swimming speed is computed from the least-cost distances between the entry and exit receivers, measured over water only; it is used to estimate the time needed to cross the Labrador Sea entrance polygon (300 km). The date of exit from the Labrador Sea is set to 1 May, for both 1SW and 2SW (Imlay et al., 2024).
- **Inputs**:
  - `R/00_raw_data/telemetry/master_dataset_smolt_manuscript_dec_10_2025 QC.csv` (not published; see the example file, section 1.6)
  - `R/01_derived_data/Year_Median_Date_StJean.fst`
  - `R/01_derived_data/Year_Median_Date_Trinite.fst`
- **Outputs**:
  - `R/01_derived_data/Year_Migration_Dates_River_Sea.fst` — Published. Columns: `river`, `cohort`, `smolting_doy`, `smolt_quant_25_doy`, `smolt_quant_75_doy`, `sobi_entry`, `sobi_exit`, `entering_lab_sea_entry`, `entering_lab_sea_exit`, `lab_sea_entry`, `lab_sea_exit_1sw`, `lab_sea_exit_2sw`.

##### `02_MigrationRoutes_Creation.R`

- **Description**: builds the five polygons of the post-smolt migration route, from the river mouths to the Labrador Sea (projection NAD83(CSRS) / Statistics Canada Lambert, EPSG:6623):
  - **Trinité corridor**: 100 km wide coastal band along the North Shore, from 5 km downstream of the river mouth to the boundary of the Strait of Belle Isle, cut south of Anticosti.
  - **St. Jean corridor**: 100 km wide buffer around a line that starts in Gaspé Bay, passes about 50 km east of Anticosti and joins the middle of the boundary of the Strait of Belle Isle.
  - **Strait of Belle Isle**: quadrilateral between the southern boundary of the strait and two anchor points to the north, on the coasts of Labrador and Newfoundland.
  - **Labrador Sea entrance**: square of about 300 km × 300 km located just north of the strait.
  - **Labrador Sea**: offshore area between 45 and 61 °N and between 40 and 64 °W, excluding land and the strait.
- **Inputs**:
  - No file. The land boundaries (Natural Earth) are downloaded by the `rnaturalearth` package.
- **Outputs**:
  - `R/00_raw_data/ReferenceMigrationLayers/TriniteMigrationCorridor.shp` — Not published (regenerated by this script)
  - `R/00_raw_data/ReferenceMigrationLayers/StJeanMigrationCorridor.shp` — Not published (regenerated by this script)
  - `R/00_raw_data/ReferenceMigrationLayers/BelleIsleStrait.shp` — Not published (regenerated by this script)
  - `R/00_raw_data/ReferenceMigrationLayers/LabradorSeaEntrance.shp` — Not published (regenerated by this script)
  - `R/00_raw_data/ReferenceMigrationLayers/LabradorSea.shp` — Not published (regenerated by this script; also used by `02_Plankton_Models.R`)

##### `03_MigrationRoute_Division_St_Jean.R`

- **Description**: divides the St. Jean migration corridor into 5 successive segments of equal length along the coast. The cutting lines connect points evenly spaced along the coastal side of the corridor to the nearest points on its offshore side. These segments make it possible to follow the progression of the fish through the corridor (see `07_Subset_By_Periods.R`).
- **Inputs**:
  - `R/00_raw_data/ReferenceMigrationLayers/StJeanMigrationCorridor.shp`
- **Outputs**:
  - `R/00_raw_data/ReferenceMigrationLayers/StJeanMigrationCorridor_Parts.shp` — Not published (regenerated by this script). One row per segment (`FID` = 1 to 5).

##### `04_MigrationRoute_Division_Trinite.R`

- **Description**: same processing as `03_MigrationRoute_Division_St_Jean.R`, applied to the Trinité corridor. The coastal and offshore sides are distinguished by their distance to the coast (< 500 m and > 5 km).
- **Inputs**:
  - `R/00_raw_data/ReferenceMigrationLayers/TriniteMigrationCorridor.shp`
- **Outputs**:
  - `R/00_raw_data/ReferenceMigrationLayers/TriniteMigrationCorridor_Parts.shp` — Not published (regenerated by this script). One row per segment (`FID` = 1 to 5).

##### `05_Download_ORAS5_LabradorSea.R`

- **Description**: downloads from the Copernicus Climate Data Store (CDS API, `ecmwfr` package) the daily level-4 ESA CCI sea surface temperatures for the Labrador Sea region, one request per month. Months already downloaded are skipped, and the script reports the missing months at the end.
- **Inputs**:
  - No file. The data are downloaded by the script.
- **Outputs**:
  - `R/00_raw_data/SST_ORAS5_Labrador/ESACCI_SST_L4_<YYYY>_<MM>.nc` — Downloadable from the source (see section 1.7)
- **Requirements**: an account on the Copernicus Climate Data Store and an API token, to be placed in a `.Renviron` file at the root of the repository (`CDS_TOKEN=your_token`). This file is excluded via `.gitignore`.

##### `06_SST_Index.R`

- **Description**: computes a daily thermal habitat index for each polygon of the migration route and each segment of the corridors. For each pixel, the temperature (converted from kelvins to °C) gives a growth index between 0 and 1 according to a piecewise linear function: 0 below T<sub>min</sub>, linear increase up to T<sub>opt1</sub>, plateau at 1 between T<sub>opt1</sub> and T<sub>opt2</sub>, linear decrease down to T<sub>max</sub>, then 0 beyond. The parameters are T<sub>min</sub> = 1 °C, T<sub>opt</sub> = 6–14 °C and T<sub>max</sub> = 18 °C for 1SW, and T<sub>min</sub> = 1 °C, T<sub>opt</sub> = 6–12 °C and T<sub>max</sub> = 15 °C for 2SW. Pixels with 50% or more ice cover are given 0. The index of each pixel is multiplied by the cell area (km²), then summed over the polygon.
- **Inputs**:
  - `R/00_raw_data/SST_ORAS5_Labrador/ESACCI_SST_L4_<YYYY>_<MM>.nc`
  - `R/00_raw_data/ReferenceMigrationLayers/LabradorSea.shp`
  - `R/00_raw_data/ReferenceMigrationLayers/LabradorSeaEntrance.shp`
  - `R/00_raw_data/ReferenceMigrationLayers/BelleIsleStrait.shp`
  - `R/00_raw_data/ReferenceMigrationLayers/StJeanMigrationCorridor_Parts.shp`
  - `R/00_raw_data/ReferenceMigrationLayers/TriniteMigrationCorridor_Parts.shp`
- **Outputs**:
  - `R/01_derived_data/SST_Index/SST_Index_<Polygon>/TH_<YYYY>_<MM>.fst` (one file per polygon and per month) — Not published (regenerated by this script). One row per day and per segment. Columns: `date`, `year`, `month`, `doy`, `FID`, `n_pixels`, `TH_1SW_sum`, `TH_1SW_mean`, `TH_2SW_sum`, `TH_2SW_mean`, `sst_mean`, `sst_median`, `sst_sd`, `sst_min`, `sst_max`, `ice_mean`.

##### `07_Subset_By_Periods.R`

- **Description**: sums the daily thermal habitat index (`TH_*_sum`) over the migration window of each cohort, for each polygon. Missing dates are replaced by the mean of the 5 neighbouring cohorts on each side, and data after 31 December 2023 are excluded. The windows are the following:
  - **River corridor**: from the smoltification date to the entry into the Strait of Belle Isle. The window is divided into 5 equal sub-periods, one per corridor segment, travelled in order.
  - **Strait of Belle Isle**: from entry into to exit from the strait.
  - **Labrador Sea entrance**: from the exit from the strait to the end of the crossing of the polygon.
  - **Labrador Sea**: from the entry into the Labrador Sea (year t) to 1 May of year t + 1 for 1SW, or of year t + 2 for 2SW. For 2SW, the 1SW thermal parameters apply until 1 May of year t + 1, then the 2SW parameters.

  In the polygons located before the Labrador Sea, the 1SW parameters apply to both groups, since all fish are then in their first year at sea.
- **Inputs**:
  - `R/01_derived_data/SST_Index/SST_Index_<Polygon>/TH_<YYYY>_<MM>.fst`
  - `R/01_derived_data/Year_Migration_Dates_River_Sea.fst`
- **Outputs**:
  - `R/03_clean_data/Sea/TH_annual_stj_1sw.fst` — Published.
  - `R/03_clean_data/Sea/TH_annual_stj_2sw.fst` — Published.
  - `R/03_clean_data/Sea/TH_annual_tri_1sw.fst` — Published.
  - `R/03_clean_data/Sea/TH_annual_tri_2sw.fst` — Published.
- **Variables created** (one row per cohort):

  | Variable | Period | Description |
  |---|---|---|
  | `StJeanCorridor_TH` (St. Jean) / `TriniteCorridor_TH` (Trinité) | Smoltification → entry into the Strait of Belle Isle | Cumulative thermal habitat index in the river corridor (1SW parameters) |
  | `BelleIsle_TH` | Entry → exit from the Strait of Belle Isle | Cumulative thermal habitat index in the strait (1SW parameters) |
  | `LabradorEntrance_TH` | Exit from the strait → end of the crossing of the Labrador Sea entrance | Cumulative thermal habitat index at the Labrador Sea entrance (1SW parameters) |
  | `LabradorSea_TH` | Entry into the Labrador Sea → 1 May of t + 1 (1SW files) or of t + 2 (2SW files) | Cumulative thermal habitat index in the Labrador Sea (1SW parameters, then 2SW after 1 May of t + 1 for the 2SW files) |
  

#### 2.3.2 Food availability

Folder: `R/02_prepare_data/SEA_Food_Availability/`

##### `01_Prepare_CPR_Data.R`

- **Description**: prepares the raw CPR data. The script renames the taxon columns (CPR codes), replaces presence-only values (0.001) with NA, converts the table to long format (one row per sample and per taxon), then converts abundance into dry biomass from the body length of each taxon (taken from the data of Imlay et al., 2024, available under a CC-BY 4.0 licence: https://doi.org/10.17605/OSF.IO/78PWT). The individual mass is computed with the formula of Richardson et al. (2006) (0.08 × L<sup>2.1</sup>), except for *Oithona* spp., for which the formula of Krylov (1968) is used (0.008 × L<sup>3</sup>).
- **Inputs**:
  - `R/00_raw_data/Plankton_LabradorSea/NWAtl_CPR_Data.xlsx`
- **Outputs**:
  - `R/01_derived_data/Food_Availability/cpr_data_raw.fst` — Not published (regenerated by this script). Columns: `Sample_Id`, `Latitude`, `Longitude`, `Midpoint_Date_Local`, `Year`, `Month`, `Chlorophyll_Index`, `taxa`, `count`, `length_mm`, `biomass_g`.

##### `02_Plankton_Models.R`

- **Description**: produces an annual zooplankton biomass index for the Labrador Sea. The script sums the biomass of all taxa per sample, keeps the samples located within the Labrador Sea polygon, then computes the mean monthly biomass (log(mean + 1) scale). A Gaussian GAM (`s(YEAR, k = 10) + s(Month, bs = "cc", k = 12) + te(mLat, mLon)`) is fitted to these monthly means; it was selected after comparison with a Tweedie model (DHARMa diagnostics and temporal trend of the residuals). The model then predicts the biomass per month and per year at 100 points drawn at random among the sampling positions. The 12 monthly values are summed for each year and each point, and the annual index is the median of the 100 points. The index is validated by leave-one-year-out cross-validation (LOO). Because consecutive years are strongly correlated, only the year of sea entry (t) is kept. The index is the same for both rivers.
- **Inputs**:
  - `R/01_derived_data/Food_Availability/cpr_data_raw.fst`
  - `R/00_raw_data/ReferenceMigrationLayers/LabradorSea.shp`
- **Outputs**:
  - `R/03_clean_data/Sea/Plankton_LabradorSea_stj_1SW.fst` — Published.
  - `R/03_clean_data/Sea/Plankton_LabradorSea_tri_1SW.fst` — Published.
  - `R/03_clean_data/Sea/Plankton_LabradorSea_stj_2SW.fst` — Published.
  - `R/03_clean_data/Sea/Plankton_LabradorSea_tri_2SW.fst` — Published.
- **Variables created**:

  | Variable | Period | Description |
  |---|---|---|
  | `Plankton_t` | Year t (1SW and 2SW) | Annual zooplankton biomass index in the Labrador Sea (sum of the 12 predicted monthly values, median of the 100 points) |
  | `Plankton_tp1` | Year t + 1 (2SW only) | Same index, for the following year |

#### 2.3.3 Temperature difference between sea and river at sea entry

Folder: `R/02_prepare_data/SEA_Temperature_River_VS_Sea/`

##### `Temperature_River_VS_Sea.R`

- **Description**: computes, for each year and each river, the difference between the sea surface temperature and the river water temperature during the smolt downstream migration period. This period runs from the 5th to the 95th percentile of the trap capture dates (`Year_Median_Date_*.fst`). The river temperature is the mean of the daily temperature simulated by CEQUEAU over this period. The sea temperature is the mean of the daily mean temperatures (ESA CCI SST) of the first segment of the river's migration corridor (`FID = 1`). The difference is computed as the sea temperature minus the river temperature: a positive value indicates a sea warmer than the river. Because this difference is measured in the year of downstream migration, it is the same for the 1SW and 2SW of a given cohort; the rows are duplicated to follow the structure of the survival tables.
- **Inputs**:
  - `R/00_raw_data/SST_ORAS5_Labrador/ESACCI_SST_L4_<YYYY>_<MM>.nc`
  - `R/00_raw_data/Temperature_Discharge_River/CEQEAU_Temperature_Discharge_1979_2024_SJ.csv`
  - `R/00_raw_data/Temperature_Discharge_River/CEQEAU_Temperature_Discharge_1979_2024_TR.csv`
  - `R/00_raw_data/ReferenceMigrationLayers/StJeanMigrationCorridor_Parts.shp`
  - `R/00_raw_data/ReferenceMigrationLayers/TriniteMigrationCorridor_Parts.shp`
  - `R/01_derived_data/Year_Median_Date_StJean.fst`
  - `R/01_derived_data/Year_Median_Date_Trinite.fst`
- **Outputs**:
  - `R/03_clean_data/Sea/DeltaT_SST_River_stj_1SW.fst` — Published.
  - `R/03_clean_data/Sea/DeltaT_SST_River_stj_2SW.fst` — Published.
  - `R/03_clean_data/Sea/DeltaT_SST_River_tri_1SW.fst` — Published.
  - `R/03_clean_data/Sea/DeltaT_SST_River_tri_2SW.fst` — Published.
- **Variables created** (one row per year):

  | Variable | Period | Description |
  |---|---|---|
  | `mean_sst_during_smolting` | Downstream migration (5th to 95th percentile of captures) | Mean sea surface temperature in the first segment of the corridor (°C) |
  | `mean_river_temp` | Downstream migration (5th to 95th percentile of captures) | Mean river water temperature (°C) |
  | `delta_T_sst_minus_river` | Downstream migration (5th to 95th percentile of captures) | Temperature difference between sea and river (°C) |
  | `n_days_window`, `n_days_valid` | — | Number of days in the period and number of days used for the sea temperature |


#### 2.3.4 Exploitation rates at sea

Folder: `R/02_prepare_data/SEA_Exploitation/`

##### `ExploitationRate_Greenland.R`

- **Description**: computes the exploitation rates at sea and aligns them with the cohort year (year of smolt downstream migration, t). The exploitation rate in Greenland is the number of North American salmon caught in Greenland divided by the Pre-Fishery abundance of non-mature 1SW of the North American stock complex, for 1982 to 2023. For 1993 and 1994, because the catches are not available in the table, the rates were read directly from Figure 5.1.3.1 of the WGNAS 2025 report. The 1SW and 2SW exploitation rates of the North American stock complex come from Imlay et al. (2024).
- **Inputs**:
  - `R/00_raw_data/Exploitation_Sea/WGNAS_2025_Table5_2_2_3_EstimatedPercentageNumbersSalmonCaughtGreenland.xlsx`
  - `R/00_raw_data/Exploitation_Sea/WGNAS_2025_Table4_3_6_1_EstimatesPFA_Jan1st_NAC.xlsx`
  - `R/00_raw_data/Exploitation_Sea/FisheryExploitationRates.rds`
- **Outputs**:
  - `R/03_clean_data/Sea/ExploitationRate_Sea.fst` — Published.
- **Variables created** (one row per cohort):

  | Variable | Period | Description |
  |---|---|---|
  | `er1NAC_lag1` | Year t + 1 | Exploitation rate of North American 1SW (Imlay et al., 2024) |
  | `er_WG_lag1` | Year t + 1 | Exploitation rate in Greenland (catches of North American salmon / PFA of non-mature 1SW) |
  | `er2NAC_lag2` | Year t + 2 | Exploitation rate of North American 2SW (Imlay et al., 2024) |

#### 2.3.5 Pre-Fishery abundance (density at sea)

Folder: `R/02_prepare_data/SEA_Density_Lab/`

##### `Density_Lab.R`

- **Description**: computes the Pre-Fishery abundance (PFA) on 1 January of the salmon present at sea, following the method of Imlay et al. (2024). The abundance during the first year at sea is the North American 1SW cohort (mature and non-mature), the only one present in the Northwest Atlantic at this stage. The abundance during the second year at sea is the sum of the North American non-mature fish and the Southern European non-mature fish, which then share the same areas (Labrador, Greenland, Iceland, Faroe Islands).
- **Inputs**:
  - `R/00_raw_data/Exploitation_Sea/WGNAS_2025_Table4_3_6_1_EstimatesPFA_Jan1st_NAC.xlsx`
  - `R/00_raw_data/Exploitation_Sea/WGNAS_2025_Table3_3_5_3_EstimatesPFA_Jan1st_SNEAC.xlsx`
- **Outputs**:
  - `R/03_clean_data/Sea/Density_Lab.fst` — Published.
- **Variables created** (one row per year):

  | Variable | Period | Description |
  |---|---|---|
  | `pfa_1sw` | First year at sea | Pre-Fishery abundance of the North American 1SW cohort (thousands of fish) |
  | `pfa_2sw` | Second year at sea | Pre-Fishery abundance of the North American and Southern European non-mature fish (thousands of fish) |

---

## 3. Modelling

Folder: `R/04_analysis/`

Scripts for dataset assembly, variable selection and modelling, described with the same template as section 2. The analyses are carried out separately for the freshwater phase (`River/`) and for the marine phase (`Sea/`).

### 3.1 Freshwater survival

Folder: `R/04_analysis/River/`

#### `01_Combine_Covariates_River.R`

- **Description**: assembles, for each river, the response variable (egg-to-smolt survival) and all the freshwater covariates into a single table, with one row per cohort (spawning year). The tables are joined by cohort (full join, to keep all years present in at least one table).
- **Inputs**:
  - `R/00_raw_data/salmon/Salmon_Monitoring_StJean_Trinite_2024.csv` (column `surv_egg_smolt_pct`)
  - `R/03_clean_data/River/River_Temperatures.fst`
  - `R/03_clean_data/River/River_Discharge.fst`
  - `R/03_clean_data/River/River_Winter_Temperatures.fst`
  - `R/03_clean_data/River/ForestActivity_Cohort.fst`
- **Outputs**:
  - St. Jean: `R/04_analysis/River/Final_Data/PreVariablesSelection_River_StJean.fst` — Published
  - Trinité: `R/04_analysis/River/Final_Data/PreVariablesSelection_River_Trinite.fst` — Published
- **Columns**:
  - `river`, `cohort` (spawning year)
  - `surv_egg_smolt`: egg-to-smolt survival (%), response variable
  - Covariates of water temperature, discharge, freezing degree-days and forestry activity (see the tables of sections 2.2.1 to 2.2.4)

  The covariate names follow these conventions:

  | Prefix or suffix | Meaning |
  |---|---|
  | `co` | Cohort-weighted variable |
  | `gr` | Growth period |
  | `em` | Emergence period |
  | `dd_u` / `dd_o` | Sum of temperatures below / above a threshold |
  | `temp` | Water temperature |
  | `disc` | Discharge |
  | `fdd_per` | Cumulative freezing degree-days from 1 October to the end of the indicated month (air temperature) |
  | `fdd_mon` | Monthly freezing degree-days (air temperature) |
  | `fdd_bimon` | Two-month freezing degree-days (air temperature) |
  | `oc`, `no`, `de`, `ja`, `fe`, `ma`, `ap` | October to April |


#### `02_Variables_Selection_Exploration.R` (river)

- **Description**: selects the covariates of egg-to-smolt survival, then checks the assumptions of beta regression for each river.
  1. A priori removal of redundant covariates:
     - discharge variance;
     - 1 − PGTI indices, perfectly correlated with the PGTI by construction;
     - all FDD metrics except `co_fdd_per_oc_ja` and `co_fdd_bimon_ja_fe`.
  2. VIF-based selection, first within each group (temperature, discharge, FDD), then across all remaining covariates. At each step, the variable with the highest VIF is removed. Correlations are examined with a PCA and a correlation matrix.
  3. Check of the assumptions of the beta regression model (logit link):
     - temporal autocorrelation of the residuals (ACF, Durbin-Watson test);
     - distribution of the response variable;
     - simulated residuals (DHARMa);
     - overdispersion.
  
- **Inputs**:
  - `R/04_analysis/River/Final_Data/PreVariablesSelection_River_StJean.fst`
  - `R/04_analysis/River/Final_Data/PreVariablesSelection_River_Trinite.fst`
- **Outputs**:

  | File | Description | Status |
  |---|---|---|
  | `R/04_analysis/River/Final_Data/Final_DataSet_STJ.fst` | Survival (proportion) and retained covariates, St. Jean | Published |
  | `R/04_analysis/River/Final_Data/Final_DataSet_TRI.fst` | Survival (proportion) and retained covariates, Trinité | Published |

- **Retained covariates**:

  | Covariate | St. Jean | Trinité |
  |---|:---:|:---:|
  | `cohort` | ✓ | ✓ |
  | `co_gr_median_temp` | ✓ | ✓ |
  | `co_gr_var_temp` | | ✓ |
  | `co_gr_dd_u_16_temp` | | ✓ |
  | `co_gr_dd_o_20_temp` | ✓ | |
  | `co_year_mean_sum_ptgi` | ✓ | ✓ |
  | `em_median_temp` | ✓ | ✓ |
  | `em_var_temp` | ✓ | ✓ |
  | `co_gr_mean_disc` | ✓ | ✓ |
  | `co_30d_mean_disch` | | ✓ |
  | `em_quant_75qt_disc` | ✓ | |
  | `co_fdd_bimon_ja_fe` | | ✓ |
  | `co_mean_forest` | ✓ | ✓ |

  Emergence discharge data are not available for the Trinité.

#### `03_River_Modeling.R`

- **Description**: models egg-to-smolt survival by beta regression (logit link) for each river, using the covariates retained by `02_Variables_Selection_Exploration.R`.
  1. Centring and scaling of the covariates. `cohort` stays in years; the scaling parameters are kept to express the figures in original units.
  2. Check of the linearity of the relationships in a beta GAM with smooths (`mgcv`): an edf ≈ 1 indicates a linear relationship.
  3. Fit of the global model, then model selection by AICc (`MuMIn::dredge`) and averaging of the top-set models (ΔAICc < 2), as full average and conditional average.
  4. For each model of the top set:
     - pseudo-R²;
     - leave-one-out cross-validation (LOO-CV), which compares the out-of-sample error and R² with their in-sample values to detect overfitting.
  5. Model-averaged marginal effects of each variable present in the top set:
     - the focal variable varies over its observed range, the other covariates are held at their mean;
     - the predictions of all top-set models are weighted by their AICc weights;
     - the confidence intervals are obtained by percentile bootstrap, because `predict.betareg()` does not provide standard errors.
- **Inputs**:
  - `R/04_analysis/River/Final_Data/Final_DataSet_STJ.fst`
  - `R/04_analysis/River/Final_Data/Final_DataSet_TRI.fst`

#### `04_Creating_Models_Final_Figures_River.R`

- **Description**: produces the article figure of the model-averaged marginal effects for egg-to-smolt survival. The script refits the final global model of the St. Jean (same covariates as `03_River_Modeling.R`), redoes the model selection by AICc and the averaging of the top-set models (ΔAICc < 2), then computes the model-averaged marginal effects with 95% bootstrap confidence intervals. The Trinité has no figure: the null model is part of its top set and no covariate is significant in the model average.
- **Inputs**:
  - `R/04_analysis/River/Final_Data/Final_DataSet_STJ.fst`
- **Outputs**:
  - `R/04_analysis/River/Combined_MarginalEffects_River.png` — Published

### 3.2 Return from sea

Folder: `R/04_analysis/Sea/`

#### `01_Combine_Covariates_Sea.R`

- **Description**: assembles the four datasets of the return-from-sea models (St. Jean and Trinité, 1SW and 2SW), with one row per cohort (year of sea entry). For each river × sea age combination, the response variable (return rate) is joined by cohort to the marine covariates: temperature difference between sea and river at sea entry, food availability (Labrador Sea plankton, year t), thermal habitat indices, smolt condition factor, Pre-Fishery abundance and exploitation rates. For 1SW, only the exploitation rate `er1NAC_lag1` is kept. For 2SW, `er_WG_lag1` and `er2NAC_lag2` are combined into a single exploitation index (`er_index`), because they measure the same phenomenon in two linked areas and are collinear.
- **Inputs**:
  - `R/00_raw_data/salmon/Salmon_Monitoring_StJean_Trinite_2024.csv` (columns `return_1sw_pct`, `return_2sw_pct` and `smolt_k`)
  - `R/03_clean_data/Sea/DeltaT_SST_River_{stj,tri}_{1SW,2SW}.fst`
  - `R/03_clean_data/Sea/Plankton_LabradorSea_{stj,tri}_{1SW,2SW}.fst`
  - `R/03_clean_data/Sea/TH_annual_{stj,tri}_{1sw,2sw}.fst`
  - `R/03_clean_data/Sea/ExploitationRate_Sea.fst`
  - `R/03_clean_data/Sea/Density_Lab.fst`
- **Outputs**:
  - `R/04_analysis/Sea/Final_Data/PreSelection_DataSet_1sw_stj.fst` — Published
  - `R/04_analysis/Sea/Final_Data/PreSelection_DataSet_2sw_stj.fst` — Published
  - `R/04_analysis/Sea/Final_Data/PreSelection_DataSet_1sw_tri.fst` — Published
  - `R/04_analysis/Sea/Final_Data/PreSelection_DataSet_2sw_tri.fst` — Published
- **Columns**:

  | Column | Files | Description |
  |---|---|---|
  | `river`, `cohort` | All | River and year of sea entry |
  | `return_1sw` / `return_2sw` | 1SW / 2SW | Smolt-to-adult return rate (%), response variable |
  | `delta_T_sst_minus_river` | All | Temperature difference between sea and river at sea entry |
  | `Plankton_t` | All | Zooplankton biomass index in the Labrador Sea, year t |
  | `StJeanCorridor_TH` / `TriniteCorridor_TH`, `BelleIsle_TH`, `LabradorEntrance_TH`, `LabradorSea_TH` | All | Thermal habitat indices per polygon |
  | `smolt_k` | All | Mean Fulton's condition factor of the smolts |
  | `pfa_1sw` / `pfa_2sw` | 1SW / 2SW | Pre-Fishery abundance (thousands of fish); removed at the variable selection step (VIF too high) |
  | `er1NAC_lag1` | 1SW | Exploitation rate of North American 1SW, year t + 1 |
  | `er_index` | 2SW | 2SW exploitation index (see below) |

- **Variables created**:

  | Variable | Period | Description |
  |---|---|---|
  | `er_index` | Years t + 1 and t + 2 | Score on the first axis of a PCA (centred and scaled variables) of the log of `er_WG_lag1` and the log of `er2NAC_lag2`, computed separately for each river. The sign is oriented so that a high value corresponds to high fishing pressure. The script checks the correlation between the two rates (criterion: r ≥ 0.70), the variance explained by the first axis (target: ≥ 70%) and the sign of the contributions of the two variables. |

#### `02_Variables_Selection_Exploration.R` (sea)

- **Description**: selects the covariates of the 1SW and 2SW return rates, then checks the assumptions of beta regression for each river.
  1. Exploration of the correlations with a PCA and a correlation matrix, first for the surface temperatures, then for all covariates. No variable is removed on the basis of correlations alone.
  2. VIF-based selection.
  3. Check of the assumptions of the beta regression model (logit link):
     - distribution of the response variable;
     - simulated residuals (DHARMa);
     - temporal autocorrelation (ACF, Durbin-Watson test);
     - overdispersion.

- **Inputs**:
  - `R/04_analysis/Sea/Final_Data/PreSelection_DataSet_1sw_stj.fst`
  - `R/04_analysis/Sea/Final_Data/PreSelection_DataSet_2sw_stj.fst`
  - `R/04_analysis/Sea/Final_Data/PreSelection_DataSet_1sw_tri.fst`
  - `R/04_analysis/Sea/Final_Data/PreSelection_DataSet_2sw_tri.fst`
- **Outputs**:

  | File | Description | Status |
  |---|---|---|
  | `R/04_analysis/Sea/Final_Data/Final_DataSet_STJ_1SW.fst` | 1SW return rate (proportion) and covariates, St. Jean | Published |
  | `R/04_analysis/Sea/Final_Data/Final_DataSet_STJ_2SW.fst` | 2SW return rate (proportion) and covariates, St. Jean | Published |
  | `R/04_analysis/Sea/Final_Data/Final_DataSet_TRI_1SW.fst` | 1SW return rate (proportion) and covariates, Trinité | Published |
  | `R/04_analysis/Sea/Final_Data/Final_DataSet_TRI_2SW.fst` | 2SW return rate (proportion) and covariates, Trinité | Published |


- **Retained covariates**:

  | Covariate | St. Jean 1SW | Trinité 1SW | St. Jean 2SW | Trinité 2SW |
  |---|:---:|:---:|:---:|:---:|
  | `cohort` | ✓ | ✓ | ✓ | ✓ |
  | `delta_T_sst_minus_river` | ✓ | ✓ | ✓ | ✓ |
  | `StJeanCorridor_TH` / `TriniteCorridor_TH` | ✓ | ✓ | ✓ | ✓ |
  | `BelleIsle_TH` | ✓ | ✓ | ✓ | ✓ |
  | `LabradorEntrance_TH` | ✓ | ✓ | ✓ | ✓ |
  | `LabradorSea_TH` | ✓ | ✓ | ✓ | ✓ |
  | `Plankton_t` | ✓ | ✓ | ✓ | ✓ |
  | `smolt_k` | ✓ | ✓ | ✓ | ✓ |
  | `er1NAC_lag1` | ✓ | ✓ | | |
  | `er_index` | | | ✓ | |

#### `03_Sea_Modeling.R`

- **Description**: models the smolt-to-adult return rates, after one winter at sea (1SW) and after two winters (2SW), by beta regression (logit link) for each river. The models use the covariates retained by `02_Variables_Selection_Exploration.R` (sea). For each of the four river × sea age combinations:
  1. Centring and scaling of the covariates. `cohort` stays in years; the scaling parameters are kept to express the figures in original units.
  2. Check of the linearity of the relationships in a beta GAM with smooths (`mgcv`).
  3. Fit of the global model, then model selection by AICc (`MuMIn::dredge`) and averaging of the top-set models (ΔAICc < 2). The relative importance of the variables (sum of AICc weights) and the confidence intervals of the coefficients are also computed.
  4. For each model of the top set:
     - pseudo-R²;
     - leave-one-out cross-validation (LOO-CV), which compares the out-of-sample error and R² with their in-sample values to detect overfitting.
  5. Model-averaged marginal effects of each variable present in the top set:
     - the focal variable varies over its observed range, the other covariates are held at their mean;
     - the predictions of all top-set models are weighted by their AICc weights;
     - the confidence intervals are obtained by percentile bootstrap, because `predict.betareg()` does not provide standard errors.
- **Inputs**:
  - `R/04_analysis/Sea/Final_Data/Final_DataSet_STJ_1SW.fst`
  - `R/04_analysis/Sea/Final_Data/Final_DataSet_STJ_2SW.fst`
  - `R/04_analysis/Sea/Final_Data/Final_DataSet_TRI_1SW.fst`
  - `R/04_analysis/Sea/Final_Data/Final_DataSet_TRI_2SW.fst`

#### `04_Creating_Models_Final_Figures_Sea.R`

- **Description**: produces the article figure of the model-averaged marginal effects for the 1SW and 2SW return rates. The script refits the four final global models (same formulas as `03_Sea_Modeling.R`), redoes the model selection by AICc and the averaging of the top-set models (ΔAICc < 2), then computes the model-averaged marginal effects with 95% bootstrap confidence intervals.
- **Inputs**:
  - `R/04_analysis/Sea/Final_Data/Final_DataSet_STJ_1SW.fst`
  - `R/04_analysis/Sea/Final_Data/Final_DataSet_STJ_2SW.fst`
  - `R/04_analysis/Sea/Final_Data/Final_DataSet_TRI_1SW.fst`
  - `R/04_analysis/Sea/Final_Data/Final_DataSet_TRI_2SW.fst`
- **Outputs**:
  - `R/04_analysis/Sea/Combined_MarginalEffects_Sea.png` — Published
  

---

## 4. Software environment

- **R**: 4.5.2
- **Packages**:

  | Package | Version |
  |---|---|
  | `betareg` | 3.2.4 |
  | `car` | 3.1.5 |
  | `cowplot` | 1.2.0 |
  | `data.table` | 1.18.2.1 |
  | `DHARMa` | 0.4.7 |
  | `dplyr` | 1.2.0 |
  | `ecmwfr` | 2.0.3 |
  | `elevatr` | 0.99.1 |
  | `factoextra` | 2.0.0 |
  | `FactoMineR` | 2.13 |
  | `fst` | 0.9.8 |
  | `future.apply` | 1.20.2 |
  | `gdistance` | 1.6.5 |
  | `ggplot2` | 4.0.2 |
  | `ggpubr` | 0.6.3 |
  | `ggrepel` | 0.9.6 |
  | `ggridges` | 0.5.7 |
  | `glmmTMB` | 1.1.14 |
  | `glmnet` | 4.1.10 |
  | `here` | 1.0.2 |
  | `lmtest` | 0.9.40 |
  | `lubridate` | 1.9.5 |
  | `lwgeom` | 0.2.15 |
  | `maps` | 3.4.3 |
  | `mgcv` | 1.9.3 |
  | `modi` | 0.1.3 |
  | `MuMIn` | 1.48.11 |
  | `ncdf4` | 1.24 |
  | `nngeo` | 0.4.8 |
  | `osmdata` | 0.3.0 |
  | `patchwork` | 1.3.2 |
  | `purrr` | 1.2.1 |
  | `raster` | 3.6.32 |
  | `readxl` | 1.4.5 |
  | `reshape2` | 1.4.5 |
  | `rnaturalearth` | 1.2.0 |
  | `rnaturalearthdata` | 1.0.0 |
  | `rstudioapi` | 0.18.0 |
  | `scales` | 1.4.0 |
  | `sf` | 1.1.0 |
  | `slider` | 0.3.3 |
  | `smoothr` | 1.3.0 |
  | `stringr` | 1.6.0 |
  | `terra` | 1.8.93 |
  | `tibble` | 3.3.1 |
  | `tidyr` | 1.3.2 |
  | `tidyverse` | 2.0.0 |
  | `units` | 1.0.0 |
  | `viridis` | 0.6.5 |
  | `whitebox` | 2.4.3 |
  | `writexl` | 1.5.4 |
  | `zoo` | 1.8.15 |