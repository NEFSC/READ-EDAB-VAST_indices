# --- Script Header ---
# Title: Create VAST Inputs for Forage Index Indicators
# Author: AST
# Date: October 2026
# Description: This script processes and integrates Northeast Fisheries Science Center (NEFSC) and
# Northeast Marine Area Monitoring and Assessment Program (NEAMAP) fish stomach content data with OISST
# sea surface temperature (SST) data. The primary output is a combined dataset
# aggregated at the station level, suitable for VAST modeling or other analyses of
# forage index indicators for the State of the Ecosystem report.

# Libraries & functions ----

library(tidyverse)
library(here)
library(dendextend)
library(sf) # For spatial data manipulation
library(raster) # For handling raster data
library(terra) # Another option for raster data (modern alternative)
library(nngeo) # For nearest-neighbor spatial joins

# load custom utility functions
source(here::here("R/utils.R"))

# Data ----

# load data
load(here::here("shared_data/2027/allfh.rmd.epu.Rdata"))

# csvs

NEAMAPstationSST <- read.csv(here::here(
  "shared_data/2027/NEAMAP SST_2007_2025.csv"
))

## read in piscivore predator list
pisccompletedf <- read.csv(here::here("forage/static/pisccomplete.csv"))

# NEAMAP data
NEAMAPblueprey_csv <- read.csv(here::here(
  "forage/static/Full Prey List_Common Names.csv"
))

# Analyses ----

## Update prey list ----

# This section identifies the key predators (piscivores) and prey species
# based on diet overlap and observation counts.

# Get prey lists from both NEFSC and NEAMAP data.
# this analysis uses a static prey list from NEAMAP
# but theoretically could change for NEFSC each year because it uses the new allfh data

### simplify identification of NEFSC prey items ----
# prey = the things bluefish eat
NEFSCprey <- allfh |>
  dplyr::filter(pdcomnam == "BLUEFISH") |>
  group_by(pynam, pycomnam, gencom) |>
  summarise(count = n()) |>
  filter(
    !pynam %in%
      c(
        "EMPTY",
        "BLOWN",
        "FISH",
        "OSTEICHTHYES",
        "ANIMAL REMAINS",
        "FISH SCALES"
      ),
    !gencom %in%
      c(
        "ARTHROPODA",
        "ANNELIDA",
        "CNIDARIA",
        "UROCHORDATA",
        "ECHINODERMATA",
        "WORMS",
        "BRACHIOPODA",
        "COMB JELLIES",
        "BRYOZOA",
        "SPONGES",
        "MISCELLANEOUS",
        "OTHER"
      )
  ) |>
  arrange(desc(count)) |>
  dplyr::mutate(pynam2 = tolower(pynam)) |>
  dplyr::rename(NEFSC = count)


### get NEAMAP prey ----

# March 2023, formally add NEAMAP to prey decisions
NEAMAPprey <- NEAMAPblueprey_csv |>
  filter(
    !SCIENTIFIC.NAME %in%
      c(
        "Actinopterygii",
        "fish scales",
        "Decapoda (megalope)",
        "unidentified material",
        "Plantae",
        "unidentified animal"
      ),
    !is.na(BLUEFISH)
  ) |>
  dplyr::select(COMMON.NAME, SCIENTIFIC.NAME, BLUEFISH) |>
  dplyr::mutate(
    pynam2 = tolower(SCIENTIFIC.NAME),
    pynam2 = stringr::str_replace(pynam2, "spp.", "sp")
  ) |>
  dplyr::rename(NEAMAP = BLUEFISH)

### combine prey lists ----

# new criteria March 2023, >20 observations NEAMAP+NEFSC, but keep mackerel
# removes the flatfish order (too broad) and unid Urophycis previously in NEAMAP
blueprey <- NEFSCprey |>
  dplyr::full_join(NEAMAPprey) |>
  dplyr::mutate(
    # NEAMAP = ifelse(is.na(NEAMAP), 0, NEAMAP),
    # NEFSC = ifelse(is.na(NEFSC), 0, NEFSC),
    total = sum(NEFSC, NEAMAP, na.rm = TRUE),
    PREY = coalesce(SCIENTIFIC.NAME, pynam),
    COMMON = coalesce(pycomnam, COMMON.NAME),
    pynam = coalesce(pynam, toupper(pynam2))
  ) |>
  dplyr::arrange(desc(total)) |>
  dplyr::filter(total > 20 | pynam == "SCOMBER SCOMBRUS") |> # >20 leaves out mackerel
  dplyr::mutate(
    COMMON = case_when(
      pynam == "ILLEX SP" ~ "Shortfin squids",
      pynam2 == "teuthida" ~ "Unidentified squids",
      TRUE ~ COMMON
    )
  ) |>
  dplyr::mutate(
    PREY = stringr::str_to_sentence(PREY),
    COMMON = stringr::str_to_sentence(COMMON)
  )

### Classify predators and prey in NEFSC data ----
# Filter the raw food habits data to include only the defined piscivore predators.

survey_predator_prey <- allfh |>
  # join feeding guild from pisccompletedf
  left_join(
    pisccompletedf,
    by = c("pdcomnam" = "COMNAME", "sizecat" = "SizeCat")
  ) |>
  # filter out rows that aren't piscivores
  filter(!is.na(feedguild)) |>
  # identify which ones are bluefish prey
  mutate(
    blueprey = case_when(
      pynam %in% blueprey$pynam ~ "blueprey",
      TRUE ~ "othprey"
    )
  )

## Aggregate NEFSC Data at the Station Level ----
# This section aggregates the detailed stomach content data to create a single row
# per survey station, calculating various summary statistics.

bluepyall_stn <- survey_predator_prey |>
  #create id linking cruise6_station
  #create season_ng spring and fall Spring=Jan-May, Fall=June-Dec
  mutate(
    id = paste0(cruise6, "_", station),
    year = as.numeric(year),
    month = as.numeric(month),
    day = as.numeric(day),
    season_ng = case_when(
      month <= 6 ~ "SPRING",
      month >= 7 ~ "FALL",
      TRUE ~ as.character(NA)
    )
  ) |>
  dplyr::select(
    year,
    month,
    day,
    season_ng,
    id,
    stratum,
    pynam,
    pyamtw,
    pywgti,
    pyvoli,
    blueprey,
    pdcomnam,
    pdid,
    pdlen,
    pdsvol,
    pdswgt,
    beglat,
    beglon,
    declat,
    declon,
    bottemp,
    surftemp,
    setdepth
  ) |>
  group_by(id) |>
  #mean blueprey g per stomach per tow: sum all blueprey g/n stomachs in tow
  mutate(
    bluepywt = case_when(blueprey == "blueprey" ~ pyamtw, TRUE ~ 0.0),
    bluepynam = case_when(blueprey == "blueprey" ~ pynam, TRUE ~ NA_character_)
  )

# Now get station data in one line
stndat <- bluepyall_stn |>
  dplyr::select(
    year,
    season_ng,
    id,
    beglat,
    beglon,
    declat,
    declon,
    bottemp,
    surftemp,
    setdepth
  ) |>
  distinct()

#pisc stomachs in tow count pdid for each pred and sum
piscstom <- bluepyall_stn |>
  group_by(id, pdcomnam) |>
  summarise(nstompd = n_distinct(pdid)) |>
  group_by(id) |>
  summarise(nstomtot = sum(nstompd))

bluepyall_stn |>
  dplyr::filter(id == "197303_154") |>
  dplyr::select(id, pdcomnam, pdid) |>
  group_by(id, pdcomnam) |>
  summarise(nstompd = n_distinct(pdid)) |>
  group_by(id) |>
  summarise(nstomtot = sum(nstompd))

#mean and var pred length per tow
pisclen <- bluepyall_stn |>
  summarise(meanpisclen = mean(pdlen), varpisclen = var(pdlen))

# Aggregated prey at station level with predator covariates
bluepyagg_stn <- bluepyall_stn |>
  summarise(
    sumbluepywt = sum(bluepywt),
    nbluepysp = n_distinct(bluepynam, na.rm = T),
    npreysp = n_distinct(pynam),
    npiscsp = n_distinct(pdcomnam)
  ) |>
  left_join(piscstom) |>
  mutate(meanbluepywt = sumbluepywt / nstomtot) |>
  left_join(pisclen) |>
  left_join(stndat)

bluepyagg_stn2 <- bluepyall_stn |>
  group_by(id) |>
  summarise(
    year = first(year),
    season_ng = first(season_ng),
    month = first(month),
    day = first(day),
    declat = first(declat),
    declon = first(declon),
    beglat = first(beglat),
    beglon = first(beglon),
    bottemp = first(bottemp),
    surftemp = first(surftemp),
    setdepth = first(setdepth),
    sumbluepywt = sum(bluepywt, na.rm = TRUE),
    nbluepysp = n_distinct(bluepynam, na.rm = TRUE),
    npreysp = n_distinct(pynam, na.rm = TRUE),
    npiscsp = n_distinct(pdcomnam, na.rm = TRUE),
    nstomtot = n_distinct(pdcomnam, pdid, na.rm = TRUE),
    meanpisclen = mean(pdlen, na.rm = TRUE),
    varpisclen = var(pdlen, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    meanbluepywt = sumbluepywt / nstomtot,
    vessel = ifelse(year < 2009, "AL", "HB")
  )

colnames <- colnames(bluepyagg_stn)
head(bluepyagg_stn) == head(bluepyagg_stn2 |> dplyr::select(colnames))

## all values are either TRUE or NA (missing lengths)

(bluepyagg_stn == (bluepyagg_stn2 |> dplyr::select(colnames))) |>
  as.numeric() |>
  sum(na.rm = TRUE)
(bluepyagg_stn == (bluepyagg_stn2 |> dplyr::select(colnames))) |>
  as.numeric() |>
  length()
(bluepyagg_stn == (bluepyagg_stn2 |> dplyr::select(colnames))) |>
  as.numeric() |>
  is.na() |>
  sum()

# current dataset, fix declon, add vessel, rename NEFSC
#nefsc_bluepyagg_stn <- readRDS(here("fhdat/bluepyagg_stn.rds")) |>
nefsc_bluepyagg_stn <- bluepyagg_stn |>
  mutate(
    declon = -declon,
    vessel = case_when(
      year < 2009 ~ "AL",
      year >= 2009 ~ "HB",
      TRUE ~ as.character(NA)
    )
  )

## Combine NEFSC and NEAMAP Datasets ----
# This section reads in the NEAMAP data and combines it with the processed NEFSC data.

# Read in NEAMAP updated input from Jim Gartland, reformat with same names
neamap_bluepreyagg_stn <- process_neamap_data(
  "shared_data/2027/NEAMAP_Mean stomach weights_Bluefish Prey_Aug2026.csv"
)

# combine NEAMAP and NEFSC
bluepyagg_stn_all <- nefsc_bluepyagg_stn |>
  bind_rows(neamap_bluepreyagg_stn)


## Add month and day to NEFSC data ----
# year and month were dropped earlier in the analysis, not sure why
# but they are added back in here

# get NEFSC station and date data
NEFSCstations <- allfh |>
  dplyr::mutate(
    id = paste0(cruise6, "_", station),
    year = as.numeric(year),
    month = as.numeric(month),
    day = as.numeric(day)
  ) |>
  dplyr::select(id, year, month, day) |>
  dplyr::distinct()

# NEAMAP station data comes from SST data
NEAMAPstations <- NEAMAPstationSST |>
  dplyr::mutate(
    id = station,
    year = as.numeric(year),
    month = as.numeric(month),
    day = as.numeric(day)
  ) |>
  dplyr::select(id, year, month, day) |>
  dplyr::distinct()

Allstations <- bind_rows(NEFSCstations, NEAMAPstations)

#station id, lat lon, year month day

# remake diethauls
# id, lat, long from bluepyagg_stn_all
diethauls <- bluepyagg_stn_all |>
  dplyr::select(id, declat, declon) |>
  # add year month day from Allstations
  dplyr::left_join(Allstations) |>
  dplyr::distinct()

# add year month day to diet data
bluepyagg_stn_all <- left_join(bluepyagg_stn_all, diethauls)

## Fill missing NEFSC SSTs with NEAMAP temperature data ----
# Add SST into NEAMAP and reintegrate into full dataset

# add NEAMAP SST to surftemp field
NEAMAPidSST <- NEAMAPstationSST |>
  mutate(id = station) |>
  dplyr::select(id, SST)

bluepyagg_stn_all <- left_join(bluepyagg_stn_all, NEAMAPidSST, by = "id") |>
  mutate(surftemp = coalesce(surftemp, SST)) |>
  dplyr::select(-SST)

## Integrate OISST Sea Surface Temperature Data ----
# This section adds OISST data to the combined dataset by finding the nearest
# SST measurement in time and space for each station.

# Download and process OISST data for the specified years.
# only run for new year(s)

process_oisst(
  years = 2025,
  out_dir = here::here("forage/static/sst"),
  varname = "sst",
  input = here::here("forage\\2027\\inputs\\sst.day.mean.2025.nc")
)

stations <- bluepyagg_stn_all |>
  dplyr::mutate(
    day = str_pad(day, 2, pad = '0'),
    month = str_pad(month, 2, pad = '0'),
    yrmody = as.numeric(paste0(year, month, day))
  ) |>
  dplyr::select(id, declon, declat, year, yrmody) |>
  na.omit() |>
  sf::st_as_sf(coords = c("declon", "declat"), crs = 4326, remove = FALSE)


#list of SST dataframes
SSTdfs <- list.files(
  here("forage/static/sst"),
  pattern = "*.rds",
  full.names = TRUE
)


dietstn_OISST <- join_oisst_to_stations(
  stations = stations,
  oisst_files = SSTdfs
)

# # join 2025 stations only
# dietstn_OISST <- join_oisst_to_stations(
#   stations = stations |>
#     dplyr::filter(year == 2025),
#   oisst_files = SSTdfs[44]
# )
#
# # last year's joined station temperatures
# historic <- readRDS(here("forage/data/dietstn_OISST_1982_2024.rds"))

# takes >10 minutes to run, save and tack on next years rather than full merge?
# saveRDS(dietstn_OISST, here("forage/2027/outputs/dietstn_OISST_1982_2025.rds"))

## Merge OISST into diet data ----
final_data <- left_join(
  bluepyagg_stn_all,
  dietstn_OISST |>
    dplyr::select(id, oisst = sst) |>
    sf::st_drop_geometry(),
  by = "id"
)

# Save the final dataset ----
saveRDS(
  final_data,
  here::here("forage/2027/outputs/test_VAST_input.rds")
)
