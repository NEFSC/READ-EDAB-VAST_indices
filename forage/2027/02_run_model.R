# This is the VAST code used in Gaichas et al 2023, but with data years
# extended from 1973-2024
# AND NO CUSTOM EXTRAPOLATION GRID SO NO 3-nmile INSHORE STRATA
# VAST attempt 2 univariate model as a script
# modified from https://github.com/James-Thorson-NOAA/VAST/wiki/Index-standardization

## 2026-10-02 AST

# Load packages
library(here)
library(dplyr)
library(VAST)

source(here::here("R/utils.R"))

## Data wrangling ----
#Read in data, separate spring and fall, and rename columns for VAST

bluepyagg_stn <- readRDS(here::here(
  "forage/2027/outputs/test_VAST_input.rds"
))

# make SST column that uses surftemp unless missing or 0
# there are 3 surftemp 0 values in the dataset, all with oisst > 15
bluepyagg_stn <- bluepyagg_stn |>
  dplyr::mutate(
    sstfill = ifelse((is.na(surftemp) | surftemp == 0), oisst, surftemp)
  )

# filter to assessment years at Tony's suggestion
# code Vessel as AL=0, HB=1, NEAMAP=2

bluepyagg_stn_fall <- bluepyagg_stn |>
  filter(season_ng == "FALL") |>
  mutate(
    AreaSwept_km2 = 1, #Elizabeth's code
    Vessel = as.numeric(as.factor(vessel)) - 1
  ) |>
  dplyr::select(
    Catch_g = meanbluepywt, #use bluepywt for individuals input in example
    Year = year,
    Vessel,
    AreaSwept_km2,
    Lat = declat,
    Lon = declon,
    meanpisclen,
    npiscsp,
    oisst, #leaves out everything before 1982
    sstfill
  ) |>
  na.omit() |>
  as.data.frame()

bluepyagg_stn_spring <- bluepyagg_stn |>
  filter(season_ng == "SPRING") |>
  mutate(
    AreaSwept_km2 = 1, #Elizabeth's code
    Vessel = as.numeric(as.factor(vessel)) - 1
  ) |>
  dplyr::select(
    Catch_g = meanbluepywt,
    Year = year,
    Vessel,
    AreaSwept_km2,
    Lat = declat,
    Lon = declon,
    meanpisclen,
    npiscsp,
    oisst, #leaves out everything before 1982
    sstfill
  ) |>
  na.omit() |>
  as.data.frame()


## Make settings  ----
# (turning off bias.correct to save time for example)
# NEFSC strata limits https://github.com/James-Thorson-NOAA/VAST/issues/302

MAB <- c(1010:1080, 1100:1120, 1600:1750, 3010:3450, 3470, 3500, 3510)
GB <- c(1090, 1130:1210, 1230, 1250, 3460, 3480, 3490, 3520:3550)
GOM <- c(1220, 1240, 1260:1290, 1360:1400, 3560:3830)
SS <- c(1300:1352, 3840:3990)

# using VAST built in grid because custom one is breaking Oct 2024
MAB2 <- FishStatsUtils::northwest_atlantic_grid |>
  dplyr::filter(stratum_number %in% MAB) |>
  dplyr::select(stratum_number) |>
  dplyr::distinct()

# Georges Bank EPU
GB2 <- FishStatsUtils::northwest_atlantic_grid |>
  dplyr::filter(stratum_number %in% GB) |>
  dplyr::select(stratum_number) |>
  dplyr::distinct()

# gulf of maine EPU (for SOE)
GOM2 <- FishStatsUtils::northwest_atlantic_grid |>
  dplyr::filter(stratum_number %in% GOM) |>
  dplyr::select(stratum_number) |>
  dplyr::distinct()

# scotian shelf EPU (for SOE)
SS2 <- FishStatsUtils::northwest_atlantic_grid |>
  dplyr::filter(EPU == "Scotian_Shelf") |>
  dplyr::select(stratum_number) |>
  dplyr::distinct()

# needed to cover the whole northwest atlantic grid--lets try without
allother2 <- FishStatsUtils::northwest_atlantic_grid |>
  dplyr::filter(
    !EPU %in%
      c('Mid_Atlantic_Bight', 'Georges_Bank', 'Gulf_of_Maine', 'Scotian_Shelf')
  ) |>
  dplyr::select(stratum_number) |>
  dplyr::distinct()

# all epus
allEPU2 <- FishStatsUtils::northwest_atlantic_grid |>
  dplyr::filter(
    EPU %in%
      c('Mid_Atlantic_Bight', 'Georges_Bank', 'Gulf_of_Maine', 'Scotian_Shelf')
  ) |>
  dplyr::select(stratum_number) |>
  dplyr::distinct()

# default configs, not really specified anyway
FieldConfig = matrix(
  "IID",
  ncol = 2,
  nrow = 3,
  dimnames = list(
    c("Omega", "Epsilon", "Beta"),
    c("Component_1", "Component_2")
  )
)


RhoConfig <- c(
  "Beta1" = 0, # temporal structure on years (intercepts)
  "Beta2" = 0,
  "Epsilon1" = 0, # temporal structure on spatio-temporal variation
  "Epsilon2" = 0
)
# 0 off (fixed effects)
# 1 independent
# 2 random walk
# 3 constant among years (fixed effect)
# 4 AR1

OverdispersionConfig <- c("eta1" = 0, "eta2" = 0)
# eta1 = vessel effects on prey encounter rate
# eta2 = vessel effects on prey weight

strata.limits <- as.list(c(
  "AllEPU" = allEPU2,
  "MAB" = MAB2,
  "GB" = GB2,
  "GOM" = GOM2,
  "SS" = SS2,
  "allother" = allother2
))

settings = make_settings(
  n_x = 500,
  Region = "northwest_atlantic",
  strata.limits = strata.limits,
  purpose = "index2",
  bias.correct = TRUE,
  OverdispersionConfig = OverdispersionConfig,
  mesh_package = 'fmesher'
)

## Run model fall ----

working_dir <- here::here(sprintf("forage/2027/outputs/fall_model"))

if (!dir.exists(working_dir)) {
  dir.create(working_dir)
}

fit <- fit_model(
  settings = settings,
  Lat_i = bluepyagg_stn_fall$Lat,
  Lon_i = bluepyagg_stn_fall$Lon,
  t_i = bluepyagg_stn_fall$Year,
  b_i = as_units(bluepyagg_stn_fall[, 'Catch_g'], 'g'),
  a_i = rep(1, nrow(bluepyagg_stn_fall)),
  v_i = bluepyagg_stn_fall$Vessel,
  Q_ik = as.matrix(bluepyagg_stn_fall[, c(
    "npiscsp",
    "meanpisclen",
    "sstfill"
  )]),
  working_dir = paste0(working_dir, "/")
)

saveRDS(fit, file = paste0(working_dir, "/fit.rds"))
fit <- readRDS(paste0(working_dir, "/fit.rds"))

# Plot results
plot(fit, working_dir = paste0(working_dir, "/"))

# extract center of gravity
cog <- extract_cog(fit)
saveRDS(cog, here::here(working_dir, paste0("fall_cog.rds")))

## Run model spring ----

working_dir <- here::here(sprintf(
  "forage/2027/outputs/spring_model"
))

if (!dir.exists(working_dir)) {
  dir.create(working_dir)
}

fit <- fit_model(
  settings = settings,
  Lat_i = bluepyagg_stn_spring[, 'Lat'],
  Lon_i = bluepyagg_stn_spring[, 'Lon'],
  t_i = bluepyagg_stn_spring[, 'Year'],
  b_i = as_units(bluepyagg_stn_spring[, 'Catch_g'], 'g'),
  a_i = rep(1, nrow(bluepyagg_stn_spring)),
  v_i = bluepyagg_stn_spring$Vessel,
  Q_ik = as.matrix(bluepyagg_stn_spring[, c(
    "npiscsp",
    "meanpisclen",
    "sstfill"
  )]),
  working_dir = paste0(working_dir, "/")
)

saveRDS(fit, file = paste0(working_dir, "/fit.rds"))
fit <- readRDS(paste0(working_dir, "/fit.rds"))

# Plot results
plot(fit, working_dir = paste0(working_dir, "/"))

# extract center of gravity
cog <- extract_cog(fit)
saveRDS(cog, here::here(working_dir, paste0("spring_cog.rds")))

## Create SOE indices ----

### abundance indices ----
stratlook_EPUonly <- data.frame(
  Stratum = c(
    "Stratum_1",
    "Stratum_2",
    "Stratum_3",
    "Stratum_4",
    "Stratum_5",
    "Stratum_6"
  ),
  Region = c("AllEPU", "MAB", "GB", "GOM", "SS", "allother")
)

SOEinputs(
  infile = "forage/2027/outputs/fall_model/Index.csv",
  season = "Fall",
  stratlook = stratlook_EPUonly,
  outfile = paste0(
    "forage/2027/outputs/fallforageindex_",
    Sys.Date(),
    ".rds"
  )
)


SOEinputs(
  infile = "forage/2027/outputs/spring_model/Index.csv",
  season = "Spring",
  stratlook = stratlook_EPUonly,
  outfile = paste0(
    "forage/2027/outputs/springforageindex_",
    Sys.Date(),
    ".rds"
  )
)

### center of gravity indices ----

SOEinputsCOG(
  infile = "forage/2027/outputs/fall_model/fall_cog.rds",
  season = "Fall",
  outfile = paste0("forage/2027/outputs/fallforagecog_", Sys.Date(), ".rds")
)


SOEinputsCOG(
  infile = "forage/2027/outputs/spring_model/spring_cog.rds",
  season = "Spring",
  outfile = paste0(
    "forage/2027/outputs/springforagecog_",
    Sys.Date(),
    ".rds"
  )
)
