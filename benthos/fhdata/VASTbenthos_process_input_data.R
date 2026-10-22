# VASTbenthos_process_input_data.R
#
# Builds station-level "benthos in benthivore stomachs" datasets (NEFSC + NEAMAP)
# with bottom temperature, used as input to the macrobenthos and megabenthos
# VAST indices for the State of the Ecosystem report.
#
# Original workflow: S. Gaichas (NOAA-EDAB/benthosindex). Maintained by M. Grezlik.
#
# Outputs (benthos/fhdata/) keep their legacy names so the VAST scripts still work:
#   {macro,mega}benall_stn.rds            prey-level NEFSC data
#   {macro,mega}benagg_stn_all.rds        station-level NEFSC + NEAMAP
#   {macro,mega}benagg_stn_all_modBT.rds  + GLORYS bottom temp  <- VAST input

library(dplyr)
library(dendextend)

## ---- Config: the only section that should need editing each year ----------

soe_year  <- 2027
data_year <- soe_year - 2   # last year of survey data in this SOE

paths <- list(
  allfh        = here::here("shared_data", soe_year, "allfh.rmd.epu.RData"),        # NEFSC stomachs (Brian Smith)
  neamap_dates = here::here("shared_data", soe_year, "NEAMAP SST_2007_2025.csv"),   # NEAMAP station dates
  neamap = c(                                                                       # NEAMAP stomachs (Jim Gartland)
    macro = here::here("benthos/data", soe_year, "NEAMAP_Mean Stomach Weights_Macrobenthos prey_wWQ_Aug2026.csv"),
    mega  = here::here("benthos/data", soe_year, "NEAMAP_Mean Stomach Weights_Megabenthos prey_wWQ_Aug2026.csv")
  ),
  dietoverlap = here::here("benthos/data-raw/tgmat.2022-02-15.csv"),  # predator diet overlap matrix
  sizedef     = here::here("benthos/data-raw/Table3.xlsx"),           # NEFSC FH size categories (Brian)
  prey        = here::here("benthos/data/prey.RData"),                # Rpath prey mappings, June 2023
  glorys_dir  = "//nefscdata/EDAB_Datasets/GLORYS/V12/PRODUCTS/annual_bottomT_vast",
  out_dir     = here::here("benthos/fhdata")
)

glorys_v2_start <- 1993   # Joe's V2 files correct an error in 1993+; V1 used before that

# Methods switches. Defaults reproduce the methods used for the 2026 SOE.
filter_benthivores <- FALSE  # TRUE = keep only benthivore predators in NEFSC data (methods question, see issue)

missing <- unlist(paths[c("allfh", "neamap_dates", "neamap", "dietoverlap", "sizedef", "prey")])
missing <- missing[!file.exists(missing)]
if (length(missing)) stop("Missing input files:\n", paste(missing, collapse = "\n"))

## ---- Predator list: benthivores from diet-overlap clustering --------------

load(paths$allfh)   # object is called `allfh`

dietoverlap <- readr::read_csv(paths$dietoverlap, show_col_types = FALSE)

dend <- as.dendrogram(hclust(dist(dietoverlap)))
dend <- rotate(dend, 1:136)
dend <- color_branches(dend, k = 6)  # Brian uses 6 categories
labels(dend) <- paste0(as.character(names(dietoverlap[-1]))[order.dendrogram(dend)],
                       "(", labels(dend), ")")
dend <- hang.dendrogram(dend, hang_height = 0.1)

# predators under the node containing the given leaves, as COMNAME / SizeCat
guild_members <- function(leaves) {
  members <- partition_leaves(dend)[[which_node(dend, leaves)]]
  data.frame(COMNAME = toupper(stringr::str_remove(members, "\\..*")),
             SizeCat = stringr::str_remove(stringr::str_extract(members, "\\..*[:upper:]+"), "\\.."))
}

piscdf  <- guild_members(c("Bluefish..S(37)", "Bluefish..M(36)", "Bluefish..L(35)"))
plankdf <- guild_members(c("Atlantic herring..S(20)", "Atlantic mackerel..S(23)", "Blueback herring..XS(34)"))
alldf   <- guild_members(c("Barndoor skate..S(26)", "Bluefish..L(35)"))

sizedef <- readxl::read_excel(paths$sizedef, range = "B3:H55") |>
  rename(XS = `Extra-Small`, S = "Small", M = "Medium", L = "Large", XL = `Extra-Large`,
         CommonName = `Common Name`, SpeciesName = `Species Name`) |>
  tidyr::pivot_longer(-c(1:2), names_to = "sizecat", values_to = "range") |>
  filter(!range == "-") |>
  # predators missing from TM216; Brian: default is S (<=20 cm), M (21-50 cm), L (>50 cm)
  bind_rows(data.frame(CommonName = c("Northern kingfish", "Buckler dory", "Fourbeard rockling",
                                      "Fourbeard rockling", "Cunner"),
                       sizecat    = c("M", "S", "S", "M", "M"),
                       range      = c("21-50", "≤20", "≤20", "21-50", "21-50"))) |>
  mutate(COMNAME = toupper(CommonName))

# benthivores = everything that is not a piscivore or planktivore
benthivores <- alldf |>
  anti_join(bind_rows(piscdf, plankdf), by = c("COMNAME", "SizeCat")) |>
  left_join(sizedef, join_by(COMNAME == COMNAME, SizeCat == sizecat)) |>
  left_join(ecodata::species_groupings, by = "COMNAME") |>
  select(COMNAME, ITISSPP, SCINAME, SizeCat = SizeCat.x, SizeRange_cm = range) |>
  distinct()

fh_nefsc <- allfh |>
  left_join(benthivores |> mutate(is_benthivore = TRUE),
            by = c("pdcomnam" = "COMNAME", "sizecat" = "SizeCat"))

if (filter_benthivores) {
  fh_nefsc <- filter(fh_nefsc, is_benthivore %in% TRUE)
}

## ---- Prey lists: Rpath Macrobenthos and Megabenthos -----------------------

load(paths$prey)  # object is called `prey`

# adjustments from modeling group discussions:
# sand/mole crabs move from Megabenthos to Macrobenthos; scallops removed from Macrobenthos
crabs_to_macro <- c("ALBUNEA PARETII", "EMERITA TALPOIDA", "RANILIA MURICATA", "HIPPIDAE")
scallops       <- c("PECTINIDAE", "PECTINIDAE SHELL", "PECTINIDAE VISCERA")

prey_names <- list(
  macro = prey |>
    filter((RPATH == "Macrobenthos" & !PYNAM %in% scallops) |
             (RPATH == "Megabenthos"  &  PYNAM %in% crabs_to_macro)) |>
    pull(PYNAM),
  mega  = prey |>
    filter(RPATH == "Megabenthos", !PYNAM %in% crabs_to_macro) |>
    pull(PYNAM)
)

## ---- NEAMAP station dates (month/day needed for the bottom temp match) ----

neamap_dates <- read.csv(paths$neamap_dates) |>
  mutate(id = station,
         year = as.numeric(year), month = as.numeric(month), day = as.numeric(day)) |>
  select(id, year, month, day) |>
  distinct()

## ---- Build station-level dataset for one prey group -----------------------
# group = "macro" or "mega". Column names follow the legacy pattern
# (e.g. meanmacrobenpywt / meanmegabenpywt) so the VAST scripts are unchanged.

build_benthos_data <- function(group) {
  
  ben  <- paste0(group, "ben")                   # "macroben" / "megaben"
  flag <- paste0(group, "benthos")               # "macrobenthos" / "megabenthos"
  nm   <- c(macro = "mac", mega = "mgb")[[group]] # NEAMAP column code
  pywt <- paste0(ben, "pywt"); pynm <- paste0(ben, "pynam")
  sumw <- paste0("sum", ben, "pywt"); meanw <- paste0("mean", ben, "pywt")
  
  # NEFSC, prey level. Season: Spring = Jan-Jun, Fall = Jul-Dec
  all_stn <- fh_nefsc |>
    mutate(!!flag := if_else(pynam %in% prey_names[[group]], ben, "othprey"),
           id = paste0(cruise6, "_", station),
           year = as.numeric(year), month = as.numeric(month), day = as.numeric(day),
           season_ng = case_when(month <= 6 ~ "SPRING", month >= 7 ~ "FALL")) |>
    select(year, month, day, season_ng, id, stratum,
           pynam, pyamtw, pywgti, pyvoli, all_of(flag),
           pdcomnam, pdid, pdlen, pdsvol, pdswgt,
           beglat, beglon, declat, declon,
           bottemp, surftemp, setdepth) |>
    group_by(id) |>
    mutate(!!pywt := if_else(.data[[flag]] == ben, pyamtw, 0),
           !!pynm := if_else(.data[[flag]] == ben, pynam, NA_character_))
  
  saveRDS(all_stn, file.path(paths$out_dir, paste0(ben, "all_stn.rds")))
  
  # NEFSC, station level: mean benthos g per stomach per tow + predator covariates
  stndat <- all_stn |>
    select(year, month, day, season_ng, id, beglat, beglon, declat, declon,
           bottemp, surftemp, setdepth) |>
    distinct()
  
  nstom <- all_stn |>
    group_by(id, pdcomnam) |>
    summarise(nstompd = n_distinct(pdid), .groups = "drop_last") |>
    summarise(nstomtot = sum(nstompd))
  
  predlen <- all_stn |>
    summarise(meanbenthivorelen = mean(pdlen), varbenthivorelen = var(pdlen))
  
  nefsc <- all_stn |>
    summarise(!!sumw := sum(.data[[pywt]]),
              !!paste0("n", ben, "pysp") := n_distinct(.data[[pynm]], na.rm = TRUE),
              npreysp = n_distinct(pynam),
              nbenthivoresp = n_distinct(pdcomnam)) |>
    left_join(nstom, by = "id") |>
    mutate(!!meanw := .data[[sumw]] / nstomtot) |>
    left_join(predlen, by = "id") |>
    left_join(stndat, by = "id") |>
    mutate(declon = -declon,
           vessel = case_when(year < 2009 ~ "AL", year >= 2009 ~ "HB"))
  
  # NEAMAP arrives aggregated; rename to match NEFSC
  neamap_names <- c(
    id = "station", nbenthivoresp = "nbenthivorespp",
    meanbenthivorelen = "meanbenthivorelen.simple", season_ng = "season",
    declat = "lat", declon = "lon", bottemp = "bWT", surftemp = "SST", setdepth = "depthm",
    setNames(paste0("sum",  nm, "preywt"), sumw),
    setNames(paste0("n",    nm, "preysp"), paste0("n", ben, "pysp")),
    setNames(paste0("mean", nm, "preywt"), meanw)
  )
  
  neamap <- readr::read_csv(paths$neamap[[group]], show_col_types = FALSE) |>
    mutate(vessel = "NEAMAP") |>
    rename(all_of(neamap_names)) |>
    left_join(neamap_dates, by = c("id", "year"))
  
  agg_all <- bind_rows(nefsc, neamap)
  saveRDS(agg_all, file.path(paths$out_dir, paste0(ben, "agg_stn_all.rds")))
  
  agg_all
}

agg <- list(macro = build_benthos_data("macro"),
            mega  = build_benthos_data("mega"))

## ---- Bottom temperature: nearest GLORYS cell on the survey date -----------
# Done once for all stations (macro and mega share stations), then joined to each.

stations <- bind_rows(agg) |>
  mutate(numdate = as.numeric(as.Date(paste0(year, "-", month, "-", day)))) |>
  select(id, declon, declat, year, numdate) |>
  na.omit() |>
  distinct() |>
  sf::st_as_sf(coords = c("declon", "declat"), crs = 4326, remove = FALSE)

if (anyDuplicated(stations$id)) stop("Station ids with conflicting location/date between macro and mega")

bt_files <- tibble::tibble(
  path = c(list.files(file.path(paths$glorys_dir, "V1"), "\\.rds$", full.names = TRUE),
           list.files(file.path(paths$glorys_dir, "V2"), "\\.rds$", full.names = TRUE))) |>
  mutate(version = basename(dirname(path)),
         year    = as.numeric(sub(".*_(\\d{4})\\.rds$", "\\1", basename(path)))) |>
  filter((version == "V1" & year <  glorys_v2_start) |
           (version == "V2" & year >= glorys_v2_start),
         year %in% stations$year)

if (!data_year %in% bt_files$year) stop("No GLORYS bottom temp file for ", data_year)

match_bt_year <- function(path, yr) {
  stn_yr <- filter(stations, year == yr)
  bt <- readRDS(path) |>
    mutate(numdate = as.numeric(date)) |>
    filter(numdate %in% stn_yr$numdate) |>
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
  
  # one nearest-neighbour join per survey day (same result as per station, far fewer calls)
  purrr::map(unique(stn_yr$numdate), \(d) {
    bt_d <- bt[bt$numdate == d, "mod_bt"]
    if (nrow(bt_d) == 0) return(NULL)
    sf::st_join(stn_yr[stn_yr$numdate == d, ], bt_d,
                join = nngeo::st_nn, k = 1, progress = FALSE) |>
      sf::st_drop_geometry() |>
      select(id, mod_bt)
  }) |>
    purrr::list_rbind()
}

mod_bt <- purrr::map2(bt_files$path, bt_files$year, match_bt_year) |>
  purrr::list_rbind()

## ---- Save VAST inputs and print checks ------------------------------------

for (group in names(agg)) {
  out <- left_join(agg[[group]], mod_bt, by = "id")
  saveRDS(out, file.path(paths$out_dir, paste0(group, "benagg_stn_all_modBT.rds")))
  
  message("\n", group, "benthos: stations by survey/season in recent years")
  print(count(out, vessel, season_ng, year) |> filter(year >= data_year - 2), n = Inf)
  
  message(group, "benthos: bottom temp coverage")
  print(out |>
          group_by(year) |>
          summarise(n = n(),
                    pct_insitu = mean(!is.na(bottemp) & bottemp != 0),
                    pct_glorys = mean(!is.na(mod_bt))) |>
          tail(5))
  
  if (!data_year %in% out$year) warning(group, "benthos has no data for ", data_year)
}