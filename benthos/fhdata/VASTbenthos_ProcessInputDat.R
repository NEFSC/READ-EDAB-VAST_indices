# Streamlined version of CreateVASTInputs.Rmd for operational updates to forage index indicators
# Modified for benthos index
# October 2024
#   This one is updating with 2025 NEFSC and NEAMAP data and OISST
#   To be used in the 2025 State of the Ecosystem report

library(tidyverse)
library(here)
library(dendextend)

# Load NEFSC stomach data received from Brian Smith

# object is called `allfh`

load(here("2025_shared_data/allfh.rmd.epu.RData"))

# load(url("https://github.com/NOAA-EDAB/forageindex/raw/main/fhdat/allfh.rmd.epu.RData"))

# # object is called `allfh`
# load(url("https://github.com/NOAA-EDAB/forageindex/raw/main/fhdat/allfh.Rdata"))
# 
# #object is called allfh21
# load(url("https://github.com/NOAA-EDAB/forageindex/raw/main/fhdat/allfh21.RData"))
# 
# #object is called allfh22
# load(url("https://github.com/NOAA-EDAB/forageindex/raw/main/fhdat/allfh22.RData"))
# 
# # bind all NEFSC stomach datasets
# allfh <- allfh %>%
#   dplyr::bind_rows(allfh21) |>
#   dplyr::bind_rows(allfh22)
# 
# # can we deload the 21 and 22 datasets?
# rm("allfh21", "allfh22")

###############################################################################
# read predator similarity info to generate predator list
# Input NEFSC food habits overlap matrix:

dietoverlap <- read_csv(here("benthos/data-raw/tgmat.2022-02-15.csv"))

# use dendextend functions to get list
d_dietoverlap <- dist(dietoverlap)

guilds <- hclust(d_dietoverlap)

#plot(guilds)

dend <- as.dendrogram(guilds)

dend <- rotate(dend, 1:136)

dend <- color_branches(dend, k=6) # Brian uses 6 categories

labels(dend) <- paste(as.character(names(dietoverlap[-1]))[order.dendrogram(dend)],
                      "(",labels(dend),")", 
                      sep = "")

dend <- hang.dendrogram(dend,hang_height=0.1)

pisc <- partition_leaves(dend)[[
  which_node(dend, c("Bluefish..S(37)", "Bluefish..M(36)", "Bluefish..L(35)"))
]]

piscdf <- data.frame("COMNAME" = toupper(str_remove(pisc, "\\..*")),
                     "SizeCat" = str_remove(str_extract(pisc, "\\..*[:upper:]+"), "\\.."))


plank <- partition_leaves(dend)[[
  which_node(dend, c("Atlantic herring..S(20)", "Atlantic mackerel..S(23)", "Blueback herring..XS(34)"))
]]

plankdf <- data.frame("COMNAME" = toupper(str_remove(plank, "\\..*")),
                      "SizeCat" = str_remove(str_extract(plank, "\\..*[:upper:]+"), "\\.."))

all <- partition_leaves(dend)[[
  which_node(dend, c("Barndoor skate..S(26)", "Bluefish..L(35)"))
]]

alldf <- data.frame("COMNAME" = toupper(str_remove(all, "\\..*")),
                    "SizeCat" = str_remove(str_extract(all, "\\..*[:upper:]+"), "\\.."))

sizedef <- readxl::read_excel(here::here("benthos/data-raw/Table3.xlsx"),
                              range = "B3:H55") |>  # table of NEFSC FH size categories from Brian
  dplyr::rename(XS = `Extra-Small`,
                S = "Small",
                M = "Medium",
                L = "Large",
                XL = `Extra-Large`,
                CommonName = `Common Name`,
                SpeciesName = `Species Name`) |>
  tidyr::pivot_longer(-c(1:2), names_to = "sizecat", values_to = "range") |>
  dplyr::filter(!range=="-") |>
  dplyr::mutate(COMNAME = toupper(CommonName)) 


# add missing by hand from word Table 1.docx 
# (pollock was misspelled, fixed in spreadsheet)
# (removed "flounder" from Windowpane in spreadsheet)
# asked Brian for these; response:
# For these predators and most likely any others not listed in TM216, 
# the default is S (<=20 cm), M (>20 and <=50 cm), and L (>50 cm).
# Northern kingfish M
# Buckler dory S
# Fourbeard rockling S
# Fourbeard rockling M
# Cunner M

extra <- data.frame(CommonName = c("Northern kingfish",
                                   "Buckler dory",
                                   "Fourbeard rockling",
                                   "Fourbeard rockling",
                                   "Cunner"),  
                    sizecat = c("M",
                                "S",
                                "S",
                                "M",
                                "M"), 
                    range = c("21-50",
                              "≤20",
                              "≤20",
                              "21-50",
                              "21-50"))  |>
  dplyr::mutate(COMNAME = toupper(CommonName)) 

sizedef <- bind_rows(sizedef, extra)

# dplyr::mutate(COMNAME = toupper(`Common Name`),
#               XS = readr::parse_number(`Extra-Small`),
#               S = readr::parse_number(Small),
#               M = readr::parse_number(Medium),
#               L = readr::parse_number(Large),
#               XL = readr::parse_number(`Extra-Large`))


# we'll call benthivores everything but piscivores and planktivores
benthivores <- alldf |>
  dplyr::anti_join(dplyr::bind_rows(piscdf, plankdf)) |>
  dplyr::left_join(sizedef, join_by(COMNAME == COMNAME, SizeCat == sizecat)) |>
  dplyr::left_join(ecodata::species_groupings, by = "COMNAME") |>
  dplyr::select(COMNAME, ITISSPP, SCINAME, SizeCat = SizeCat.x, SizeRange_cm = range) |>
  dplyr::distinct()

fh.nefsc.benthivore.complete <- allfh %>%
  #filter(pynam != "EMPTY") %>%
  left_join(benthivores, by = c("pdcomnam" = "COMNAME",
                                   "sizecat" = "SizeCat")) 


##############################################################################
# Get prey list from Rpath megabenthos and macrobenthos categories

# Rpath prey list, object is called prey
load(here("benthos/data/prey.RData")) #original mappings from June 2023

# make adjustments to list based on modeling group discussions
addtomacro <- prey |>
  dplyr::filter(RPATH == "Megabenthos") |> 
  dplyr::filter(PYNAM %in% c("ALBUNEA PARETII", "EMERITA TALPOIDA", "RANILIA MURICATA", "HIPPIDAE"))

megaben <- prey |>
  dplyr::filter(RPATH == "Megabenthos") |>
  dplyr::filter(!PYNAM %in% c("ALBUNEA PARETII", "EMERITA TALPOIDA", "RANILIA MURICATA", "HIPPIDAE"))

macroben <- prey |>
  dplyr::filter(RPATH == "Macrobenthos") |> 
  dplyr::filter(!PYNAM %in% c("PECTINIDAE", "PECTINIDAE SHELL", "PECTINIDAE VISCERA")) |>
  dplyr::bind_rows(addtomacro)

macrobenfh <- allfh %>%
  dplyr::left_join(macroben %>% 
                     dplyr::select(PYNAM, RPATH) %>%
                     setNames(., tolower(names(.)))) %>%
  dplyr::filter(!is.na(rpath))

megabenfh <- allfh %>%
  dplyr::left_join(megaben %>% 
                     dplyr::select(PYNAM, RPATH) %>%
                     setNames(., tolower(names(.)))) %>%
  dplyr::filter(!is.na(rpath))

macrobenITIS <- macroben |>
  dplyr::select(PYNAM, PYCOMNAM, PYSPP) |>
  dplyr::left_join(ecodata::species_groupings, join_by(PYNAM == SCINAME)) |>
  dplyr::select(PYNAM, PYCOMNAM, COMNAME, PYSPP, ITISSPP) |>
  as.data.frame()

megabenITIS <- megaben |>
  dplyr::select(PYNAM, PYCOMNAM, PYSPP) |>
  dplyr::left_join(ecodata::species_groupings, join_by(PYNAM == SCINAME)) |>
  dplyr::select(PYNAM, PYCOMNAM, COMNAME, PYSPP, ITISSPP) |>
  as.data.frame()


fh.nefsc.benthivore.complete.macrobenthos <- fh.nefsc.benthivore.complete %>%
  mutate(macrobenthos = case_when(pynam %in% macrobenITIS$PYNAM ~ "macroben",
                              TRUE ~ "othprey"))

fh.nefsc.benthivore.complete.megabenthos <- fh.nefsc.benthivore.complete %>%
  mutate(megabenthos = case_when(pynam %in% megabenITIS$PYNAM ~ "megaben",
                                  TRUE ~ "othprey"))


###############################################################################
# Make the NEFSC macrobenthos dataset aggregating prey based on prey list
# lets keep the month and day info for the merge with modeled bottom temperature!

macrobenall_stn <- fh.nefsc.benthivore.complete.macrobenthos %>%
  #create id linking cruise6_station
  #create season_ng spring and fall Spring=Jan-May, Fall=June-Dec
  mutate(id = paste0(cruise6, "_", station),
         year = as.numeric(year),
         month = as.numeric(month),
         day = as.numeric(day),
         season_ng = case_when(month <= 6 ~ "SPRING",
                               month >= 7 ~ "FALL",
                               TRUE ~ as.character(NA))
  ) %>%
  dplyr::select(year, month, day, season_ng, id, stratum,
                pynam, pyamtw, pywgti, pyvoli, macrobenthos, 
                pdcomnam, pdid, pdlen, pdsvol, pdswgt, 
                beglat, beglon, declat, declon, 
                bottemp, surftemp, setdepth) %>%
  group_by(id) %>%
  #mean macrobenthos g per stomach per tow: sum all macrobenthos g/n stomachs in tow
  mutate(macrobenpywt = case_when(macrobenthos == "macroben" ~ pyamtw,
                              TRUE ~ 0.0),
         macrobenpynam = case_when(macrobenthos == "macroben" ~ pynam,
                               TRUE ~ NA_character_)) 

# Optional: save at prey disaggregated stage for paper
saveRDS(macrobenall_stn, here("benthos/fhdata/macrobenall_stn.rds"))

# Now get station data in one line
stndat <- macrobenall_stn %>%
  dplyr::select(year, month, day, season_ng, id, 
                beglat, beglon, declat, declon, 
                bottemp, surftemp, setdepth) %>%
  distinct()

#benthivore stomachs in tow count pdid for each pred and sum
benthivorestom <- macrobenall_stn %>%
  group_by(id, pdcomnam) %>%
  summarise(nstompd = n_distinct(pdid)) %>%
  group_by(id) %>%
  summarise(nstomtot = sum(nstompd))

#mean and var pred length per tow
benthivorelen <-  macrobenall_stn %>%
  summarise(meanbenthivorelen = mean(pdlen),
            varbenthivorelen = var(pdlen))

# Aggregated prey at station level with predator covariates
macrobenagg_stn <- macrobenall_stn %>%
  summarise(summacrobenpywt = sum(macrobenpywt),
            nmacrobenpysp = n_distinct(macrobenpynam, na.rm = T),
            npreysp = n_distinct(pynam),
            nbenthivoresp = n_distinct(pdcomnam)) %>%
  left_join(benthivorestom) %>%
  mutate(meanmacrobenpywt = summacrobenpywt/nstomtot) %>%
  left_join(benthivorelen) %>%
  left_join(stndat)

# save at same stage as before, writing over old file
#saveRDS(bluepyagg_stn, here("fhdat/bluepyagg_stn.rds"))

# current dataset, fix declon, add vessel, rename NEFSC
#nefsc_bluepyagg_stn <- readRDS(here("fhdat/bluepyagg_stn.rds")) %>%
nefsc_macrobenagg_stn <- macrobenagg_stn %>%
  mutate(declon = -declon,
         vessel = case_when(year<2009 ~ "AL",
                            year>=2009 ~ "HB", 
                            TRUE ~ as.character(NA)))

##############################################################################
# Add NEAMAP macrobenthos to make full aggregated stomach dataset

# Read in NEAMAP updated input from Jim Gartland, reformat with same names
neamap_macrobenthos_stn <- read_csv(here("2025_shared_data/NEAMAP_Mean Stomach Weights_Macrobenthos prey_wWQ_sept2025.csv")) %>%  
  mutate(vessel = "NEAMAP") %>%
  rename(id = station,
         summacrobenpywt = summacpreywt,
         nmacrobenpysp = nmacpreysp,
         #npreysp = ,
         nbenthivoresp = nbenthivorespp,
         nstomtot = nstomtot, 
         meanmacrobenpywt = meanmacpreywt,
         meanbenthivorelen = meanbenthivorelen.simple, 
         #varbenthivorelen = ,
         season_ng = season,
         declat  = lat,
         declon = lon,
         bottemp = bWT,
         surftemp = SST, # new NEAMAP already contains SST
         setdepth = depthm) 

# Add NEAMAP month and day information

NEAMAPstationSST <- read.csv(here("2025_shared_data/NEAMAP SST_2007_2024.csv"))

# NEAMAPstationSST <- read.csv("https://raw.githubusercontent.com/NOAA-EDAB/forageindex/main/fhdat/NEAMAP%20SST_2007_2022.csv")
# 
# NEAMAPstationSST23 <- read.csv("https://github.com/NOAA-EDAB/forageindex/raw/refs/heads/main/fhdat/NEAMAP%20SST_2023.csv")
# 
# NEAMAPstationSST <- dplyr::bind_rows(NEAMAPstationSST, NEAMAPstationSST23)


NEAMAPstations <- NEAMAPstationSST %>%
  dplyr::mutate(id = station,
                year = as.numeric(year),
                month = as.numeric(month),
                day = as.numeric(day)
                ) %>%
  dplyr::select(id, year, month, day) |>
  dplyr::distinct()

neamap_macrobenthos_stn <- dplyr::left_join(neamap_macrobenthos_stn, NEAMAPstations)


# combine NEAMAP and NEFSC
macrobenagg_stn_all <-  nefsc_macrobenagg_stn %>%
  bind_rows(neamap_macrobenthos_stn) 

# Save before SST integration step
saveRDS(macrobenagg_stn_all, here("benthos/fhdata/macrobenagg_stn_all.rds"))

###############################################################################
# Now make the megabenthos dataset using the same steps, NEFSC then NEAMAP
# These will be SEPARATE UNIVARIATE MODELS to start so making them separate

###############################################################################
# Make the NEFSC macrobenthos dataset aggregating prey based on prey list

megabenall_stn <- fh.nefsc.benthivore.complete.megabenthos %>%
  #create id linking cruise6_station
  #create season_ng spring and fall Spring=Jan-May, Fall=June-Dec
  mutate(id = paste0(cruise6, "_", station),
         year = as.numeric(year),
         month = as.numeric(month),
         day = as.numeric(day),
         season_ng = case_when(month <= 6 ~ "SPRING",
                               month >= 7 ~ "FALL",
                               TRUE ~ as.character(NA))
  ) %>%
  dplyr::select(year, month, day, season_ng, id, stratum,
                pynam, pyamtw, pywgti, pyvoli, megabenthos, 
                pdcomnam, pdid, pdlen, pdsvol, pdswgt, 
                beglat, beglon, declat, declon, 
                bottemp, surftemp, setdepth) %>%
  group_by(id) %>%
  #mean megabenthos g per stomach per tow: sum all megabenthos g/n stomachs in tow
  mutate(megabenpywt = case_when(megabenthos == "megaben" ~ pyamtw,
                                  TRUE ~ 0.0),
         megabenpynam = case_when(megabenthos == "megaben" ~ pynam,
                                   TRUE ~ NA_character_)) 

# Optional: save at prey disaggregated stage for paper
saveRDS(megabenall_stn, here("benthos/fhdata/megabenall_stn.rds"))

# Now get station data in one line
stndat <- megabenall_stn %>%
  dplyr::select(year, month, day, season_ng, id, 
                beglat, beglon, declat, declon, 
                bottemp, surftemp, setdepth) %>%
  distinct()

#benthivore stomachs in tow count pdid for each pred and sum
benthivorestom <- megabenall_stn %>%
  group_by(id, pdcomnam) %>%
  summarise(nstompd = n_distinct(pdid)) %>%
  group_by(id) %>%
  summarise(nstomtot = sum(nstompd))

#mean and var pred length per tow
benthivorelen <-  megabenall_stn %>%
  summarise(meanbenthivorelen = mean(pdlen),
            varbenthivorelen = var(pdlen))

# Aggregated prey at station level with predator covariates
megabenagg_stn <- megabenall_stn %>%
  summarise(summegabenpywt = sum(megabenpywt),
            nmegabenpysp = n_distinct(megabenpynam, na.rm = T),
            npreysp = n_distinct(pynam),
            nbenthivoresp = n_distinct(pdcomnam)) %>%
  left_join(benthivorestom) %>%
  mutate(meanmegabenpywt = summegabenpywt/nstomtot) %>%
  left_join(benthivorelen) %>%
  left_join(stndat)

# save at same stage as before, writing over old file
#saveRDS(bluepyagg_stn, here("fhdat/bluepyagg_stn.rds"))

# current dataset, fix declon, add vessel, rename NEFSC
#nefsc_bluepyagg_stn <- readRDS(here("fhdat/bluepyagg_stn.rds")) %>%
nefsc_megabenagg_stn <- megabenagg_stn %>%
  mutate(declon = -declon,
         vessel = case_when(year<2009 ~ "AL",
                            year>=2009 ~ "HB", 
                            TRUE ~ as.character(NA)))

##############################################################################
# Add NEAMAP megabenthos to make full aggregated stomach dataset

# Read in NEAMAP updated input from Jim Gartland, reformat with same names
neamap_megabenthos_stn <- read_csv(here("2025_shared_data/NEAMAP_Mean Stomach Weights_Megabenthos prey_wWQ_sep2025.csv")) %>%  
  mutate(vessel = "NEAMAP") %>%
  rename(id = station,
         summegabenpywt = summgbpreywt,
         nmegabenpysp = nmgbpreysp,
         #npreysp = ,
         nbenthivoresp = nbenthivorespp,
         nstomtot = nstomtot, 
         meanmegabenpywt = meanmgbpreywt,
         meanbenthivorelen = meanbenthivorelen.simple, 
         #varbenthivorelen = ,
         season_ng = season,
         declat  = lat,
         declon = lon,
         bottemp = bWT,
         surftemp = SST, # new NEAMAP already contains SST
         setdepth = depthm) 

# Add NEAMAP month and day information
# done above but uncomment if running separately
# NEAMAPstationSST <- read.csv("https://raw.githubusercontent.com/NOAA-EDAB/forageindex/main/fhdat/NEAMAP%20SST_2007_2022.csv")
# 
# NEAMAPstations <- NEAMAPstationSST %>%
#   dplyr::mutate(id = station,
#                 year = as.numeric(year),
#                 month = as.numeric(month),
#                 day = as.numeric(day)
#   ) %>%
#   dplyr::select(id, year, month, day) |>
#   dplyr::distinct()

neamap_megabenthos_stn <- dplyr::left_join(neamap_megabenthos_stn, NEAMAPstations)

# combine NEAMAP and NEFSC
megabenagg_stn_all <-  nefsc_megabenagg_stn %>%
  bind_rows(neamap_megabenthos_stn) 

# Save before SST integration step
saveRDS(megabenagg_stn_all, here("benthos/fhdata/megabenagg_stn_all.rds"))

###############################################################################

# STOP HERE for initial models, explore without modeled temp covariates

# ###############################################################################

# functions that read in bottom temperature nc files are in 
# https://noaa-edab.github.io/benthosindex/BottomTempFill.html
# assuming that those files exist and are in the folder data-raw/bottomtemp/bt_data
# the following will merge them with the data above

library(sf)
library(raster)
library(terra)
library(nngeo)

#######################################################################
# Macrobenthos

macrobenagg_stn_all <- readRDS(here("benthos/fhdata/macrobenagg_stn_all.rds"))

stations <- macrobenagg_stn_all %>%
  #dplyr::mutate(day = str_pad(day, 2, pad='0'),
  #              month = str_pad(month, 2, pad='0'),
  #              yrmody = as.numeric(paste0(year, month, day))) %>%
  # consider this instead, may not need to pad strings? already date fields in the bt data
  dplyr::mutate(date = as.Date(paste0(year,"-", month,"-", day)), numdate = as.numeric(date)) |>
  dplyr::select(id, declon, declat, year, numdate) %>%
  na.omit() %>%
  sf::st_as_sf(coords=c("declon","declat"), crs=4326, remove=FALSE)



#list of SST dataframes
# BTdfs <- list.files(here("benthos/data-raw/bottomtemp/bt_data/"), pattern = "*.rds")

# 12/05/2025. Joe created a workflow to update GLORYS data. Reading in from there now -MTG
# 
# BTdfs <- list.files(("~/EDAB_Datasets/GLORYS/glorys_bottomT/annual_bottomT_vast/"), pattern = "*.rds")

# 12/17/2025 Joe created sub folder V2 to fix an error with years 1993-2025
# making sure the V2 files are used for those years

# Base directories
base_dir <- "~/EDAB_Datasets/GLORYS/glorys_bottomT/annual_bottomT_vast/"
v1_dir <- file.path(base_dir, "V1")
v2_dir   <- file.path(base_dir, "V2")

# List files
BTdfs_v1 <- list.files(v1_dir, pattern = "\\.rds$")
BTdfs_v2   <- list.files(v2_dir, pattern = "\\.rds$")

# Helper function to extract year from either filename style
get_year <- function(x) {
  as.numeric(sub(".*_(\\d{4})\\.rds$", "\\1", x))
}

# Extract years
v1_years <- get_year(BTdfs_v1)
v2_years   <- get_year(BTdfs_v2)

# Keep original data only up to 1992
BTdfs_keep_orig <- BTdfs_v1[v1_years <= 1992]

# Use V2 data for 1993–2025
BTdfs_use_v2 <- BTdfs_v2[v2_years >= 1993 & v2_years <= 2025]

# Combine and sort by year
BTdfs <- c(BTdfs_keep_orig, BTdfs_use_v2)
BTdfs <- BTdfs[order(get_year(BTdfs))]


dietstn_mod_bt <- tibble()

# rewriting the for loop to use the V2 subfolder when needed


for (df in BTdfs) {
  
  yr <- get_year(df)
  
  # choose correct file location
  btdf <- readRDS(
    if (yr >= 1993 && yr <= 2024) {
      file.path(v2_dir, df)
    } else {
      file.path(v1_dir, df)
    }
  )
  
  if (!yr %in% unique(stations$year)) next
  
  # keep only bluefish dates in SST year
  stationsyr <- stations %>%
    dplyr::filter(year == yr)
  
  if (nrow(stationsyr) == 0) next
  
  # keep only modeled bt days in survey dataset
  btdf_survdays <- btdf %>%
    dplyr::mutate(numdate = as.numeric(date)) %>%
    dplyr::filter(numdate %in% stationsyr$numdate) %>%
    dplyr::mutate(
      year   = as.numeric(year),
      month = as.numeric(month),
      day   = as.numeric(day),
      declon = longitude,
      declat = latitude
    ) %>%
    dplyr::select(-longitude, -latitude) %>%
    sf::st_as_sf(coords = c("declon", "declat"),
                 crs = 4326, remove = FALSE)
  
  if (nrow(btdf_survdays) == 0) next
  
  yrdietmodBT <- do.call(
    rbind,
    lapply(split(stationsyr, seq_len(nrow(stationsyr))), function(x) {
      
      y <- btdf_survdays[btdf_survdays$numdate == unique(x$numdate), ]
      if (nrow(y) == 0) return(NULL)
      
      sf::st_join(
        x, y,
        join = st_nn,
        k = 1,
        progress = FALSE
      )
    })
  )
  
  dietstn_mod_bt <- rbind(dietstn_mod_bt, yrdietmodBT)
}

# old version below

# for(df in BTdfs){
#   btdf <- readRDS(paste0("~/EDAB_Datasets/GLORYS/glorys_bottomT/annual_bottomT_vast/", df))
#   
#   if(unique(btdf$year) %in% unique(stations$year)){
#     # keep only bluefish dates in SST year
#     stationsyr <- stations %>%
#       filter(year == unique(btdf$year))
#     
#     # keep only modeled bt days in survey dataset
#     btdf_survdays <- btdf %>%
#       dplyr::mutate(numdate = as.numeric(date))%>%
#       dplyr::filter(numdate %in% unique(stationsyr$numdate)) %>%
#       dplyr::mutate(year = as.numeric(year),
#                     month = as.numeric(month),
#                     day = as.numeric(day),
#                     declon = longitude,
#                     declat = latitude) %>%
#       dplyr::select(-longitude, -latitude) %>%
#       sf::st_as_sf(coords=c("declon","declat"), crs=4326, remove=FALSE)
#     
#     # now join by nearest neighbor and date
#     
#     #https://stackoverflow.com/questions/71959927/spatial-join-two-data-frames-by-nearest-feature-and-date-in-r
#     
#     yrdietmodBT <- do.call(
#       rbind,
#       lapply(split(stationsyr, seq_len(nrow(stationsyr))), function(x) {
#         
#         y <- btdf_survdays[btdf_survdays$numdate == unique(x$numdate), ]
#         
#         if (nrow(y) == 0) {
#           return(NULL)
#         }
#         
#         sf::st_join(
#           x,
#           y,
#           join = st_nn,
#           k = 1,
#           progress = FALSE
#         )
#       })
#     )
#     
#     
#     #   #datatable solution--works but doesnt seem faster?
#     #    df1 <- data.table(stationsyr)
#     #
#     #  .nearest_samedate <- function(x) {
#     #    st_join(st_as_sf(x), sstdf_survdays[sstdf_survdays$yrmody == unique(x$yrmody),], join = st_nearest_feature)
#     #  }
#     # #
#     #  yrdietOISST <- df1[, .nearest_samedate(.SD), by = list(1:nrow(df1))]
#     
#     dietstn_mod_bt <- rbind(dietstn_mod_bt,  yrdietmodBT)
#   }
# }

#saveRDS(dietstn_OISST, here("data-raw/dietstn_OISST.rds"))

# Now join with OISST dataset

#bluepyagg_stn_all <- readRDS(here("fhdat/bluepyagg_stn_all.rds"))
#dietstn_OISST <- readRDS(here("data-raw/dietstn_OISST.rds"))


dietstn_mod_bt_merge <- dietstn_mod_bt %>%
  dplyr::rename(declon = declon.x,
                declat = declat.x,
                year = year.x) %>%
  dplyr::select(id, mod_bt) %>%
  sf::st_drop_geometry()

macrobenagg_stn_all_modBT <- left_join(macrobenagg_stn_all, dietstn_mod_bt_merge)

saveRDS(macrobenagg_stn_all_modBT, here::here("benthos/fhdata/macrobenagg_stn_all_modBT.rds"))


#######################################################################
# Megabenthos

megabenagg_stn_all <- readRDS(here("benthos/fhdata/megabenagg_stn_all.rds"))

stations <- megabenagg_stn_all %>%
  #dplyr::mutate(day = str_pad(day, 2, pad='0'),
  #              month = str_pad(month, 2, pad='0'),
  #              yrmody = as.numeric(paste0(year, month, day))) %>%
  # consider this instead, may not need to pad strings? already date fields in the bt data
  dplyr::mutate(date = as.Date(paste0(year,"-", month,"-", day)), numdate = as.numeric(date)) |>
  dplyr::select(id, declon, declat, year, numdate) %>%
  na.omit() %>%
  sf::st_as_sf(coords=c("declon","declat"), crs=4326, remove=FALSE)



#list of SST dataframes
# BTdfs <- list.files(here("benthos/data-raw/bottomtemp/bt_data/"), pattern = "*.rds")

# 12/17/2025 Joe created sub folder V2 to fix an error with years 1993-2025
# making sure the V2 files are used for those years

# Base directories
base_dir <- "~/EDAB_Datasets/GLORYS/glorys_bottomT/annual_bottomT_vast/"
v1_dir <- file.path(base_dir, "V1")
v2_dir   <- file.path(base_dir, "V2")

# List files
BTdfs_v1 <- list.files(v1_dir, pattern = "\\.rds$")
BTdfs_v2   <- list.files(v2_dir, pattern = "\\.rds$")

# Helper function to extract year from either filename style
get_year <- function(x) {
  as.numeric(sub(".*_(\\d{4})\\.rds$", "\\1", x))
}

# Extract years
v1_years <- get_year(BTdfs_v1)
v2_years   <- get_year(BTdfs_v2)

# Keep original data only up to 1992
BTdfs_keep_orig <- BTdfs_v1[v1_years <= 1992]

# Use V2 data for 1993–2025
BTdfs_use_v2 <- BTdfs_v2[v2_years >= 1993 & v2_years <= 2025]

# Combine and sort by year
BTdfs <- c(BTdfs_keep_orig, BTdfs_use_v2)
BTdfs <- BTdfs[order(get_year(BTdfs))]


dietstn_mod_bt <- tibble()

# rewriting the for loop to use the V2 subfolder when needed


for (df in BTdfs) {
  
  yr <- get_year(df)
  
  # choose correct file location
  btdf <- readRDS(
    if (yr >= 1993 && yr <= 2024) {
      file.path(v2_dir, df)
    } else {
      file.path(v1_dir, df)
    }
  )
  
  if (!yr %in% unique(stations$year)) next
  
  # keep only bluefish dates in SST year
  stationsyr <- stations %>%
    dplyr::filter(year == yr)
  
  if (nrow(stationsyr) == 0) next
  
  # keep only modeled bt days in survey dataset
  btdf_survdays <- btdf %>%
    dplyr::mutate(numdate = as.numeric(date)) %>%
    dplyr::filter(numdate %in% stationsyr$numdate) %>%
    dplyr::mutate(
      year   = as.numeric(year),
      month = as.numeric(month),
      day   = as.numeric(day),
      declon = longitude,
      declat = latitude
    ) %>%
    dplyr::select(-longitude, -latitude) %>%
    sf::st_as_sf(coords = c("declon", "declat"),
                 crs = 4326, remove = FALSE)
  
  if (nrow(btdf_survdays) == 0) next
  
  yrdietmodBT <- do.call(
    rbind,
    lapply(split(stationsyr, seq_len(nrow(stationsyr))), function(x) {
      
      y <- btdf_survdays[btdf_survdays$numdate == unique(x$numdate), ]
      if (nrow(y) == 0) return(NULL)
      
      sf::st_join(
        x, y,
        join = st_nn,
        k = 1,
        progress = FALSE
      )
    })
  )
  
  dietstn_mod_bt <- rbind(dietstn_mod_bt, yrdietmodBT)
}

# old version below

# for(df in BTdfs){
#   btdf <- readRDS(paste0(here("benthos/data-raw/bottomtemp/bt_data/", df)))
#   
#   if(unique(btdf$year) %in% unique(stations$year)){
#     # keep only bluefish dates in SST year
#     stationsyr <- stations %>%
#       filter(year == unique(btdf$year))
#     
#     # keep only modeled bt days in survey dataset
#     btdf_survdays <- btdf %>%
#       dplyr::mutate(numdate = as.numeric(date))%>%
#       dplyr::filter(numdate %in% unique(stationsyr$numdate)) %>%
#       dplyr::mutate(year = as.numeric(year),
#                     month = as.numeric(month),
#                     day = as.numeric(day),
#                     declon = longitude,
#                     declat = latitude) %>%
#       dplyr::select(-longitude, -latitude) %>%
#       sf::st_as_sf(coords=c("declon","declat"), crs=4326, remove=FALSE)
#     
#     # now join by nearest neighbor and date
#     
#     #https://stackoverflow.com/questions/71959927/spatial-join-two-data-frames-by-nearest-feature-and-date-in-r
#     
#     yrdietmodBT <- do.call('rbind', lapply(split(stationsyr, 1:nrow(stationsyr)), function(x) {
#       sf::st_join(x, btdf_survdays[btdf_survdays$numdate == unique(x$ numdate),],
#                   #join = st_nearest_feature
#                   join = st_nn, k = 1, progress = FALSE
#       )
#     }))
#     
#     #   #datatable solution--works but doesnt seem faster?
#     #    df1 <- data.table(stationsyr)
#     #
#     #  .nearest_samedate <- function(x) {
#     #    st_join(st_as_sf(x), sstdf_survdays[sstdf_survdays$yrmody == unique(x$yrmody),], join = st_nearest_feature)
#     #  }
#     # #
#     #  yrdietOISST <- df1[, .nearest_samedate(.SD), by = list(1:nrow(df1))]
#     
#     dietstn_mod_bt <- rbind(dietstn_mod_bt,  yrdietmodBT)
#   }
# }

#saveRDS(dietstn_OISST, here("data-raw/dietstn_OISST.rds"))

# Now join with OISST dataset

#bluepyagg_stn_all <- readRDS(here("fhdat/bluepyagg_stn_all.rds"))
#dietstn_OISST <- readRDS(here("data-raw/dietstn_OISST.rds"))


dietstn_mod_bt_merge <- dietstn_mod_bt %>%
  dplyr::rename(declon = declon.x,
                declat = declat.x,
                year = year.x) %>%
  dplyr::select(id, mod_bt) %>%
  sf::st_drop_geometry()

megabenagg_stn_all_modBT <- left_join(megabenagg_stn_all, dietstn_mod_bt_merge)

saveRDS(megabenagg_stn_all_modBT, here::here("benthos/fhdata/megabenagg_stn_all_modBT.rds"))
