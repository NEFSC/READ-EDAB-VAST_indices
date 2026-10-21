# MTG 03/03/2026
# Creating this script to generate cogout.rds which is used to generate the SOE COG indicator
# Chunk comes from Results.Rmd
# That script only creates cogout.rds for macrobenthos spring
# repeating for all combinations of Rpath group and season

library(FishStatsUtils)
library(purrr)
library(here)

# get directories from pyindex ----------
outdir  <- here::here("benthos/pyindex")
moddirs <- list.dirs(outdir, recursive = FALSE, full.names = TRUE)

# generalize function from Results.Rmd -------------

getcogVAST <- function(d.name){
  
  message("Processing: ", d.name)
  
  fit <- VAST::reload_model(
    readRDS(file.path(d.name, "fit.rds"))
  )
  
  # create test directory if it doesn't exist
  dir.create(file.path(d.name, "test"), showWarnings = FALSE)
  
  cogout <- FishStatsUtils::plot_range_index(
    Sdreport      = fit$parameter_estimates$SD, 
    Report        = fit$Report, 
    TmbData       = fit$data_list,
    year_labels   = as.numeric(fit$year_labels),
    years_to_plot = fit$years_to_plot,
    Znames        = colnames(fit$data_list$Z_xm)
  )
  
  # save inside the SAME model directory
  saveRDS(
    cogout,
    file.path(d.name, "cogout.rds")
  )
}

# Run function for every combination of Rpath group and Season --------

purrr::walk(moddirs, getcogVAST)
