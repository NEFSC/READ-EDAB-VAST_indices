# format center of gravity output for ecodata submission and SOE

library(dplyr)
library(ggplot2)
library(tidyr)

SOEinputsCOG <- function(infile, season, taxa, outfile) {
  
  cogout <- readRDS(infile)
  
  benthoscog <- as.data.frame(cogout$COG_Table) |>
    dplyr::mutate(direction = ifelse(m==1, "Eastward", "Northward")) |>
    dplyr::select("Time" = Year,
                  "Center of Gravity" = COG_hat, 
                  "Center of Gravity SE" = SE,
                  direction) |>
    tidyr::pivot_longer(c("Center of Gravity", "Center of Gravity SE"), 
                        names_to = "Var", values_to = "Value") |>
    #direction into Var
    tidyr::unite(Var, direction:Var, sep = " ") |>
    dplyr::mutate(Units = "km",
                  EPU = "ALLEPU") |>
    dplyr::select(Time, Var, Value, EPU, Units)
  
  benthoscog$Var <- stringr::str_c(season, taxa, benthoscog$Var, sep = " ")
  
  #readr::write_csv(benthosindex, outfile)
  saveRDS(benthoscog, outfile)
  
  
}

# make data files
SOEinputsCOG(infile = "benthos/pyindex/macrobenthos_fall_500_cov/cogout.rds",
             season = "Fall", 
             taxa = "Macrobenthos",
             outfile = "benthos/pyindex/fallmacrobenthoscog.rds")

SOEinputsCOG(infile = "benthos/pyindex/macrobenthos_spring_500_cov/cogout.rds",
             season = "Spring", 
             taxa = "Macrobenthos",
             outfile = "benthos/pyindex/springmacrobenthoscog.rds")

SOEinputsCOG(infile = "benthos/pyindex/megabenthos_fall_500_cov/cogout.rds",
             season = "Fall", 
             taxa = "Megabenthos",
             outfile = "benthos/pyindex/fallmegabenthoscog.rds")

SOEinputsCOG(infile = "benthos/pyindex/megabenthos_spring_500_cov/cogout.rds",
             season = "Spring", 
             taxa = "Megabenthos",
             outfile = "benthos/pyindex/springmegabenthoscog.rds")
