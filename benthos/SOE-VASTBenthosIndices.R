# create rds for ecodata input
# aim for similar structure to other ecodata datasets


library(dplyr)
library(ggplot2)
library(tidyr)

SOEinputs <- function(infile, season, taxa, outfile) {
  
  splitoutput <- read.csv(infile)
  
  # warning, hardcoded. obviously
  stratlook <- data.frame(Stratum = c("Stratum_1",
                                      "Stratum_2",
                                      "Stratum_3",
                                      "Stratum_4",
                                      "Stratum_5"),
                          Region  = c("AllEPU", 
                                      "MAB",
                                      "GB",
                                      "GOM",
                                      "SS"))
  
  benthosindex <- splitoutput %>%
    left_join(stratlook) %>%
    dplyr::select(Time, 
                  EPU = Region, 
                  "Biomass Index Estimate" = Estimate, 
                  "Biomass Index Estimate SE" = Std..Error.for.Estimate) %>%
    tidyr::pivot_longer(c("Biomass Index Estimate", "Biomass Index Estimate SE"), 
                        names_to = "Var", values_to = "Value") %>%
    dplyr::filter(EPU %in% c("MAB", "GB", "GOM", "AllEPU")) %>%
    dplyr::mutate(Units = "relative grams per stomach") %>%
    dplyr::select(Time, Var, Value, EPU, Units)
  
  benthosindex$Var <- stringr::str_c(season, taxa, benthosindex$Var, sep = " ")
  
  saveRDS(benthosindex, outfile)
  
} 


# make data files
SOEinputs(infile = "benthos/pyindex/macrobenthos_fall_500_cov/Index.csv",
          season = "Fall", 
          taxa = "Macrobenthos",
          outfile = "benthos/pyindex/fallmacrobenthosindex.rds")

SOEinputs(infile = "benthos/pyindex/macrobenthos_spring_500_cov/Index.csv",
          season = "Spring", 
          taxa = "Macrobenthos",
          outfile = "benthos/pyindex/springmacrobenthosindex.rds")

SOEinputs(infile = "benthos/pyindex/megabenthos_fall_500_cov/Index.csv",
          season = "Fall", 
          taxa = "Megabenthos",
          outfile = "benthos/pyindex/fallmegabenthosindex.rds")

SOEinputs(infile = "benthos/pyindex/megabenthos_spring_500_cov/Index.csv",
          season = "Spring", 
          taxa = "Megabenthos",
          outfile = "benthos/pyindex/springmegabenthosindex.rds")


# test plot
# foragewide <- forageindex %>%
#   pivot_wider(names_from = Var, values_from = Value)
# 
# 
# ggplot(foragewide, aes(x=Time, y=`Forage Fish Biomass Estimate`, colour = EPU)) +
#   geom_errorbar(aes(ymin=`Forage Fish Biomass Estimate`+`Forage Fish Biomass Estimate SE`, 
#                     ymax=`Forage Fish Biomass Estimate`-`Forage Fish Biomass Estimate SE`))+
#   geom_point()+
#   geom_line()

macrowide <- fallmacrobenthosindex |> tidyr::pivot_wider(names_from = Var, values_from = Value)

ggplot2::ggplot(macrowide, ggplot2::aes(x=Time, y=`Fall Macrobenthos Biomass Index Estimate`, colour=EPU)) + 
  ggplot2::geom_ribbon(ggplot2::aes(ymin=`Fall Macrobenthos Biomass Index Estimate`-`Fall Macrobenthos Biomass Index Estimate SE`, 
                                      ymax=`Fall Macrobenthos Biomass Index Estimate`+`Fall Macrobenthos Biomass Index Estimate SE`,
                                    fill=EPU),
                      alpha=0.5)+
  ggplot2::geom_point()+
  ggplot2::geom_line() #+
  #ggplot2::facet_wrap(~EPU)

