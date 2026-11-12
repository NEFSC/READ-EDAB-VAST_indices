# VASTbenthos_run_models.R
#
# Fits the univariate VAST models for the benthos indicator (State of the Ecosystem):
# macrobenthos and megabenthos x fall and spring. Replaces the two legacy scripts
# VASTunivariate_{macro,mega}benthos_allsurvs.R with identical data prep and model
# settings, plus a convergence check and center-of-gravity output.
#
# modified from https://github.com/James-Thorson-NOAA/VAST/wiki/Index-standardization
# Original workflow: S. Gaichas (NOAA-EDAB/benthosindex). Maintained by M. Grezlik.
#
# Input:  benthos/fhdata/{macro,mega}benagg_stn_all_modBT.rds  (from VASTbenthos_process_input_data.R)
# Output: benthos/pyindex/{group}benthos_{season}_500_cov/
#           fit.rds (gitignored), Index.csv, cogout.rds, parameter_estimates, plots
#
# Run all four models in sequence (~96 h): source this file.
# Run one model: set `vast_jobs` before sourcing, e.g.
#   vast_jobs <- data.frame(group = "macro", season = "fall"); source(<this file>)
# Run in parallel: see the launcher at the bottom of this file.

library(dplyr)
library(VAST)

## ---- Config ----------------------------------------------------------------

first_year <- 1980    # models use years > 1979
n_knots    <- 500

# TRUE = quick check that the script runs end to end (recent years only, no bias correction,
# output to *_test folders). Takes minutes, not hours. Results are NOT for the SOE.
# Keeps 500 knots: fewer knots pushes logkappa2 onto its bound and VAST stops.
# Can also be set before sourcing: test_mode <- TRUE; source(<this file>)
if (!exists("test_mode")) test_mode <- FALSE

# "" = current version; "benthivores_only" = predator-filter version (inputs from benthos/fhdata/benthivores_only)
if (!exists("variant")) variant <- ""
in_dir  <- file.path(here::here("benthos/fhdata"),  variant)
out_dir <- file.path(here::here("benthos/pyindex"), variant)

# TRUE = center catchability covariates at their mean (index is then predicted at average
# predator length / number of species / bottom temp instead of 0). Same likelihood, better
# convergence; rescales the index by a constant only. Can be set before sourcing.
if (!exists("center_covariates")) center_covariates <- FALSE

## ---- Data prep (identical to the legacy scripts) ---------------------------

prep_vast_data <- function(group, season) {
  ben <- paste0(group, "ben")
  
  readRDS(file.path(in_dir, paste0(ben, "agg_stn_all_modBT.rds"))) |>
    # bottom temp: in situ unless missing or 0, otherwise GLORYS
    mutate(btfill = ifelse(is.na(bottemp) | bottemp == 0, mod_bt, bottemp)) |>
    filter(season_ng == toupper(season), year >= if (test_mode) 2015 else first_year) |>
    mutate(AreaSwept_km2 = 1,
           Vessel = as.numeric(as.factor(vessel)) - 1) |>   # AL = 0, HB = 1, NEAMAP = 2
    select(Catch_g = all_of(paste0("mean", ben, "pywt")),
           Year = year, Vessel, AreaSwept_km2,
           Lat = declat, Lon = declon,
           meanbenthivorelen, nbenthivoresp,
           mod_bt,   # legacy: keeping mod_bt here drops stations with in situ BT but no GLORYS match
           btfill) |>  #   (methods question in the 2027 issue; kept for consistency with 2026 SOE)
    na.omit() |>
    mutate(across(c(meanbenthivorelen, nbenthivoresp, btfill),
                  \(x) if (center_covariates) x - mean(x) else x)) |>
    as.data.frame()
}

## ---- Spatial strata: EPUs as NEFSC strata (south of Hatteras excluded) -----
# NEFSC strata limits https://github.com/James-Thorson-NOAA/VAST/issues/302

epu_strata <- list(
  MAB = c(1010:1080, 1100:1120, 1600:1750, 3010:3450, 3470, 3500, 3510),
  GB  = c(1090, 1130:1210, 1230, 1250, 3460, 3480, 3490, 3520:3550),
  GOM = c(1220, 1240, 1260:1290, 1360:1400, 3560:3830),
  SS  = c(1300:1352, 3840:3990)
)

grid_strata <- function(strata) {
  FishStatsUtils::northwest_atlantic_grid |>
    filter(stratum_number %in% strata) |>
    select(stratum_number) |>
    distinct()
}

# order matters: Stratum_1..5 in Index.csv = AllEPU, MAB, GB, GOM, SS
strata.limits <- as.list(c("AllEPU" = grid_strata(unlist(epu_strata)),
                           "MAB"    = grid_strata(epu_strata$MAB),
                           "GB"     = grid_strata(epu_strata$GB),
                           "GOM"    = grid_strata(epu_strata$GOM),
                           "SS"     = grid_strata(epu_strata$SS)))

## ---- Model settings (identical to the legacy scripts) ----------------------

# spatial, spatio-temporal and year effects on for both components
FieldConfig <- matrix("IID", ncol = 2, nrow = 3,
                      dimnames = list(c("Omega", "Epsilon", "Beta"), c("Component_1", "Component_2")))

# 0 = year intercepts and spatio-temporal fields independent among years (fixed effects)
RhoConfig <- c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 0, "Epsilon2" = 0)

# vessel effects off (catchability covariates performed better in model selection)
OverdispersionConfig <- c("eta1" = 0, "eta2" = 0)

settings <- make_settings(n_x = n_knots,
                          Region = "northwest_atlantic",
                          Version = "VAST_v14_0_1",   # needed to prevent error from newer dev version number
                          strata.limits = strata.limits,
                          purpose = "index2",
                          bias.correct = !test_mode,
                          use_anisotropy = TRUE,
                          fine_scale = TRUE,
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          OverdispersionConfig = OverdispersionConfig,
                          mesh_package = "fmesher")

## ---- Convergence check -----------------------------------------------------

check_vast_fit <- function(fit, dir, max_grad = 1e-4) {
  pe <- fit$parameter_estimates
  problems <- c(
    if (pe$max_gradient > max_grad) paste("max gradient =", signif(pe$max_gradient, 3)),
    if (!isTRUE(pe$SD$pdHess)) "Hessian not positive definite",
    if (any(is.nan(summary(pe$SD, "fixed")[, "Std. Error"]))) "NaN standard errors",
    if (!is.null(pe$diagnostics)) {
      d <- pe$diagnostics
      at_bound <- d$Param[abs(d$MLE - d$Lower) < 1e-6 | abs(d$MLE - d$Upper) < 1e-6]
      if (length(at_bound)) paste("parameters at bounds:", paste(at_bound, collapse = ", "))
    }
  )
  if (length(problems)) {
    msg <- paste0("CHECK BEFORE USING THIS FIT:\n- ", paste(problems, collapse = "\n- "))
    warning(msg, call. = FALSE)
    writeLines(msg, file.path(dir, "CONVERGENCE_WARNING.txt"))
    return(invisible(FALSE))
  }
  message("Convergence checks passed (max gradient ", signif(pe$max_gradient, 3), ")")
  invisible(TRUE)
}

## ---- Fit one model ---------------------------------------------------------

run_benthos_vast <- function(group, season) {
  dat <- prep_vast_data(group, season)
  
  working_dir <- file.path(out_dir, paste0(group, "benthos_", season, "_", n_knots, "_cov",
                                           if (test_mode) "_test"))
  dir.create(working_dir, showWarnings = FALSE, recursive = TRUE)
  wd <- paste0(working_dir, "/")
  
  message("\n=== ", group, "benthos ", season, ": ", nrow(dat), " stations, ",
          min(dat$Year), "-", max(dat$Year), " | ", format(Sys.time()), " ===")
  
  fit <- fit_model(settings = settings,
                   Lat_i = dat$Lat,
                   Lon_i = dat$Lon,
                   t_i   = dat$Year,
                   b_i   = as_units(dat[, "Catch_g"], "g"),
                   a_i   = rep(1, nrow(dat)),
                   v_i   = dat$Vessel,
                   Q_ik  = as.matrix(dat[, c("meanbenthivorelen", "nbenthivoresp", "btfill")]),
                   working_dir = wd,
                   newtonsteps = 3,                                        # fit_model argument: Newton steps after the 2nd optimizer run
                   control     = list(eval.max = 600, iter.max = 300, trace = 1))  # caps only the 2nd run, which starts at the 1st run's optimum)
  
  saveRDS(fit, file.path(working_dir, "fit.rds"))   # saved before the check so a failed fit can be inspected
  check_vast_fit(fit, working_dir)
  
  # center of gravity for the SOE (previously a separate step, get_cogout.R)
  cogout <- FishStatsUtils::plot_range_index(Sdreport      = fit$parameter_estimates$SD,
                                             Report        = fit$Report,
                                             TmbData       = fit$data_list,
                                             year_labels   = as.numeric(fit$year_labels),
                                             years_to_plot = fit$years_to_plot,
                                             Znames        = colnames(fit$data_list$Z_xm),
                                             PlotDir       = wd)
  saveRDS(cogout, file.path(working_dir, "cogout.rds"))
  
  tryCatch(plot(fit, working_dir = wd),   # Index.csv and diagnostic plots
           error = \(e) message("Plotting stopped early (outputs above are saved): ", conditionMessage(e)))
  
  message("=== done ", group, "benthos ", season, " | ", format(Sys.time()), " ===")
  invisible(fit)
}

## ---- Run -------------------------------------------------------------------

if (!exists("vast_jobs")) {
  vast_jobs <- expand.grid(group = c("macro", "mega"), season = c("fall", "spring"),
                           stringsAsFactors = FALSE)
}

for (i in seq_len(nrow(vast_jobs))) {
  run_benthos_vast(vast_jobs$group[i], vast_jobs$season[i])
}

## ---- Parallel launcher (run from the console, do not source) ---------------
# Each model runs in its own background R process with its own log file.
# Check RAM (Task Manager) after starting the first one; each fit can use several GB.
#
# launch <- function(group, season, variant = "") {
#   log <- here::here("benthos/pyindex", paste0(group, "_", season, if (nzchar(variant)) paste0("_", variant), ".log"))
#   callr::r_bg(function(g, s, v, script) {
#                 vast_jobs <- data.frame(group = g, season = s)
#                 variant <- v
#                 source(script, local = TRUE)
#               },
#               args   = list(g = group, s = season, v = variant,
#                             script = here::here("benthos/VASTscripts/VASTbenthos_run_models.R")),
#               stdout = log, stderr = "2>&1", wd = here::here())
# }