#' Process NEAMAP Trawl Data
#'
#' Reads a CSV file from the NEAMAP trawl survey, renames columns to a
#' standardized format, and adds a `vessel` identifier.
#'
#' @param filepath A character string specifying the path to the CSV file.
#'   It is recommended to use a path relative to the project root.
#'
#' @return A data frame with processed NEAMAP data. The returned data frame
#'   will have standardized column names, including `vessel`, `id`, `sumbluepywt`,
#'   `nbluepysp`, `npiscsp`, `nstomtot`, `meanbluepywt`, `meanpisclen`,
#'   `season_ng`, `declat`, `declon`, `bottemp`, and `setdepth`.
#'
#' @examples
#' \dontrun{
#' # Assuming 'neamap_data.csv' is in a 'data' directory at the project root
#' neamap_processed <- process_neamap_data("data/neamap_data.csv")
#' head(neamap_processed)
#' }

process_neamap_data <- function(filepath) {
  output <- read.csv(here::here(filepath)) %>%
    dplyr::mutate(vessel = "NEAMAP") %>%
    dplyr::rename(
      id = station,
      sumbluepywt = sumbluepreywt,
      nbluepysp = nbfpreyspp,
      npiscsp = npiscspp,
      nstomtot = nstomtot,
      meanbluepywt = meanbluepreywt,
      meanpisclen = meanpisclen.simple,
      season_ng = season,
      declat = lat,
      declon = lon,
      bottemp = bWT,
      setdepth = depthm
    )
  return(output)
}

#' @title Download and Process OISST Data to a Tidy Data Frame
#'
#' @description This function downloads OISST data for specified years,
#'   processes it into a raster, and then converts the raster into a
#'   tidy data frame, saving it as an RDS file. It includes checks to
#'   skip steps if files already exist.
#'
#' @param years A numeric vector of years to download and process.
#' @param varname A character string specifying the variable name to extract
#'   from the NetCDF file (e.g., "sst").
#' @param nc_to_raster The function to convert NetCDF files to rasters.
#' @param raster_to_sstdf The function to convert rasters to tidy data frames.
#'
#' @return The function does not return a value, but it saves processed data
#'   frames to the "data-raw/gridded/sst_data" directory.
#'
download_and_process_oisst <- function(
  years,
  varname,
  nc_to_raster,
  raster_to_sstdf,
  out_dir
) {
  # Create the target directory if it doesn't exist.
  dir.create(
    here::here(out_dir),
    recursive = TRUE,
    showWarnings = FALSE
  )

  # Increase the download timeout to prevent errors on large files.
  options(timeout = max(300, getOption("timeout")))

  for (i in years) {
    # Define file paths and URLs for the current year.
    nc_filename <- here::here(paste0(i, ".nc"))
    grd_filename <- here::here(
      out_dir,
      paste0("test_", i, ".grd")
    )
    rds_filename <- here::here(
      out_dir,
      paste0("sst", i, ".rds")
    )

    url <- paste0(
      "https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/sst.day.mean.",
      i,
      ".nc"
    )

    # Step 1: Download and process the NetCDF to a raster file.
    if (!file.exists(grd_filename)) {
      message(paste("Downloading and processing data for year:", i))
      download.file(url, destfile = nc_filename)
      temp_raster <- nc_to_raster(nc = nc_filename, varname = varname)
      raster::writeRaster(
        temp_raster,
        filename = grd_filename,
        overwrite = TRUE
      )
      unlink(nc_filename)
      message(paste("Finished processing and saving raster for", i))
    } else {
      message(paste("Raster file for", i, "already exists. Skipping download."))
    }

    # Step 2: Convert the raster to a tidy data frame and save as an RDS file.
    if (file.exists(grd_filename) && !file.exists(rds_filename)) {
      message(paste("Converting raster to data frame for year:", i))
      temp_raster <- raster::brick(grd_filename)
      sst_df <- raster_to_sstdf(brick = temp_raster)
      saveRDS(sst_df, rds_filename)
      message(paste("Converted raster to data frame and saved RDS for", i))
    } else {
      message(paste(
        "RDS file for",
        i,
        "already exists or raster not found. Skipping conversion."
      ))
    }
  }
}

process_oisst <- function(years, out_dir, varname, input) {
  # Create the target directory if it doesn't exist.
  dir.create(
    here::here(out_dir),
    recursive = TRUE,
    showWarnings = FALSE
  )

  for (i in years) {
    grd_filename <- here::here(
      out_dir,
      paste0("test_", i, ".grd")
    )
    rds_filename <- here::here(
      out_dir,
      paste0("sst", i, ".rds")
    )

    # Step 1: get NetCDF from network drive and process to a raster file.
    if (!file.exists(grd_filename)) {
      message(paste("Accessing and processing data for year:", i))

      temp_raster <- nc_to_raster(nc = input, varname = varname)
      raster::writeRaster(
        temp_raster,
        filename = grd_filename,
        overwrite = TRUE
      )
      message(paste("Finished processing and saving raster for", i))
    } else {
      message(paste("Raster file for", i, "already exists. Skipping download."))
    }

    # Step 2: Convert the raster to a tidy data frame and save as an RDS file.
    if (file.exists(grd_filename) && !file.exists(rds_filename)) {
      message(paste("Converting raster to data frame for year:", i))
      temp_raster <- raster::brick(grd_filename)
      sst_df <- raster_to_sstdf(brick = temp_raster)
      saveRDS(sst_df, rds_filename)
      message(paste("Converted raster to data frame and saved RDS for", i))
    } else {
      message(paste(
        "RDS file for",
        i,
        "already exists or raster not found. Skipping conversion."
      ))
    }
  }
}

#' Convert NetCDF Data to a Cropped RasterBrick
#'
#' Reads a specified variable from a NetCDF file, converts it into a
#' RasterBrick object, sets its spatial extent and Coordinate Reference System
#' (CRS), and then crops the data to a specified region.
#'
#' @param nc A character string representing the file path to the NetCDF (.nc) file.
#' @param varname A character string specifying the name of the variable to
#'   extract from the NetCDF file.
#' @param extent A numeric vector of four elements: c(xmin, xmax, ymin, ymax)
#'   to set the overall spatial extent of the raster. Defaults to global extent.
#' @param crop An object of class `raster::Extent` or a numeric vector of four
#'   elements to define the cropping boundaries. Defaults to a region
#'   over the North Atlantic.
#' @param show_images A logical value. If TRUE, plots of the full and cropped
#'   datasets will be displayed. Defaults to FALSE.
#'
#' @return A `RasterBrick` object containing the data from the specified variable,
#'   with its extent, CRS, and cropped to the defined area.
#'
#' @details
#' This function is designed to handle NetCDF files, a common format for storing
#' multi-dimensional scientific data. It uses the `raster` package to efficiently
#' read and process the gridded data. The default `extent` and `crop` arguments
#' are tailored for data over the North Atlantic. The CRS is set to a common
#' projection for that region.
#'
#' @examples
#' \dontrun{
#' # Assuming you have a NetCDF file named 'my_data.nc' with a variable 'temp'
#' # in the working directory.
#'
#' # Process the data and crop it to the default area
#' my_raster <- nc_to_raster("my_data.nc", "temp")
#'
#' # Process the data, crop to a custom area, and show plots
#' custom_area <- raster::extent(c(290, 310, 35, 55))
#' my_custom_raster <- nc_to_raster("my_data.nc", "temp", crop = custom_area, show_images = TRUE)
#' }
#'
#' @importFrom raster brick crs extent crop plot
#' @importFrom grDevices dev.off
#' @importFrom graphics par

# Bastille function from https://github.com/kimberly-bastille/ecopull/blob/main/R/utils.R

nc_to_raster <- function(
  nc,
  varname,
  extent = c(0, 360, -90, 90),
  crop = raster::extent(280, 300, 30, 50),
  show_images = FALSE
) {
  message("Reading .nc as brick...")

  r <- raster::brick(nc, varname = varname)

  message("Setting CRS...")
  raster::crs(
    r
  ) <- "+proj=longlat +lat_1=35 +lat_2=45 +lat_0=40 +lon_0=-77 +x_0=0 +y_0=0 +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"

  # not sure if this is necessary?
  raster::extent(r) <- raster::extent(extent)

  if (show_images) {
    par(mfrow = c(1, 2))
    raster::plot(r, 1, sub = "Full dataset")
  }

  message("Cropping data...")
  ne_data <- raster::crop(r, crop)
  #ne_data <- raster::rotate(ne_data) add here for future pulls

  if (show_images) {
    raster::plot(ne_data, 1, sub = "Cropped dataset")
    par(mfrow = c(1, 1))
  }

  message("Done!")

  return(ne_data)
}

#' Convert a RasterBrick to a Sea Surface Temperature Data Frame
#'
#' Processes a `RasterBrick` object, typically containing sea surface temperature
#' (SST) data, by optionally rotating the raster, cropping it to a specific
#' geographic extent, and then converting it into a tidy data frame.
#'
#' @param brick A `RasterBrick` object to be processed. This object should
#'   contain raster layers, with layer names in the format "XYYYY.MM.DD".
#' @param rotate A logical value. If `TRUE`, the raster will be rotated
#'   to correct for a 0-360 degree longitude representation. Defaults to `TRUE`.
#'
#' @return A data frame with four columns: `Lon`, `Lat`, and `sst`. The `sst`
#'   column contains the sea surface temperature values, while the other columns
#'   are used to represent the spatial coordinates of each temperature point.
#'
#' @details
#' This function is specifically designed to handle raster data where the layer
#' names encode the date and the longitude values are in a 0-360 degree range.
#' The function first rotates the raster to a standard -180 to 180 degree
#' longitude range if `rotate` is set to `TRUE`. It then crops the raster to a
#' fixed extent of -77 to -65 longitude and 35 to 45 latitude, corresponding to a
#' region off the US Northeast coast. Finally, it uses `rasterToPoints` to
#' convert the spatial data into a data frame and `pivot_longer` to transform it
#' into a tidy format.
#'
#' @importFrom raster rotate crop rasterToPoints extent
#' @importFrom dplyr rename
#' @importFrom tidyr pivot_longer
#'
#' @examples
#' \dontrun{
#' # Assuming 'sst_brick' is a RasterBrick object with the expected structure
#' # and layer names.
#'
#' # Convert the raster and apply the rotation
#' sst_data_df <- raster_to_sstdf(sst_brick)
#' head(sst_data_df)
#'
#' # Convert the raster without applying rotation
#' sst_data_no_rot <- raster_to_sstdf(sst_brick, rotate = FALSE)
#' }

# function to convert to dataframe based on
# https://towardsdatascience.com/transforming-spatial-data-to-tabular-data-in-r-4dab139f311f

raster_to_sstdf <- function(brick, rotate = TRUE) {
  if (rotate) {
    brick_r <- raster::rotate(brick)
  }
  brick_r <- raster::crop(brick_r, raster::extent(-77, -65, 35, 45))
  sstdf <- as.data.frame(raster::rasterToPoints(brick_r, spatial = TRUE))
  sstdf <- sstdf %>%
    dplyr::rename(Lon = x, Lat = y) %>%
    tidyr::pivot_longer(
      cols = starts_with("X"),
      names_to = c("year", "month", "day"),
      names_prefix = "X",
      names_sep = "\\.",
      values_to = "sst",
    )
  return(sstdf)
}

#' @title Join OISST Data to Survey Station Locations
#'
#' @description This function performs a spatial join between a data frame
#'   of survey station locations and OISST sea surface temperature (SST) data.
#'   It iterates through OISST data files and finds the nearest SST observation
#'   for each station on the same date.
#'
#' @param stations An `sf` data frame of survey stations with `year`, `yrmody`,
#'   and geometry columns.
#' @param oisst_files A character vector of file paths to the OISST data frames.
#'
#' @return An `sf` data frame with the joined OISST data for all stations,
#'   or an empty tibble if no data is found.

join_oisst_to_stations <- function(stations, oisst_files) {
  # Initialize an empty list to store the results from each year.
  yearly_results <- list()

  # Loop through each OISST data file.
  for (df_path in oisst_files) {
    message(paste("Processing file:", basename(df_path)))

    # Read the OISST data for the current year.
    sstdf <- readRDS(df_path)

    # Filter stations to match the current OISST year.
    stations_yr <- stations %>%
      dplyr::filter(year == unique(sstdf$year))

    # Check if there are stations for this year before proceeding.
    if (nrow(stations_yr) > 0) {
      message("Prepping OISST data for spatial join")

      # Prepare the OISST data for the spatial join.
      sstdf_survdays <- sstdf %>%
        dplyr::mutate(
          yrmody = as.numeric(paste0(year, month, day)),
          declon = Lon,
          declat = Lat
        ) %>%
        dplyr::filter(yrmody %in% unique(stations_yr$yrmody)) %>%
        dplyr::select(-Lon, -Lat) %>%
        sf::st_as_sf(coords = c("declon", "declat"), crs = 4326, remove = FALSE)

      message("Performing spatial join for nearest-neighbor SST values")

      # Use purrr::map_dfr to perform joins and bind rows in one step.
      yr_combined <- unique(stations_yr$yrmody) %>%
        purrr::map_dfr(function(date) {
          stations_on_date <- stations_yr %>% dplyr::filter(yrmody == date)
          sst_on_date <- sstdf_survdays %>% dplyr::filter(yrmody == date)

          if (nrow(sst_on_date) > 0) {
            sf::st_join(
              stations_on_date,
              sst_on_date,
              join = nngeo::st_nn,
              k = 1,
              progress = FALSE
            )
          } else {
            return(NULL)
          }
        })

      # Add the combined results for the current year to the list.
      yearly_results[[length(yearly_results) + 1]] <- yr_combined
    }
  }

  # After the loop, combine all the yearly results into one final data frame.
  joined_data <- dplyr::bind_rows(yearly_results)

  return(joined_data)
}

#' Identify Mismatched Values Between Two Data Frames
#'
#' Compares two data frames row by row to identify columns where values do not
#' match. It returns a list containing the names of the mismatched columns
#' and the indices of rows where all values are identical. The function assumes
#' both data frames have the same dimensions and column names.
#'
#' @param data1 A data frame.
#' @param data2 A data frame to be compared against `data1`.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{`cols_mismatched`}{A character vector of unique column names
#'       that have at least one mismatched value.}
#'     \item{`rows_w_same_vals`}{An integer vector of row indices where all
#'       values in `data1` and `data2` were identical.}
#'   }
#'
#' @details
#' The function iterates through each row of the input data frames. For each
#' row, it compares the values in `data1` to the corresponding values in `data2`.
#' To ensure a proper comparison, it reorders the columns of `data2` to match
#' the column order of `data1`. The comparison is performed for non-NA values.
#' Mismatches are flagged, and the column names are collected. Additionally,
#' rows where all values match are also identified.
#'
#' @importFrom dplyr select filter
#' @importFrom tibble tibble
#'
#' @examples
#' # Create two sample data frames for comparison
#' df1 <- tibble::tibble(
#'   id = 1:5,
#'   value1 = c(10, 20, 30, 40, 50),
#'   category = c("A", "B", "C", "D", "E")
#' )
#'
#' df2 <- tibble::tibble(
#'   id = 1:5,
#'   value1 = c(10, 21, 30, 40, 50), # Mismatch in row 2
#'   category = c("A", "B", "X", "D", "E") # Mismatch in row 3
#' )
#'
#' # Use the function to find mismatches
#' results <- show_mismatched_vals(df1, df2)
#'
#' # View the output
#' print(results)
#'
#' # Expected output:
#' # $cols_mismatched
#' # [1] "value1"   "category"
#' #
#' # $rows_w_same_vals
#' # [1] 1 4 5

# a function to look at rows that are persisting after anti_join
show_mismatched_vals <- function(data1, data2) {
  output <- c()
  rows_w_same <- c()
  for (i in 1:nrow(data1)) {
    this1 <- data1[i, ] |> t()
    this2 <- data2[i, ] |>
      # reorder so columns match
      dplyr::select(colnames(this1)) |>
      # dplyr::select(-X) |>
      t()

    combined_dat <- tibble::tibble(
      names = colnames(data1),
      data1 = this1,
      data2 = this2
    ) |>
      dplyr::filter(!(is.na(data1) & is.na(data2)))

    # different_col <- which(is.na(combined_dat$data1 == combined_dat$data2))
    different_col <- which(combined_dat$data1 != combined_dat$data2)
    if (length(different_col) == 0) {
      rows_w_same <- c(rows_w_same, i)
    }

    output <- c(output, combined_dat[different_col, "names"])
  }

  out <- list(cols_mismatched = unique(output), rows_w_same_vals = rows_w_same)
  return({
    out
  })
}

#' Process and Format Forage Fish Data for State of the Ecosystem Reports
#'
#' Reads an input CSV file containing forage fish data, joins it with a
#' stratification lookup table, reshapes the data into a long format, filters it
#' by specific ecosystem production units (EPUs), and saves the processed output
#' as an RDS file.
#'
#' @param infile A character string specifying the file path to the input CSV
#'   file. This file should contain columns like `Time`, `Region`, `Estimate`,
#'   and `Std..Error.for.Estimate`.
#' @param season A character string indicating the season (e.g., "fall" or
#'   "spring"). This value is prepended to the variable names in the output.
#' @param stratlook A data frame used as a lookup table to join with the input
#'   data. It must contain columns that can be used to join with the input file
#'   (e.g., `Region` or similar).
#' @param outfile A character string specifying the file path where the
#'   processed output will be saved as an RDS file.
#'
#' @return The function invisibly returns the processed data frame, but its
#'   primary action is to save the data frame to the specified `outfile` as an
#'   RDS file.
#'
#' @details
#' The function performs several key data manipulation steps using the `dplyr`
#' and `tidyr` packages. It first joins the input data with the provided
#' `stratlook` data frame. It then selects and renames columns to a more
#' descriptive format. The `pivot_longer` function is used to convert the
#' "Forage Fish Biomass Estimate" and its standard error into a long format,
#' which is often more suitable for plotting and analysis. The data is filtered
#' to include only specific EPUs (MAB, GB, GOM, AllEPU), and a `season` prefix
#' is added to the variable names. Finally, the resulting tidy data frame is
#' saved.
#'
#' @importFrom dplyr left_join select filter mutate
#' @importFrom tidyr pivot_longer
#' @importFrom stringr str_c
#' @importFrom readr read_csv
#'
#' @examples
#' \dontrun{
#' # Assuming 'raw_data.csv' and 'strat_lookup.csv' exist with the correct columns
#'
#' # Create a sample `stratlook` data frame for the example
#' strat_lookup <- data.frame(
#'   Region = c("MAB", "GB", "GOM", "AllEPU"),
#'   Other_Info = c("A", "B", "C", "D")
#' )
#'
#' # Create a dummy input CSV file
#' write.csv(data.frame(
#'   Time = c(2020, 2021),
#'   Region = c("MAB", "GB"),
#'   Estimate = c(1.5, 2.3),
#'   Std..Error.for.Estimate = c(0.1, 0.2)
#' ), "dummy_data.csv", row.names = FALSE)
#'
#' # Run the function for the 'fall' season and save the output
#' SOEinputs(
#'   infile = "dummy_data.csv",
#'   season = "fall",
#'   stratlook = strat_lookup,
#'   outfile = "fall_forage_output.rds"
#' )
#'
#' # To verify, you can read the saved RDS file
#' processed_data <- readRDS("fall_forage_output.rds")
#' print(processed_data)
#' }
# function to wrangle the SOE indices
SOEinputs <- function(infile, season, stratlook, outfile) {
  splitoutput <- read.csv(infile)

  forageindex <- splitoutput %>%
    left_join(stratlook) %>%
    dplyr::select(
      Time,
      EPU = Region,
      "Forage Fish Biomass Estimate" = Estimate,
      "Forage Fish Biomass Estimate SE" = Std..Error.for.Estimate
    ) %>%
    tidyr::pivot_longer(
      c("Forage Fish Biomass Estimate", "Forage Fish Biomass Estimate SE"),
      names_to = "Var",
      values_to = "Value"
    ) %>%
    dplyr::filter(EPU %in% c("MAB", "GB", "GOM", "AllEPU")) %>%
    dplyr::mutate(Units = "relative grams per stomach") %>%
    dplyr::select(Time, Var, Value, EPU, Units)

  forageindex$Var <- stringr::str_c(season, forageindex$Var, sep = " ")

  #readr::write_csv(forageindex, outfile)
  saveRDS(forageindex, outfile)
}

# I can definitely help you with that.
# I'll write the documentation for your extract_cog function using the Roxygen2 format,
# which is the standard for R packages.

# Overview of the Documentation
# The documentation will be structured to be clear and helpful for anyone using the function. It will include:

# Description: A brief summary of what the function does.

# Arguments (@param): An explanation of the single input, model_fit.

# Details (@details): A more in-depth look at the internal workings of the function, particularly its reliance on the FishStatsUtils package.

# Value (@return): A description of what the function returns.

# Examples (@examples): A runnable example to show the function in action.

# Documentation for extract_cog
# Here is the Roxygen2 documentation. You can copy this code and paste it directly above your function in the R script.

# R

#'Extract Center of Gravity from a VAST Model Fit
#'
#' This function calculates the center of gravity (COG) for a spatio-temporal
#' model fit object from the VAST (Vector Autoregressive Spatio-Temporal)
#' modeling framework. It uses the `FishStatsUtils` package to perform the
#' calculation and returns the results as a data frame. It also creates and
#' then deletes a temporary directory for plotting.
#'
#' @param model_fit A list object containing a fitted VAST model. This object
#'   is typically the output of a VAST modeling run and must contain the
#'   following elements:
#'   \itemize{
#'     \item `parameter_estimates`: A list with an `SD` element.
#'     \item `Report`: The TMB report object.
#'     \item `data_list`: The TMB data list.
#'     \item `year_labels`: A vector of year labels.
#'     \item `years_to_plot`: A vector of years for plotting.
#'   }
#'
#' @return A data frame containing the calculated center of gravity. This data
#'   frame includes columns for the x and y coordinates of the COG, along with
#'   their standard errors and the year.
#'
#' @details
#' The function relies on `FishStatsUtils::plot_range_index()` to perform the
#' COG calculation. While this function also produces plots, they are not
#' returned by `extract_cog()`. Instead, a temporary directory is created for
#' these plots and then immediately deleted (`unlink`). This allows the function
#' to use the COG calculation functionality of `plot_range_index()` without
#' keeping the associated plot files.
#'
#' @importFrom here here
#' @importFrom FishStatsUtils plot_range_index
#' @importFrom utils unzip
#'
#' @examples
#' \dontrun{
#' # Assuming 'model_fit' is a pre-existing list object from a VAST run
#' # with the required components.
#'
#' # Load a sample VAST model fit object (VAST provides these for examples)
#' my_model_fit <- list(
#'   parameter_estimates = list(SD = SD),
#'   Report = Report,
#'   data_list = TmbData,
#'   year_labels = 1982:2018,
#'   years_to_plot = 1982:2018
#' )
#'
#' # Extract the center of gravity
#' cog_data <- extract_cog(my_model_fit)
#'
#' # View the first few rows of the results
#' head(cog_data)
#' }

extract_cog <- function(model_fit) {
  dir.create(here::here("temp"))
  cog <- FishStatsUtils::plot_range_index(
    Sdreport = model_fit$parameter_estimates$SD,
    Report = model_fit$Report,
    TmbData = model_fit$data_list,
    year_labels = as.numeric(model_fit$year_labels),
    years_to_plot = model_fit$years_to_plot,
    Znames = colnames(model_fit$data_list$Z_xm),
    PlotDir = here::here("temp")
  ) #already have plots, will delete this directory
  unlink(here::here("temp"), recursive = TRUE)
  return(cog)
}

#' Process and Format Forage Fish Center of Gravity Data
#'
#' Reads an RDS file containing center of gravity (COG) data, processes it by
#' renaming columns, converting it to a tidy (long) format, and adding
#' descriptive columns. The function then saves the processed data frame as a
#' new RDS file.
#'
#' @param infile A character string specifying the file path to the input RDS
#'   file. This file should contain a data frame named `COG_Table` within the
#'   list output of `FishStatsUtils::plot_range_index()`.
#' @param season A character string indicating the season (e.g., "fall" or
#'   "spring"). This value is prepended to the variable names in the output.
#' @param outfile A character string specifying the file path where the
#'   processed output will be saved as an RDS file.
#'
#' @return The function invisibly returns the processed data frame, but its
#'   primary action is to save the data frame to the specified `outfile` as an
#'   RDS file.
#'
#' @details
#' The function is designed to take the raw COG output from `FishStatsUtils::plot_range_index()`
#' and format it for reporting. It first accesses the `COG_Table` within the
#' input list. It then uses `dplyr::mutate` to create a `direction` column
#' based on the `m` variable, which distinguishes between the eastward (`m=1`) and
#' northward (`m=2`) components of the COG. The `tidyr::pivot_longer` and
#' `tidyr::unite` functions are used to reshape the data, making it tidy
#' with a single `Var` column that combines the direction and variable name.
#' Finally, it adds `Units` and `EPU` columns before saving the result.
#'
#' @importFrom dplyr mutate select
#' @importFrom tidyr pivot_longer unite
#' @importFrom stringr str_c
#' @importFrom readr read_csv
#'
#' @examples
#' \dontrun{
#' # Create a dummy COG list and save it as an RDS file
#' dummy_cog <- list(COG_Table = data.frame(
#'   Year = c(2020, 2020, 2021, 2021),
#'   m = c(1, 2, 1, 2), # 1=Eastward, 2=Northward
#'   COG_hat = c(10, 20, 15, 25),
#'   SE = c(1, 2, 1.5, 2.5)
#' ))
#' saveRDS(dummy_cog, "dummy_cog.rds")
#'
#' # Run the function for the 'fall' season and save the output
#' SOEinputsCOG(
#'   infile = "dummy_cog.rds",
#'   season = "fall",
#'   outfile = "fall_cog_output.rds"
#' )
#'
#' # To verify, you can read the saved RDS file
#' processed_data <- readRDS("fall_cog_output.rds")
#' print(processed_data)
#' }

SOEinputsCOG <- function(infile, season, outfile) {
  cogout <- readRDS(infile)

  foragecog <- as.data.frame(cogout$COG_Table) |>
    dplyr::mutate(direction = ifelse(m == 1, "Eastward", "Northward")) |>
    dplyr::select(
      "Time" = Year,
      "Forage Fish Center of Gravity" = COG_hat,
      "Forage Fish Center of Gravity SE" = SE,
      direction
    ) |>
    tidyr::pivot_longer(
      c("Forage Fish Center of Gravity", "Forage Fish Center of Gravity SE"),
      names_to = "Var",
      values_to = "Value"
    ) |>
    #direction into Var
    tidyr::unite(Var, direction:Var, sep = " ") |>
    dplyr::mutate(Units = "km", EPU = "ALLEPU") |>
    dplyr::select(Time, Var, Value, EPU, Units)

  foragecog$Var <- stringr::str_c(season, foragecog$Var, sep = " ")

  #readr::write_csv(forageindex, outfile)
  saveRDS(foragecog, outfile)
}
