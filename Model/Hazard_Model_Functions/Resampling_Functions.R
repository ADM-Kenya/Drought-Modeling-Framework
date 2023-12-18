#==============================================================================#
#' Script for Re-Sampling the  Model Data                                      #
#'                                                                             #
#' This code is used for re-sampling the NDVI/ NDII/ and LST Anomalies as well #
#' the Tamsat SPI3 Precipitation and Copernicus Land Cover Data to the same    #
#' spatial resolution (~ 0.01 degree of latitude/ ~ 1.113 km).                 #
#'                                                                             #
#' The user needs to specify:                                                  #
#' @param anomalies_dir The path to directory containing the NDVI, NDII, and   #
#'                      LST Anomalies.                                         #
#' @param spi3_dir The path to the directory containing the Tamsat SPI3        #
#'                 Precipitation data.                                         #
#' @param lc_file The Copernicus Land Cover file (land cover mask).            #
#' @param output_dir The directory to store the resampled and masked anomalies #
#'                   and SPI3 data.                                            #
#==============================================================================#



# Install and load required packages

install.packages("terra")
library(terra)



########## Paths and Variables ##########

### PATHS:
anomalies_dir <- "E:/Maxi/07_MODIS/05_MODIS_S3_Baselines_Anomalies_Sarah/"
spi3_dir <- "E:/Maxi/05_tamsat/spi_output/"
lc_file <- "E:/Maxi/04_copernicus_landcover/copernicusLC_cropHerbShrubMask_100m_Kenya.tif"

output_dir <- "E:/Maxi/06_Drought_Model/02_inputData_Sarah/05_inputData_copernicusLC_cropHerbShrub/"



########## Functions ##########

#' Extract Date from Filename
#' 
#' This function extracts the date from the filename.
#' 
#' @param filename The filename to extract the date from.

extractDate <- function(filename){
  # regular expression pattern to match the date in the file name
  date_pattern <- "\\d{4}-\\d{2}"
  
  # extract the date
  date <- regmatches(filename, regexpr(date_pattern, filename))
  
  # check if the date was found in the filename
  if (length(date) == 1) {
    return(date)
  } else {
    return('Date not found')
  }
}


#' Resample and Mask Data
#' 
#' This function resamples the model input data (NDVI/ NDII/ LST Anomalies and 
#' SPI3 Precipitation) to the same spatial resolution (~ 0.01 degree of latitude)
#' and masks the resampled data with a Copernicus land cover map.
#' 

resample_LCmask <- function(anomalies_dir, spi3_dir, lc_file, output_dir){
  
  ########## Import Data ##########
  
  ### ANOMALIES:
  ndviAnomalies_files <- list.files(path = anomalies_dir, pattern = "NDVI_Anomalies.nc$", 
                                    recursive = T, full.names = T)
  ndiiAnomalies_files <- list.files(path = anomalies_dir, pattern = "NDII_Anomalies.nc$", 
                                    recursive = T, full.names = T)
  lstAnomalies_files <- list.files(path = anomalies_dir, pattern = "LST_Anomalies.nc$", 
                                   recursive = T, full.names = T)
  
  ndviAnomalies_stack <- rast(ndviAnomalies_files)
  ndiiAnomalies_stack <- rast(ndiiAnomalies_files)
  lstAnomalies_stack <- rast(lstAnomalies_files)
  
  ### SPI3:
  spi3_files <- list.files(path = spi3_dir, pattern = "\\d{4}-\\d{2}\\_SPI_3.tif$", full.names = T)
  
  # Remove files before January 2001
  spi3_files <- spi3_files[-c(1:216)]
  
  spi3_stack <- rast(spi3_files)
  
  ### COPERNICUS LAND COVER:
  lc_rast <- rast(lc_file)
  
  ########## Extract Date ##########
  
  # The date needs to be extracted to set the z-dimension later
  date <- NULL
  
  # Extract the date from the SPI3 precipitation files
  for (i in 1:length(spi3_files)) {
    date[i] <- extractDate(spi3_files[i])
    date[i] <- paste0(date[i], "-01")
  }
  
  date <- as.Date(date)
  
  monthly_dates <- as.numeric(format(date, "%Y%m"))
  
  monthly_dates_un <- unique(monthly_dates)
  
  ########## Resampling ##########
  
  ### SPI:
  
  # Define an aggregation factor to resample TAMSAT's spatial resolution from 
  # degree to square km (original resolution 0.0375, approx. 4km)
  aggr_fact <- 3.75
  
  spi3_stack_res <- terra::disagg(spi3_stack, fact = aggr_fact)
  
  ### ANOMALIES:
  
  ndviAnomalies_stack_res <- terra::resample(ndviAnomalies_stack, spi3_stack_res, method = "bilinear")
  ndiiAnomalies_stack_res <- terra::resample(ndiiAnomalies_stack, spi3_stack_res, method = "bilinear")
  lstAnomalies_stack_res <- terra::resample(lstAnomalies_stack, spi3_stack_res, method = "bilinear")
  
  ### LANDCOVER:
  
  lc_res <- terra::resample(lc_rast, spi3_stack_res, method = "near")
  lc_res[lc_res != 1] <- NA
  
  ########## Land Cover Masking ##########
  
  # Mask anomalies with Copernicus land cover mask
  ndviAnomalies_stack_res_lc <- terra::mask(ndviAnomalies_stack_res, lc_res)
  ndiiAnomalies_stack_res_lc <- terra::mask(ndiiAnomalies_stack_res, lc_res)
  lstAnomalies_stack_res_lc <- terra::mask(lstAnomalies_stack_res, lc_res)
  
  # Mask SPI3 precipitation data with Copernicus land cover mask
  spi3_stack_res_lc <- terra::mask(spi3_stack_res, lc_res)
  
  ########## Write Data ##########
  
  ### LAND COVER MASK:
  lc_res_path <- sub(".tif", "_res.tif", lc_file)
  
  #terra::writeRaster(lc_res, filename = lc_res_path,
  #                   filetype = "GTiff", overwrite = TRUE)
  
  ### ANOMALIES:
  
  # Set z-dimension
  time(ndviAnomalies_stack_res_lc) <- monthly_dates_un
  time(ndiiAnomalies_stack_res_lc) <- monthly_dates_un
  time(lstAnomalies_stack_res_lc) <- monthly_dates_un
  
  raster_data <- list(ndviAnomalies_stack_res_lc, ndiiAnomalies_stack_res_lc, lstAnomalies_stack_res_lc)
  
  # Create metadata (names for netCDF file)
  names <- c("NDVI_Anomalies_res_lc", "NDII_Anomalies_res_lc", "LST_Anomalies_res_lc")
  
  var_names <- c("NDVI Anomalies", "NDII Anomalies", "LST Anomalies")
  
  long_names <- c("Normalized Difference Vegetation Index (NDVI): Monthly Index Anomalies",
                  "Normalized Difference Infrared Index (NDII): Monthly Index Anomalies",
                  "Land Surface Temperature (LST): Monthly Anomalies")
  
  # Write and save files in output folder (one netCDF file for every index)
  for(i in 1:length(raster_data)) {
    writeCDF(raster_data[[i]], filename = paste(output_dir, names[i], ".nc", sep = ""), 
             varname = var_names[i], unit = "None", longname = long_names[i], zname = "Month")
  }
  
  ### SPI:
  
  # Set z-dimension
  time(spi3_stack_res_lc) <- monthly_dates_un
  
  # Create metadata (names for netCDF file)
  name <- "SPI3_res_lc"
  var_name <- "SPI3"
  long_name <- "3 monthly Standardized Precipitation Index (SPI3)"
  
  # Write and save file in output folder
  writeCDF(spi3_stack_res_lc, filename = paste(output_dir, name, ".nc", sep = ""), 
           varname = var_name, unit = "None", longname = long_name, zname = "Month")
  
}
resample_LCmask(anomalies_dir, spi3_dir, lc_file, output_dir)
