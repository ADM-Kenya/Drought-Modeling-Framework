#==============================================================================#
#' Script to re-sample the monthly NDVI, NDII, and LST Sentinel-3 composites,  #
#' combine them with the MODIS composites and calculate the Baseline Mean/ STD #
#' and monthly Index Anomalies                                                 #
#'                                                                             #
#' This code is used to re-sample the monthly Sentinel-3 Index composites to   #
#' match the spatial resolution of the  MODIS composites. Then, the S3 and     #
#' MODIS composites are merged before the baseline mean and std are calculated #
#' over the entire date range of MODIS and S3 files. The baseline mean and std #
#' are then used to calculate the monthly index anomalies for the three        # 
#' variables.                                                                  #
#'                                                                             #
#' The user needs to specify:                                                  #
#' @param temp_dir The path to the temporary directory.                        #
#' @param anom_out_dir The path to the directory to store the anomalies files. #
#' @param baselines_out_dir The path to the directory to store the baseline    #
#'                          files.
#' @param S3_dir The path to the Sentinel-3 monthly calibrated NDVI, NDII, and # 
#'               LST composites.                                               #
#' @param MODIS_dir The path to the MODIS monthly NDVI, NDII, and LST          #
#'                  composites.                                                #
#==============================================================================#


########## Functions ##########

# ---------------------------------------------------------------------------- #
#' Extract Date from Sentinel-3 monthly Composites
#' 
#' This function extracts the dates from a list of Sentinel-3 monthly composites.
#' 
#' @param filename The file name of a Sentinel-3 monthly composite
#' 
#' @return None
#' 
#' @examples
#' filename: "NDVI_20190303_0179_022_320_3060"
#' return: "2021-09"
# ---------------------------------------------------------------------------- #

extractDate <- function(filename){
  # Regular expression pattern to match the date in the file name
  date_pattern <- "\\d{4}-\\d{2}"
  
  # Extract the date
  date <- regmatches(filename, regexpr(date_pattern, filename))
  
  # Check if the date was found in the file name
  if (length(date) == 1) {
    return(date)
  } else {
    return('Date not found')
  }
}

# ---------------------------------------------------------------------------- #
#' Create Sequence of Dates fitted to monthly Composites
#' 
#' This function creates a monthly sequence of dates starting with the date of 
#' the first composite and ending with the date of the last.Additionally, it
#' create a sequences of indices for the months.
#' 
#' @param S3_dir The directory containing the Sentinel-3 monthly composites.
#' 
#' @return A list containing the monthly dates and the indices.
#' 
#' @examples
#' S3_dir: "E:/Maxi/01_sentinel/05_output"
#' return:
#'   monthly_dates_un: c("202109", "202110", ...)
#'   month_indices: c(01, 02, 03, ...)
# ---------------------------------------------------------------------------- #

create_Dates_Indices <- function(S3_dir){
  # Create sequence of dates for the MODIS files
  start_date <- as.Date("2001-01-01")
  end_date <- as.Date("2020-12-01")
  
  MODIS_date <- seq.Date(start_date, end_date, by = "1 month")
  
  # List the S3 files in the directory
  S3_files <- list.files(S3_dir)
  
  # Empty list to store the dates
  S3_date <- NULL
  
  # Extract dates from S3 files
  for (i in 1:length(S3_files)) {
    S3_date[i] <- extractDate(S3_files[i])
    S3_date[i] <- paste0(S3_date[i], "-01") 
  }
  
  # Remove duplicates and convert to YYYYmm format
  S3_date <- unique(S3_date)
  
  MODIS_date <- as.Date(MODIS_date)
  S3_date <- as.Date(S3_date)
  
  MODIS_monthly_date <- as.numeric(format(MODIS_date, "%Y%m"))
  S3_monthly_date <- as.numeric(format(S3_date, "%Y%m"))
  
  # Combine the MODIS and S3 dates
  monthly_dates_combined <- c(MODIS_monthly_date, S3_monthly_date)
  
  # Remove duplicates
  monthly_dates_un <- unique(monthly_dates_combined)
  
  # Create monthly indices
  month_indices <- as.numeric(substr(monthly_dates_un, 5, 6))
  
 return(list(monthly_dates_un, month_indices))
}

# ---------------------------------------------------------------------------- #
#' Create SpatRaster stack for S3 composite files
#' 
#' @param S3_file_list A list containing the paths to the S3 monthly composites.
#' @param modis_ref SpatRaster object of the MODIS monthly composite 
#'                  used as reference for resampling

#' @return  A SpatRaster stack of all S3 monthly composites 
#'          defined in the 'S3_file_list' object
# ---------------------------------------------------------------------------- #

create_resampled_stack <- function(S3_file_list, modis_ref, index_name) {
  # Initialize an empty list to store resampled rasters
  resampled_rasters <- list()
  layer_names <- vector("character", length(S3_file_list))  # Prepare a character vector to hold layer names
  for (i in seq_along(S3_file_list)) {
    file <- S3_file_list[i]
    s3_raster <- rast(file)
    resampled_raster <- terra::resample(s3_raster, modis_ref, method="bilinear")
    # print(resampled_raster)
    
    # Extract year and month from the file path
    date_part <- sub(".*_(\\d{4}-\\d{2})_.*", "\\1", basename(file))
    layer_names[i] <- paste(index_name, date_part, sep="_")  # Assign the extracted date part to the layer names vector
    
    resampled_rasters <- c(resampled_rasters, list(resampled_raster))
    # print(resampled_rasters)
  }
  raster_stack <- rast(resampled_rasters)
  names(raster_stack) <- layer_names  # Set names of layers in the SpatRasterStack
  return(raster_stack)
}

# ---------------------------------------------------------------------------- #
#' Create SpatRaster stack for S3 composite files
#' 
#' Function to rename layers based on a continuous 
#' monthly sequence starting from a specific date
#' 
#' @param raster SpatRaster object
#' @param varname name for the variable
#' @param start_date start_date of the time-series (first layer)
#' @param end_date end date of the time-series (last layer)
#' 
#' @return A renamed SpatRaster object
# ---------------------------------------------------------------------------- #

rename_MODIS_layers_date_month <- function(raster, varname, start_date, end_date) {
  date_seq <- seq(as.Date(start_date), as.Date(end_date), by="month")
  year_month_labels <- format(date_seq, "%Y-%m")
  new_names <- paste(varname, year_month_labels, sep="_")
  names(raster) <- new_names
  return(raster)
}

# ---------------------------------------------------------------------------- #
#' Resample Sentinel-3 Composites and merge with MODIS Composites
#' 
#' This function resamples the Sentinel-3 monthly composites of every index to
#' the same spatial resolution of the MODIS monthly composites of the 
#' corresponding index. Subsequently, the S3 and MODIS composites of the same
#' index are merged to a combined raster stack.
#' 
#' @param S3_dir The directory containing the Sentinel-3 monthly composites.
#' @param MODIS_dir The directory containing the MODIS monthly composites.
#' 
#' @return A list containing three combined raster stacks. One for every variable.
#' 
#' @examples
#' S3_dir: "E:/Maxi/01_sentinel/05_output"
#' MODIS_dir: "E:/Maxi/07_MODIS/02_monthlyComposites"
#' return: ndvi_combined, ndii_combined, lst_combined
# ---------------------------------------------------------------------------- #

resampling_merging_S3_MODIS <- function(S3_dir, MODIS_dir){
  
  # Load MODIS monthly composites (max NDVI, max LST, median NDII)
  MODIS_ndvi <- rast(list.files(MODIS_dir, "NDVI", full.names = TRUE))
  MODIS_ndii <- rast(list.files(MODIS_dir, "NDII", full.names = TRUE))
  MODIS_lst <- rast(list.files(MODIS_dir, "LST", full.names = TRUE))
  
  # Rename MODIS layers
  MODIS_ndvi <- rename_MODIS_layers_date_month(MODIS_ndvi, "NDVI", "2001-01-01", "2020-12-01")
  MODIS_ndii <- rename_MODIS_layers_date_month(MODIS_ndii, "NDII", "2001-01-01", "2020-12-01")
  MODIS_lst <- rename_MODIS_layers_date_month(MODIS_lst, "LST", "2001-01-01", "2020-12-01")
  
  # List monthly S3 composites (max NDVI, max LST, median NDII)
  S3_ndvi_files <- list.files(S3_dir, pattern = "NDVI", full.names = TRUE)
  S3_ndii_files <- list.files(S3_dir, pattern = "NDII", full.names = TRUE)
  S3_lst_files <- list.files(S3_dir, pattern = "LST", full.names = TRUE)

  ### RESAMPLE:
  
  # Create resampled stacks for NDVI, NDII, and LST
  S3_ndvi_stack <- create_resampled_stack(S3_ndvi_files, MODIS_ndvi, "NDVI")
  S3_ndii_stack <- create_resampled_stack(S3_ndii_files, MODIS_ndii, "NDII")
  S3_lst_stack <- create_resampled_stack(S3_lst_files, MODIS_lst, "LST")

  ### MERGE:
  
  # Logical conditions for filtering based on dates
  merge_layers <- function(modis_stack, s3_stack) {
    modis_dates <- gsub("^[^_]*_", "", names(MODIS_ndvi)) # Remove prefix
    s3_dates <- gsub("^[^_]*_", "", names(S3_ndvi_stack))      # Remove prefix
    
    # print(s3_dates)
    
    # Identify layers by year
    modis_pre_2020 <- modis_stack[[c(!grepl("^2020", modis_dates))]]
    modis_2020_03_04 <- modis_stack[[c(grepl("2020-(03|04)", modis_dates))]]
    
    s3_2020_01_02 <- s3_stack[[c(grepl("2020-(01|02)", s3_dates))]]
    s3_post_2020_03_04 <- s3_stack[[c(!grepl("2020-(01|02)", s3_dates))]]

    # Combine the stacks
    combined <- c(modis_pre_2020, s3_2020_01_02, modis_2020_03_04, s3_post_2020_03_04)
    return(combined)
  }
  
  # Combine the raster stacks
  ndvi_combined <- merge_layers(MODIS_ndvi, S3_ndvi_stack)
  ndii_combined <- merge_layers(MODIS_ndii, S3_ndii_stack)
  lst_combined <- merge_layers(MODIS_lst, S3_lst_stack)

  return(list(ndvi_combined, ndii_combined, lst_combined))
}

# ---------------------------------------------------------------------------- #
#' Calculate Baseline Mean/ STD and Anomalies for NDVI, NDII, and LST
#' 
#' This function calculates the baseline mean and std for the combined MODIS and 
#' S3 raster stacks of NDVI, NDII, and LST. Then it uses these baselines to 
#' compute the monthly anomalies of the three variables and saves the results
#' in the defined output folder.
#' 
#' @param MODIS_dir The directory containing the MODIS monthly composites.
#' @param S3_dir The directory containing the Sentinel-3 monthly composites.
#' @param anom_out_dir The directory where the anomalies output files will be saved.
#' @param baselines_out_dir The directory where the anomalies output files will 
#'                          be saved.
#' 
#' @return None
#' 
#' @examples
#' MODIS_dir: "E:/Maxi/07_MODIS/02_monthlyComposites"
#' S3_dir: "E:/Maxi/01_sentinel/05_output"
#' out_dir: "E:/Maxi/07_MODIS/03_Baselines"
# ---------------------------------------------------------------------------- #

# function to calculate the baseline mean and std
MODIS_S3_calculateBaselinesAnomalies <- function(MODIS_dir, S3_dir, anom_out_dir, baselines_out_dir){

  ### LOAD, RESAMPLE, MERGE:
  MODIS_dir <- MODIS_monthlyComposites_dir
  S3_dir <- S3_calibratedMonthlyComposites_dir
  anom_out_dir <- MODIS_S3_monthlyAnomalies_dir
  baselines_out_dir <- MODIS_S3_Baselines_dir
  
  # Load S3 data, resample files to the same spatial resolution as the MODIS files and merge both
  ResampleMerge <- resampling_merging_S3_MODIS(S3_dir, MODIS_dir)
 
  #Read data
  ndvi_combined <- rast(unlist(ResampleMerge[1]))
  ndii_combined <- rast(unlist(ResampleMerge[2]))
  lst_combined <- rast(unlist(ResampleMerge[3]))
  
  print("Resampling and merging done!")

  # write results in a temporary directory
  temp_dir <- tempdir()
  
  writeRaster(ndvi_combined, paste0(temp_dir, "/ndvi_combined.tif"), overwrite = T)
  writeRaster(ndii_combined, paste0(temp_dir, "/ndii_combined.tif"), overwrite = T)
  writeRaster(lst_combined, paste0(temp_dir, "/lst_combined.tif"), overwrite = T)

  # Remove objects in workspace to reduced the used memory
  rm(ResampleMerge)
  rm(ndvi_combined)
  rm(ndii_combined)
  rm(lst_combined)
  gc()
  
  # Read the data back in from temporary directory
  ndvi_combined <- rast(paste0(temp_dir, "/ndvi_combined.tif"))
  ndii_combined <- rast(paste0(temp_dir, "/ndii_combined.tif"))
  lst_combined <- rast(paste0(temp_dir, "/lst_combined.tif"))
  
  print(ndvi_combined)
  print(ndii_combined)
  print(lst_combined)
  
  print("Saved resampled + merged MODIS_S3 stacks in 'rtemp' folder done!")

  ### EXTRACT DATES:

  # Create sequence of dates fitted to file dates
  DateIndices_list <- create_Dates_Indices(S3_dir)
  monthly_dates_un <- unlist(DateIndices_list[1])
  month_indices <- unlist(DateIndices_list[2])

  ### CALCULATE BASELINES:

  # Set z-dimension
  time(ndvi_combined) <- month_indices
  time(ndii_combined) <- month_indices
  time(lst_combined) <- month_indices

  # Calculate baseline mean and std for the different variables
  mean_ndvi <- tapp(ndvi_combined, month_indices, mean, na.rm = T)
  std_ndvi <- tapp(ndvi_combined, month_indices, sd, na.rm = T)

  mean_ndii <- tapp(ndii_combined, month_indices, mean, na.rm = T)
  std_ndii <- tapp(ndii_combined, month_indices, sd, na.rm = T)

  mean_lst <- tapp(lst_combined, month_indices, mean, na.rm = T)
  std_lst <- tapp(lst_combined, month_indices, sd, na.rm = T)

  # Set z-dimension
  time(mean_ndvi) <- unique(month_indices)
  time(std_ndvi) <- unique(month_indices)

  time(mean_ndii) <- unique(month_indices)
  time(std_ndii) <- unique(month_indices)

  time(mean_lst) <- unique(month_indices)
  time(std_lst) <- unique(month_indices)

  print("Baselines computed!")

  ### CALCULATE ANOMALIES:

  ndvi_anom <- list()
  ndii_anom <- list()
  lst_anom <- list()

  # Calculate anomalies for the different variables
  for (i in 1:nlyr(ndvi_combined)) {
    ndvi_anom[[i]] <- (ndvi_combined[[i]] - mean_ndvi[[paste("X", month_indices[i], sep = "")]]) /
      std_ndvi[[paste("X", month_indices[i], sep = "")]]
  }

  for (i in 1:nlyr(ndii_combined)) {
    ndii_anom[[i]] <- (ndii_combined[[i]] - mean_ndii[[paste("X", month_indices[i], sep = "")]]) /
      std_ndii[[paste("X", month_indices[i], sep = "")]]
  }

  for (i in 1:nlyr(lst_combined)) {
    lst_anom[[i]] <- (lst_combined[[i]] - mean_lst[[paste("X", month_indices[i], sep = "")]]) /
      std_lst[[paste("X", month_indices[i], sep = "")]]
  }

  # Convert lists to raster files
  ndvi_anom <- rast(ndvi_anom)
  ndii_anom <- rast(ndii_anom)
  lst_anom <- rast(lst_anom)

  # Set z-dimension
  time(ndvi_anom) <- monthly_dates_un
  time(ndii_anom) <- monthly_dates_un
  time(lst_anom) <- monthly_dates_un
  
  print(object.size(ndvi_anom))
  print(object.size(ndii_anom))
  print(object.size(lst_anom))
  
  print(lst_anom)
  
  print("Anomalies computed!")
  
  ### SAVE RESULTS

  # List with all data to be written
  anom_data <-  list(ndvi_anom, ndii_anom, lst_anom)
  baselines_data <- list(mean_ndvi, std_ndvi, mean_ndii, std_ndii, mean_lst, std_lst)

  # Set output file path and long variable names
  names_anom <-  c("NDVI_Anomalies", "NDII_Anomalies", "LST_Anomalies")
  names_baselines <- c("Baseline_Mean_NDVI", "Baseline_STD_NDVI", "Baseline_Mean_NDII",
                       "Baseline_STD_NDII", "Baseline_Mean_LST", "Baseline_STD_LST")

  varnames_anom <- c("NDVI Anomalies", "NDII Anomalies", "LST Anomalies")
  varnames_baselines <- c("Baseline Mean NDVI", "Baseline STD NDVI", "Baseline Mean NDII",
                          "Baseline STD NDII", "Baseline Mean LST", "Baseline STD LST")

  longnames_anom <- c("Normalized Difference Vegetation Index (NDVI): monthly Index Anomalies",
                  "Normalized Difference Infrared Index (NDII): monthly Index Anomalies",
                  "Land Surface Temperature (LST): monthly Anomalies")
  longnames_baselines <- c("Baseline Mean Normalized Difference Vegetation Index (NDVI)",
                           "Baseline STD Normalized Difference Vegetation Index (NDVI)",
                           "Baseline Mean Normalized Difference Infrared Index (NDII)",
                           "Baseline STD Normalized Difference Infrared Index (NDII)",
                           "Baseline Mean Land Surface Temperature (LST)",
                           "Baseline STD Land Surface Temperature (LST)")

  # Write results to netCDF files
  for (i in 1:length(anom_data)) {
    writeCDF(anom_data[[i]], filename = paste(anom_out_dir, "/", names_anom[i], ".nc", sep = ""),
             varname = varnames_anom[i], unit = "None", longname = longnames_anom[i], zname = "Month")
  }
  
  for (i in 1:length(baselines_data)) {
    writeCDF(baselines_data[[i]], filename = paste(baselines_out_dir, "/", names_baselines[i], ".nc", sep = ""),
             varname = varnames_baselines[i], unit = "None", longname = longnames_baselines[i], zname = "Month")
  }
  print("Baselines and Index Anomalies written!")
}