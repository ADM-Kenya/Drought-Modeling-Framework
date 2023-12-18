#==============================================================================#
#' Script to resample the monthly NDVI, NDII, and LST Sentinel-3 Composites    #
#' and calculate the Baseline Mean/ STD and monthly Index Anomalies            #                                     
#'                                                                             #
#' This code is used to resample the monthly Sentinel-3 composites to match    #
#' the spatial resolution of the baseline raster files. Then, the baseline     #
#' mean and std are used to calculate the monthly index anomalies for the      #
#' three variables.                                                            #
#'                                                                             #
#' The user needs to specify:                                                  #
#' @param temp_dir The path to the temporary directory.                        #
#' @param out_dir The path to the directory to store the baseline and anomalies#
#'                files.                                                       #
#' @param baselines_dir The path to the mean and std baseline files.           #
#' @param S3_dir The path to the Sentinel-3 monthly calibrated NDVI, NDII, and # 
#'               LST composites.                                               #
#==============================================================================#



# Install and load the required packages.
packages = c("furrr", "terra", "R.utils", "Rcpp")

for(i in 1:length(packages)) {
  if(packages[i] %in% rownames(installed.packages()) == FALSE) {install.packages(packages[i])}
}

lapply(packages, library, character.only = T)



########## Paths ##########

out_dir = "E:/Maxi/01_sentinel/09_monthlyAnomalies_Sarah"

S3_dir <- "E:/Maxi/01_sentinel/08_calibratedOutput_Sarah"
baselines_dir <- "E:/Maxi/07_MODIS/05_MODIS_S3_Baselines_Anomalies_Sarah"



########## Functions ##########

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

#' Create Sequence of Dates fitted to monthly Composites
#' 
#' This function creates a monthly sequence of dates starting with the date of 
#' the first composite and ending with the date of the last. Additionally, it
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

create_Dates_Indices <- function(S3_dir){
  # List the S3 files in the directory
  S3_NDVI_files <- list.files(S3_dir, pattern = "NDVI")
  S3_NDII_files <- list.files(S3_dir, pattern = "NDII")
  
  # Empty list to store the dates
  S3_NDVI_date <- NULL
  S3_NDII_date <- NULL
  
  # Extract dates from S3 files
  for (i in 1:length(S3_NDVI_files)) {
    S3_NDVI_date[i] <- extractDate(S3_NDVI_files[i])
    S3_NDVI_date[i] <- paste0(S3_NDVI_date[i], "-01") 
  }
  
  for (i in 1:length(S3_NDII_files)) {
    S3_NDII_date[i] <- extractDate(S3_NDII_files[i])
    S3_NDII_date[i] <- paste0(S3_NDII_date[i], "-01") 
  }
  
  # Remove duplicates and convert to YYYYmm format
  S3_NDVI_date <- as.Date(S3_NDVI_date)
  S3_NDII_date <- as.Date(S3_NDII_date)
  
  S3_NDVI_monthly_date <- as.numeric(format(S3_NDVI_date, "%Y%m"))
  S3_NDII_monthly_date <- as.numeric(format(S3_NDII_date, "%Y%m"))
  
  # Remove duplicates
  monthly_NDVI_dates_un <- unique(S3_NDVI_monthly_date)
  monthly_NDII_dates_un <- unique(S3_NDII_monthly_date)
  
  # Create monthly indices
  month_NDVI_indices <- as.numeric(substr(monthly_NDVI_dates_un, 5, 6))
  month_NDII_indices <- as.numeric(substr(monthly_NDII_dates_un, 5, 6))
  
  return(list(monthly_NDVI_dates_un, monthly_NDII_dates_un, month_NDVI_indices, month_NDII_indices))
}

#' Calculate monthly Anomalies for NDVI, NDII, and LST
#' 
#' This function calculates baseline mean and std to compute the monthly 
#' anomalies of ndvi, ndii, and lst and saves the results
#' in the defined output folder.
#' @param baselines_dir The directory conatining the baseline mean and std files.
#' @param S3_dir The directory containing the Sentinel-3 monthly composites.
#' @param out_dir The directory where the output files will be saved.
#' 
#' @return None
#' 
#' @examples
#' MODIS_dir: "E:/Maxi/07_MODIS/02_monthlyComposites"
#' S3_dir: "E:/Maxi/01_sentinel/05_output"
#' out_dir: "E:/Maxi/07_MODIS/03_Baselines"

# function to calculate the baseline mean and std 
S3_calculateAnomalies <- function(baselines_dir, S3_dir, out_dir){
  
  # Import the baseline mean and std files
  mean_ndvi <- rast(list.files(baselines_dir, pattern = "Baseline_Mean_NDVI", full.names = T))
  std_ndvi <- rast(list.files(baselines_dir, pattern = "Baseline_STD_NDVI", full.names = T))
  mean_ndii <- rast(list.files(baselines_dir, pattern = "Baseline_Mean_NDII", full.names = T))
  std_ndii <- rast(list.files(baselines_dir, pattern = "Baseline_STD_NDII", full.names = T))
  mean_lst <- rast(list.files(baselines_dir, pattern = "Baseline_Mean_LST", full.names = T))
  std_lst <- rast(list.files(baselines_dir, pattern = "Baseline_STD_LST", full.names = T))
  
  ### LOAD, RESAMPLE, MERGE:
  
  # List monthly S3 composites (max NDVI, max LST, median NDII)
  S3_files <- list.files(S3_dir, full.names = TRUE)
  
  for (i in S3_files) {
    S3_ndvi_files <- list.files(S3_dir, pattern = "NDVI", full.names = TRUE)
    S3_ndii_files <- list.files(S3_dir, pattern = "NDII", full.names = TRUE)
    S3_lst_files <- list.files(S3_dir, pattern = "LST", full.names = TRUE)
  }
  
  # Create three empty lists to fill with the resampled raster files
  ndvi_stack <- list()
  ndii_stack <- list()
  lst_stack <- list()
  
  # Resample S3 monthly composites to the extent of the baseline files
  for (i in 1:length(S3_ndvi_files)) {
    ndvi_raster <- rast(S3_ndvi_files[i])
    ndvi_resampled <- terra::resample(ndvi_raster, mean_ndvi, method = "bilinear")
    ndvi_stack <- append(ndvi_stack, ndvi_resampled)
  }
  
  for (i in 1:length(S3_ndii_files)) {
    ndii_raster <- rast(S3_ndii_files[i])
    ndii_resampled <- terra::resample(ndii_raster, mean_ndii, method = "bilinear")
    ndii_stack <- append(ndii_stack, ndii_resampled)
  }
  
  for (i in 1:length(S3_lst_files)) {
    lst_raster <- rast(S3_lst_files[i])
    lst_resampled <- terra::resample(lst_raster, mean_lst, method = "bilinear")
    lst_stack <- append(lst_stack, lst_resampled)
  }
  
  print("Resampling done!")
  
  ### EXTRACT DATES:
  
  # Create sequence of dates fitted to file dates
  DateIndices_list <- create_Dates_Indices(S3_dir)
  
  monthly_NDVI_dates_un <- unlist(DateIndices_list[1])
  monthly_NDII_dates_un <- unlist(DateIndices_list[2])
  
  month_NDVI_indices <- unlist(DateIndices_list[3])
  month_NDII_indices <- unlist(DateIndices_list[4])

  ### CALCULATE ANOMALIES:

  ndvi_anom <- list()
  ndii_anom <- list()
  lst_anom <- list()

  # Calculate anomalies for the three different variables
  for (i in 1:nlyr(ndvi_stack)) {
    ndvi_anom[[i]] <- (ndvi_stack[[i]] - mean_ndvi[[month_NDVI_indices[i]]]) /
      std_ndvi[[month_NDVI_indices[i]]]
  }

  for (i in 1:nlyr(ndii_stack)) {
    ndii_anom[[i]] <- (ndii_stack[[i]] - mean_ndii[[month_NDII_indices[i]]]) /
      std_ndii[[month_NDII_indices[i]]]
  }

  for (i in 1:nlyr(lst_stack)) {
    lst_anom[[i]] <- (lst_stack[[i]] - mean_lst[[month_NDVI_indices[i]]]) /
      std_lst[[month_NDVI_indices[i]]]
  }

  # Convert lists to raster files
  ndvi_anom <- rast(ndvi_anom)
  ndii_anom <- rast(ndii_anom)
  lst_anom <- rast(lst_anom)
  
  # Change names of raster files
  ndvi_anom_names <- NULL
  ndii_anom_names <- NULL
  lst_anom_names <- NULL
  
  for(i in 1:nlyr(ndvi_anom)){
    ndvi_anom_names[i] <- paste0("NDVI Anomalies_Month=", monthly_NDVI_dates_un[i])
  } 
  for(i in 1:nlyr(ndii_anom)){
    ndii_anom_names[i] <- paste0("NDII Anomalies_Month=", monthly_NDII_dates_un[i])
  } 
  for(i in 1:nlyr(lst_anom)){
    lst_anom_names[i] <- paste0("LST Anomalies_Month=", monthly_NDVI_dates_un[i])
  } 
  
  names(ndvi_anom) <- ndvi_anom_names
  names(ndii_anom) <- ndii_anom_names
  names(lst_anom) <- lst_anom_names

  # Set z-dimension
  time(ndvi_anom) <- monthly_NDVI_dates_un
  time(ndii_anom) <- monthly_NDII_dates_un
  time(lst_anom) <- monthly_NDVI_dates_un

  print("Anomalies computed!")
  
  ### COMBINE OLD AND NEW ANOMALIES
  ndvi_anom_old <- list.files(out_dir, pattern = "NDVI", full.names = T)
  ndii_anom_old <- list.files(out_dir, pattern = "NDII", full.names = T)
  lst_anom_old <- list.files(out_dir, pattern = "LST", full.names = T)
  
  # Combine with newly computed anomalies and check if there are duplicates
  if(length(ndvi_anom_old) > 0){
    ndvi_anom_old <- rast(ndvi_anom_old)
    ndvi_unique_names <- unique(c(names(ndvi_anom_old), ndvi_anom_names))
    for (i in ndvi_unique_names) { if(!i %in% names(ndvi_anom_old)){
      ndvi_anom <- ndvi_anom_old
      ndvi_anom[[i]] <- ndvi_anom_new[[i]]}
    } 
  }
  
  if(length(ndii_anom_old) > 0){
    ndii_anom_old <- rast(ndii_anom_old)
    ndii_unique_names <- unique(c(names(ndii_anom_old), ndii_anom_names))
    for (i in ndii_unique_names) { if(!i %in% names(ndii_anom_old)){
      ndii_anom <- ndii_anom_old
      ndii_anom[[i]] <- ndii_anom[[i]]}
    }
  }
  
  if(length(lst_anom_old) > 0){
    lst_anom_old <- rast(lst_anom_old)
    lst_unique_names <- unique(c(names(lst_anom_old), lst_anom_names))
    for (i in lst_unique_names) { if(!i %in% names(lst_anom_old)){
      lst_anom <- lst_anom_old
      lst_anom[[i]] <- lst_anom[[i]]}
    }    
  }
  
  print("Old and new Anomalies combined!")
  
  ### SAVE RESULTS
  raster_data <- list(ndvi_anom, ndii_anom, lst_anom)
  
  # Set output file path and long variable names
  names <-  c("NDVI_Anomalies", "NDII_Anomalies", "LST_Anomalies")

  var_names <- c("NDVI Anomalies", "NDII Anomalies", "LST Anomalies")

  long_names <- c("Normalized Difference Vegetation Index (NDVI): monthly Index Anomalies",
                  "Normalized Difference Infrared Index (NDII): monthly Index Anomalies",
                  "Land Surface Temperature (LST): monthly Anomalies")

  # Write results to netCDF files
  for (i in 1:length(raster_data)) {
    writeCDF(raster_data[[i]], filename = paste(out_dir, "/", names[i], ".nc", sep = ""),
             varname = var_names[i], unit = "None", longname = long_names[i], zname = "Month", overwrite = T)
  }
  print("Anomalies written!")
}



########## Apply Function ##########

S3_calculateAnomalies(baselines_dir, S3_dir, out_dir)
