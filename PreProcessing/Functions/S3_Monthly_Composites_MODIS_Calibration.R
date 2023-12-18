#==============================================================================# 
#' Script to create and calibrate monthly Sentinel-3 NDVI, NDII, and LST       #
#' composites.                                                                 #
#'                                                                             #
#' The code is used to create monthly composites of the Land Surface           #
#' Temperature (LST), the Normalized Difference Vegetation Index (NDVI) and the#
#' Normalized Difference Infrared Index (NDII) from the raw daily Sentinel-3   #
#' files. The monthly NDVI and LST composites are computed based on the monthly#
#' maxima and the NDII based on the median value. Additionally, the resulting  #
#' monthly composites are calibrated to match with the composites from the     #
#' MODIS Sensor.                                                               #
#'                                                                             #
#' The user needs to specify:                                                  #
#'                                                                             #
#' @param gdal_path The path to the GDAL installation directory.               #
#' @param lst_dir The path to the raw S3 LST files.                            #
#' @param ndvi_dir The directory containing the NDVI files.                    #
#' @param ndii_dir The directory containing the NDII files.                    #
#' @param output_dir The output directory to store the monthly composites.     #
#' @param calibratedOutput_dir The path to the directory to store the          #
#'                             calibrated monthly composites.                  #
#' @param date_range The range of dates for which the monthly composites are   #
#'                   created.                                                  #
#==============================================================================#  



# Install and load all required packages.
install.packages("terra")
install.packages("tools")

library(terra)
library(tools)


########## Paths and Variables ##########

### PATHS:
gdal_path <- "C:\\Program Files\\QGIS 3.22.5\\bin"
ndvi_dir <- "E:/Maxi/Rwanda_S3_LST/S3_LST_2020-2023_preprocessed/NDVI"
ndii_dir <- "E:/Maxi/01_sentinel/04_processed_download/NDII"
lst_dir <- "E:/Maxi/Rwanda_S3_LST/S3_LST_2020-2023_preprocessed/LST"
output_dir <- "E:/Maxi/Rwanda_S3_LST/S3_LST_2020-2023_monthlyComposites"
calibratedOutput_dir <- "E:/Maxi/01_sentinel/08_calibratedOutput_Sarah_TEST/"


### VARIABLES:
date_range <- c("2020-01-01", "2023-10-31")



########## Functions ##########

#' Extract date from file path.
#'
#' This nested function extracts the date portion from a file path.
#'
#' @param filepath A character string representing the full file path.
#' @return The date portion extracted from the file name, in the Date format.

extract_date <- function(filepath) {
  # Extract the basename from the full file path
  basename <- basename(filepath)     
  
  # Extract the date portion
  date <- substr(basename, regexpr("_", basename) + 1, regexpr("_", basename) + 10)        
  
  # Convert date to Date format
  as.Date(date, format = "%Y%m%d")                                              
}


#' Filter strings by date range.
#'
#' This function filters a vector of strings based on a specified date range.
#'
#' @param strings A vector of strings to filter.
#' @param date_range A character vector containing the start and end dates of 
#'                   the desired range in the format "YYYY-MM-DD".
#' @return A vector of strings that fall within the specified date range.
#' 
#' @examples
#' strings <- c("E:/Thalis/04_processed_download/LST/LST_20220801_0179_069_006_2880.tif", 
#'              "E:/Thalis/04_processed_download/LST/LST_20220802_0179_069_020_3060.tif")
#' date_range <- c("2022-02-01", "2023-04-30")
#' filter_by_date_range(strings, date_range)
#' @export

filter_by_date_range <- function(strings, date_range) {
  start_date <- as.Date(date_range[1])
  end_date <- as.Date(date_range[2])
  
  # Extract dates from file strings
  file_dates <- unname(sapply(strings, extract_date))                           
  filtered_strings <- strings[file_dates >= start_date & file_dates <= end_date]
  return(filtered_strings)
}


#' Extract year and month from file path.
#'
#' This internal function extracts the year and month from a given file path. 
#' It assumes that the basename of the file path follows the format "LST_YYYYMMDD_..."
#' where "YYYY" represents the year, "MM" represents the month, 
#' and the rest of the characters are irrelevant for this extraction.
#'
#' @param filepath A character string representing a file path
#' @return A character string representing the year and month in the format "YYYY-MM"

extract_year_month <- function(filepath) {
  # Extract the basename from the full file path
  basename <- basename(filepath)  
  
  # Extract the year portion
  year <- substr(basename, regexpr("_", basename) + 1, regexpr("_", basename) + 4)
  
  # Extract the month portion
  month <- substr(basename, regexpr("_", basename) + 5, regexpr("_", basename) + 6)  
  
  # Combine year and month
  paste(year, month, sep = "-")                                                
}


#' Divide strings by year and month.
#'
#' This function takes a vector of strings representing file paths and divides 
#' them into separate lists based on the year and month present in each file path.
#'
#' @param strings A character vector of file paths
#' @return A named list where each element contains a sublist of strings grouped 
#'         by year and month
#' 
#' @examples
#' strings <- c("E:/Thalis/04_processed_download/LST/LST_20220801_0179_069_006_2880.tif", 
#'              "E:/Thalis/04_processed_download/LST/LST_20220802_0179_069_020_3060.tif")
#' divide_strings_by_year_month(strings)
#'
#' @export

divide_strings_by_year_month <- function(strings) {
  year_month_lists <- list()
  for (string in strings) {
    year_month <- extract_year_month(string)
    if (!(year_month %in% names(year_month_lists))) {
      year_month_lists[[year_month]] <- list()
    }
    year_month_lists[[year_month]] <- c(year_month_lists[[year_month]], string)
  }
  names(year_month_lists) <- names(year_month_lists)
  return(year_month_lists)
}


#' Extract the identifier from a filepath.
#'
#' This function extracts the identifier from a given filepath. The identifier is
#' assumed to be the substring between the last "_" character and the ".tif" extension.
#'
#' @param filepath The full file path from which to extract the identifier
#' @return The extracted identifier

extract_identifier <- function(filepath) {
  # Extract the basename from the full file path
  basename <- basename(filepath)                                               
  substr(basename, regexpr("_", basename) + 1, nchar(basename) - 4)
}
 

#' Find the intersection of common identifiers between two lists of filepaths.
#'
#' This function takes two lists of filepaths and extracts the identifiers from 
#' each filepath. It then finds the common identifiers between the two lists and 
#' returns filtered versions of the input lists, containing only the filepaths 
#' with common identifiers.
#' 
#' @param list1 The first list of filepaths
#' @param list2 The second list of filepaths
#' @return A list containing the filtered versions of list1 and list2, with only 
#'         the filepaths having common identifiers
#'         
#' @examples
#' list1 <- c("E:/Thalis/04_processed_download/LST/LST_20220801_0179_069_006_2880.tif",
#'            "E:/Thalis/04_processed_download/LST/LST_20220802_0179_069_020_3060.tif")
#' list2 <- c("E:/Thalis/04_processed_download/LST/LST_20220802_0179_069_020_3060.tif",
#'            "E:/Thalis/04_processed_download/LST/LST_20220803_0179_069_024_3240.tif")
#' find_intersection(list1, list2)
#' @export
 
find_intersection <- function(list1, list2) {
  # Extract identifiers from list1 and list2
  identifiers1 <- sapply(list1, extract_identifier)
  identifiers2 <- sapply(list2, extract_identifier)
   
  # Find common identifiers
  common_identifiers <- intersect(identifiers1, identifiers2)
   
  # Filter list1 and list2 based on common identifiers
  filtered_list1 <- list1[sapply(identifiers1, function(id) id %in% common_identifiers)]
  filtered_list2 <- list2[sapply(identifiers2, function(id) id %in% common_identifiers)]
   
  return(list(filtered_list1, filtered_list2))
}


#' Build a virtual raster stack using GDAL
#'
#' This function builds a virtual raster stack using GDAL's gdalbuildvrt utility. 
#' It takes a list of input file paths and creates a virtual raster (.vrt) file 
#' that represents a stack of the input files.
#'
#' @param gdal_path The path to the GDAL installation directory
#' @param list A character vector of input file paths
#' @param list_source_dir The source directory where the input file paths are located
#' @param output_dir The directory where the output virtual raster file will be saved
#' @return None
#' 
#' @examples
#' gdal_path <- "C:/Program Files/GDAL"
#' list <- c("E:/Thalis/04_processed_download/LST/LST_20220801_0179_069_006_2880.tif",
#'           "E:/Thalis/04_processed_download/LST/LST_20220802_0179_069_020_3060.tif")
#' list_source_dir <- "path/to/source/directory"
#' output_dir <- "path/to/output/directory"
#' build_virtual_raster_stack(gdal_path, list, list_source_dir, output_dir)
#' @export

build_virtual_raster_stack <- function(gdal_path, list, list_source_dir, output_dir) {
  # Set the gdal_path variable
  gdal_path <- gdal_path
  
  # Build the full path to gdalbuildvrt
  buildvirtpath <- paste(gdal_path, "gdalbuildvrt.exe", sep = "\\")
  
  # Set the output directory variable
  output_dir <- output_dir
  
  # Create the path to the input file list text file
  path_to_input_list_txt_file <- file.path(output_dir, "input_file_list.txt")
  
  # Create the path to the output virtual raster file
  path_to_vrt_output_file <- file.path(output_dir, paste0(basename(list_source_dir), 
                                                          "_stack.vrt"))
  
  # Fill the text file with the paths of the found layers from which the VRT file 
  # should be created
  cat(list, sep = "\n", file = path_to_input_list_txt_file)
  
  # Parameters for gdalbuildvrt
  buildvrt_params <- c(
    "-separate",
    "-input_file_list", shQuote(path_to_input_list_txt_file),
    "-overwrite",
    shQuote(path_to_vrt_output_file)
  )
  
  # Create the gdalbuildvrt command
  buildvrt_cmd <- c(shQuote(buildvirtpath), buildvrt_params)
  
  # Execute the gdalbuildvrt command
  system(paste(buildvrt_cmd, collapse = " "))
  
  # Remove the list with input files
  file.remove(path_to_input_list_txt_file)
}


#' Create Monthly Composites of LST and NDVI
#'
#' This function creates monthly composites of the Land Surface Temperature (LST) 
#' and the Normalized Difference Vegetation Index (NDVI) from two list of input 
#' files. The monthly composites are based on the maximum index values of NDVI
#' and LST.
#'
#' @param gdal_path The path to the GDAL installation directory.
#' @param lst_dir The directory containing the LST files.
#' @param ndvi_dir The directory containing the NDVI files.
#' @param date_range  A character vector containing the start and end dates of 
#'                    the desired range in the format "YYYY-MM-DD".
#' @param output_dir The output directory where the composite files will be saved.
#'
#' @return None
#' @export
#'
#' @examples
#' gdal_path <- "C:\\Program Files\\QGIS 3.22.5\\bin"
#' lst_dir <- "E:/Thalis/04_processed_download/LST"
#' ndvi_dir <- "E:/Thalis/04_processed_download/NDVI"
#' date_range <- c("2022-02-01", "2023-04-30")
#' output_dir <- "E:/Thalis/05_output"
#' create_monthly_lst_ndvi_composites(gdal_path, lst_dir, ndvi_dir, date_range, output_dir)

createMonthlyComposites_ndvi_lst <- function(gdal_path, lst_dir, ndvi_dir, date_range, output_dir){
  gdal_path <- gdal_path
  
  # List all lst and ndvi files
  lst_list <- list.files(lst_dir , pattern = "*.tif$", full.names = TRUE)
  ndvi_list <- list.files(ndvi_dir , pattern = "*.tif$", full.names = TRUE)
  
  # Filter them by the date range
  lst_list <- filter_by_date_range(lst_list, date_range)
  ndvi_list <- filter_by_date_range(ndvi_list, date_range)
  
  if (length(lst_list) > 0 && length(ndvi_list) > 0){
    # Make sure that there is the same information in both datasets!
    intersection <- find_intersection(lst_list, ndvi_list)
    lst_list <- intersection[[1]]
    ndvi_list <- intersection[[2]]
    
    # Divide by year-month
    divided_lst_lists <- divide_strings_by_year_month(lst_list)
    divided_ndvi_lists <- divide_strings_by_year_month(ndvi_list)
    
    # Path to where the zip files will be extracted to 
    # (original downloads will be remained)
    output_dir <- output_dir
    
    path_to_lst_vrt_output_file <- file.path(output_dir, "LST_stack.vrt")
    path_to_ndvi_vrt_output_file <- file.path(output_dir, "NDVI_stack.vrt")
    
    # Iterate over the divided lists
    for (i in 1:length(divided_lst_lists)) {
      # Build virtual raster stacks for LST and NDVI
      build_virtual_raster_stack(gdal_path, unlist(divided_lst_lists[[i]], recursive = F), 
                                 lst_dir, output_dir)
      build_virtual_raster_stack(gdal_path, unlist(divided_ndvi_lists[[i]], recursive = F), 
                                 ndvi_dir, output_dir)
      
      # Check if the virtual raster files exist
      if (!file.exists(path_to_lst_vrt_output_file)) {
        stop("LST_stack.vrt does not exist!")
      }
      print(paste("LST_stack.vrt created! It contains", 
                  length(divided_lst_lists[[i]]), "layer.", sep = " "))
      
      if (!file.exists(path_to_ndvi_vrt_output_file)) {
        stop("NDVI_stack.vrt does not exist!")
      }
      print(paste("NDVI_stack.vrt created! It contains", 
                  length(divided_ndvi_lists[[i]]), "layer.", sep = " "))
      
      # Load the NDVI stack and find the maximum value
      ndvi_stack <- rast(path_to_ndvi_vrt_output_file)
      max_ndvi_composite <- app(ndvi_stack, fun = max, na.rm = TRUE)
      
      # Load the LST stack and find the maximum value
      lst_stack <- rast(path_to_lst_vrt_output_file)
      max_lst_composite <- app(lst_stack, fun = max, na.rm = TRUE)
      
      # Define the output file paths for the composite rasters
      path_to_lst_composite_output_file <- file.path(output_dir, paste0("LST_composite_", names(divided_lst_lists)[i], ".tif"))
      path_to_ndvi_composite_output_file <- file.path(output_dir, paste0("NDVI_composite_", names(divided_ndvi_lists)[i], ".tif"))
      
      # Write the composite rasters to disk
      terra::writeRaster(max_lst_composite, path_to_lst_composite_output_file, filetype = "GTiff", overwrite = TRUE)
      terra::writeRaster(max_ndvi_composite, path_to_ndvi_composite_output_file, filetype = "GTiff", overwrite = TRUE)
      
      print(paste0("Successfully created", paste0(" LST_composite_", names(divided_lst_lists)[i], ".tif"), "!"))
      print(paste0("Successfully created", paste0(" NDVI_composite_", names(divided_ndvi_lists)[i], ".tif"), "!"))
      
      file.remove(path_to_lst_vrt_output_file)
      file.remove(path_to_ndvi_vrt_output_file)
    }
  } else {
    print("No matching LST and NDVI tifs were found for the specified period. Therefore, no composites can be created.")
  }
}


#' Create Monthly Composites of NDII
#'
#' This function creates monthly composites of the NDII based on the monthly 
#' median index value from a list of input files.
#'
#' @param gdal_path The path to the GDAL installation directory.
#' @param ndii_dir The directory containing the NDII files.
#' @param date_range  A character vector containing the start and end dates of 
#'                    the desired range in the format "YYYY-MM-DD".
#' @param output_dir The output directory where the composite files will be saved.
#'
#' @return None

createMonthlyComposites_ndii <- function(gdal_path, ndii_dir, date_range, output_dir){
  gdal_path <- gdal_path
  ndii_list <- list.files(ndii_dir, pattern = "*.tif$", full.names = TRUE)
  ndii_list <- filter_by_date_range(ndii_list, date_range)
  divided_ndii_lists <- divide_strings_by_year_month(ndii_list)
  output_dir <- output_dir
  path_to_ndii_vrt_output_file <- file.path(output_dir, "NDII_stack.vrt")
    
  # Iterate over the divided lists
  for (i in 1:length(divided_ndii_lists)) {
    build_virtual_raster_stack(gdal_path, unlist(divided_ndii_lists[[i]], recursive = F), 
                               ndii_dir, output_dir)
    if (!file.exists(path_to_ndii_vrt_output_file)) {
      print("No matching NDII tifs were found for the specified period. Therefore, no composites can be created.")
      stop("NDII_stack.vrt does not exist!")
    }
    
    print(paste("NDII_stack.vrt created! It contains", length(divided_ndvi_lists[[i]]), 
                "layer.", sep = " "))
    ndii_stack <- rast(path_to_ndii_vrt_output_file)
    median_ndii_composite <- app(ndii_stack, fun = median, na.rm = TRUE)
    path_to_ndii_composite_output_file <- file.path(output_dir, paste0("NDII_composite_", names(divided_ndii_lists)[i], ".tif"))
    terra::writeRaster(median_ndii_composite, path_to_ndii_composite_output_file, filetype = "GTiff", overwrite = TRUE)
    print(paste0("Successfully created", paste0(" NDII_composite_", names(divided_ndii_lists)[i], ".tif"), "!"))
    file.remove(path_to_ndii_vrt_output_file)
  }
}


#' Calibrate Monthly Composites of LST and NDVI
#'
#' This function calibrates the monthly LST and NDVI composites from two lists  
#' of input files to account for difference between the S3 and MODIS sensor.
#'
#' @param composites_dir The directory containing the monthly LST and NDVI composites.
#' @param calibratedOutput_dir The output directory where the calibrated composite 
#'                             files will be saved.
#'
#' @return None

calibrateMonthlyComposites_ndvi_lst <- function(composites_dir, calibratedOutput_dir){
  
  # List all lst and ndvi composites
  lst_list <- list.files(composites_dir, pattern = "LST_composite_\\d{4}-\\d{2}\\.tif$", full.names = TRUE)
  ndvi_list <- list.files(composites_dir, pattern = "NDVI_composite_\\d{4}-\\d{2}\\.tif$", full.names = TRUE)
  
  # Iterate over the lists of composites
  for (i in 1:length(lst_list)) {
    # Load and calibrate composites
    ndvi_calibrated <- rast(ndvi_list[i]) + 0.05
    lst_calibrated <- rast(lst_list[i]) - 0.76
    
    # Create file names
    ndvi_file_name <- paste0(file_path_sans_ext(basename(ndvi_list[i])), 'calibrated')
    lst_file_name <- paste0(file_path_sans_ext(basename(lst_list[i])), 'calibrated')
    
    # Create output path
    output_dir <- calibratedOutput_dir
    
    ndvi_file_path <- paste0(output_dir, "/", ndvi_file_name)
    lst_file_path <- paste0(output_dir, "/", lst_file_name)
    
    # Write calibrated raster files
    writeRaster(ndvi_calibrated, paste0(ndvi_file_path, '.tiff'), overwrite = TRUE)
    writeRaster(lst_calibrated, paste0(lst_file_path, '.tiff'), overwrite = TRUE)
  }
}

#' Calibrate Monthly Composites of NDII
#'
#' This function calibrates the monthly NDII composites from two lists  of input 
#' files to account for difference between the S2 and MODIS sensor.
#'
#' @param composites_dir The directory containing the monthly NDII composites.
#' @param calibratedOutput_dir The output directory where the calibrated composite 
#'                             files will be saved.
#'
#' @return None

calibrateMonthlyComposites_ndii <- function(composites_dir, calibratedOutput_dir){
  # List all ndii composites
  ndii_list <- list.files(composites_dir, pattern = "NDII_composite_\\d{4}-\\d{2}\\.tif$", full.names = TRUE)
  
  # Iterate over the lists of composites
  for (i in 1:length(ndii_list)) {
    ndii_calibrated <- rast(ndii_list[i]) * 0.94 - 0.03
    ndii_file_name <- paste0(file_path_sans_ext(basename(ndii_list[i])), 'calibrated')
    output_dir <- calibratedOutput_dir
    ndii_file_path <- paste0(output_dir, ndii_file_name)
    writeRaster(ndii_calibrated, paste0(ndii_file_path, '.tiff'), overwrite = TRUE)
  }
}