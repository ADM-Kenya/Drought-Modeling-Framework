#==============================================================================#
#' Script to pre-process Sentinel-3 SY_2_SYN and SL_2_LST products             #
#'                                                                             #
#' This code is used to extract the downloaded zip SL_2_LST and SY_2_SYN       #
#' products. Then it processes each product to obtain the desired NDVI, LST    #
#' and NDII layers, together with the associated cloud cover information. It   #
#' then crops the layers to the shape of Kenya and masks out all cloud covered #
#' pixels. In the last step, the pre-processed NDVI, LST, and NDII .tif files  #
#' are stored in separate folders. The names of these files start with the     #
#' the layer name ("NDVI", "LST", "NDII"), followed by the sensing time and    #
#' the instance id of the corresponding Sentinel-3 product.                    #
#'                                                                             #
#' IMPORTANT: both function (for the LST and SYN products) use the focal       #
#' function from the raster package because the terra version was not          #
#' working.                                                                    #
#'                                                                             #
#' The user needs to specify:                                                  #
#' @param LST_dirpath The path to the raw Sentinel-3 SL_2_LST products.        #
#' @param SYN_dirpath The path to the raw Sentinel-3 SY_2_SYN products.        #
#' @param aoi_filepath The path to the shapefile of Kenya.                     #
#' @param output_dir The path to a output folder where three subfolders (NDVI, #
#'                   LST, NDII) for the processed products will be created and #
#'                   the files will be stored.                                 #
#'=============================================================================#


###############################  Functions  ####################################

# ---------------------------------------------------------------------------- #
#' Create Directory
#'
#' This function creates a directory at the specified path if it does not already exist.
#'
#' @param directory_path A character string specifying the path of the directory 
#'                       to be created.
#'
#' @return None
#'
#' @examples
#' create_directory("E:/Thalis/03_download/S3_SL_2_LST___20230524_095647")
#' 
#' @export
# ---------------------------------------------------------------------------- #

create_directory <- function(directory_path) {
  if (!dir.exists(directory_path)) {
    dir.create(directory_path)
    cat("Directory created:", directory_path, "\n")
  } else {
    cat("Directory already exists:", directory_path, "\n")
  }
}

# ---------------------------------------------------------------------------- #
#' Function, that processes sentinel-3 SL_2_LST products.
#'
#' First, the function extracts the downloaded zip SL_2_LST products, if necessary. 
#' Then it processes each product to obtain the desired LST and NDVI layers together 
#' with the associated cloud cover information. It then crops the layers to the 
#' given shapefile and mask out all the pixels covered by clouds. For each Sentinel-3 
#' product, one NDVI and one LST .tif are stored in the corresponding folders, if 
#' they contain information within the extent of the shapefile. The names of the 
#' .tif files start with the layer name ("LST" or "NDVI"), followed by the sensing 
#' starting time and the instance id of the corresponding sentinel product.
#'
#' @param LST_dirpath A character, path to directory containing zipped 
#'                    or unzipped sentinel-3 SL_2_LST products.
#' @param aoi_filepath An (optional) character, path to the shapefile 
#'                     to which the data is to be clipped.
#' @param output_dir A character, path to a folder where the 
#'                   processed products will be stored.
#' 
#' @return All processed NDVI and LST tifs in the output directory.
# ---------------------------------------------------------------------------- #

process_S3_SL_2_LST_products <- function(LST_dirpath, aoi_filepath, output_dir) {
  
  # Set the number of cores to be used for parallel execution
  num_cores <- detectCores()/2

  # Initialize parallel backend
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  ########## Data Preparation ##########
  
  # Set path to the directory containing S3_SL_2_LST zip files
  directory_path <- LST_dirpath

  # Path to where the zip files will be extracted to
  output_dir <- output_dir

  # Define output directory paths
  lst_dir <- file.path(output_dir, "LST")
  ndvi_dir <- file.path(output_dir, "NDVI")

  # Check if directories exist and create them if necessary
  create_directory(output_dir)
  create_directory(lst_dir)
  create_directory(ndvi_dir)

  # List all zip files of the input directory containing the sentinel-3 SL_2_LST products
  zip_files <- list.files(directory_path, pattern = ".zip$", full.names = TRUE)

  # If the directory contains any zip files..
  if (length(zip_files) > 0) {
    
    # .. check if all the zip files start with the correct pattern
    if(!all(grepl("^S3[AB]_SL_2_LST", basename(zip_files)))) {
      stop("One or more of the zip files are not sentinel-3 SL_2_LST products!")
    }
    cat("Extracting zip files before processing all files...\n")
    
    # .. and use foreach to extract each zip file parallel
    results <- foreach(zip_file = zip_files, .packages = "terra", .combine = 'c') %dopar% {
      # Initialize an empty list to collect results or status
      process_status <- list()
      
      # Extract name of data folder
      intermediate_path <- sub("\\.zip$", x = basename(zip_file), replacement = ".SEN3")
      
      # Extract the zip file
      unzip(
        zip_file, 
        exdir = directory_path, 
        files = c(
          paste0(intermediate_path, "/geodetic_in.nc"),
          paste0(intermediate_path, "/LST_ancillary_ds.nc"),
          paste0(intermediate_path, "/LST_in.nc"),
          paste0(intermediate_path, "/flags_in.nc")
        )
      )
      
      # Define path to the extracted flags_in.nc file
      flags_in_path <- file.path(directory_path, intermediate_path, "flags_in.nc")

      # Check and process only if flags_in.nc exists
      if(file.exists(flags_in_path)) {
        # Read the pointing_in dataset
        pointing_in_data <- rast(flags_in_path, subds = "pointing_in")
        
        # If the value 128 is found in pointing_in dataset, delete the folder
        if (any(values(pointing_in_data) == 128, na.rm = TRUE)) {
          # Delete the extracted folder
          unlink(file.path(directory_path, intermediate_path), recursive = TRUE)
          # Record status
          process_status <- c(process_status, paste(intermediate_path, "skipped and removed due to shifting"))
        } else {
          # Proceed with further processing if needed, or mark as processed
          process_status <- c(process_status, paste(intermediate_path, "processed"))
        }
      } else {
        process_status <- c(process_status, paste(intermediate_path, "flags_in.nc not found"))
      }

      # Delete the original zip file 
      file.remove(zip_file)
      
      return(print(process_status))
    }
    cat("Finished extracting zip files!", "\n")
    print(results)
  
  } else {
    cat("No zip files to extract!", "\n")
  }
  
  # -------------------------------------------------------------------------- #
  # The passed directory now contains only extracted zip files
  
  # Create list of extracted directories we want to iterate over
  dir_list <- list.dirs(path = directory_path, recursive = FALSE) 

  #' Extract the instance ID from a given string.
  #'
  #' This function extracts the instance ID from a string based on specific patterns.
  #'
  #' @param string The input string from which the instance ID will be extracted.
  #' @return The instance ID extracted from the string, or NA if no match is found.
  #' @examples
  #' string1 <- "S3A_SL_2_LST____20160906T065849_20160906T083948_20180929T155649_6059_008_220______LR1_R_NT_003.SEN3"
  #' string2 <- "S3A_SL_2_LST____20230130T073307_20230130T073607_20230131T164133_0180_095_049_2880_PS1_O_NT_004"
  #' instance_id1 <- extract_instance_id(string1)
  #' instance_id2 <- extract_instance_id(string2)
  #' print(instance_id1)  # "6059_008_220_____"
  #' print(instance_id2)  # "0180_095_049_2880"
  #'
  #' @export
  # -------------------------------------------------------------------------- #
  
  extract_instance_id <- function(string) {
    # Pattern for instrument data products disseminated in "stripes"
    pattern_stripes <- "_(\\d{4}_\\d{3}_\\d{3}_____)"

    # Pattern for instrument data products disseminated in "frames"
    pattern_frames <- "_(\\d{4}_\\d{3}_\\d{3}_\\d{4})"

    # Check if the string matches the pattern for "stripes"
    if (grepl(pattern_stripes, string)) {
      instance_id <- regmatches(string, regexpr(pattern_stripes, string))
      return(sub("_", "", instance_id))
    }

    # Check if the string matches the pattern for "frames"
    if (grepl(pattern_frames, string)) {
      instance_id <- regmatches(string, regexpr(pattern_frames, string))
      return(sub("_", "", instance_id))
    }

    # Return NA if no match is found
    return(NA)
  }

  count_new_tifs <- 0
  cat("Processing files...", "\n")
  # Loop through extracted zip files and process data to get the LST and NDVI layers
  new_tifs <- foreach(subdir = dir_list, .packages = "terra", .combine = c) %dopar% {

    # Check if the products to be generated already exist
    # Generate name for the TIF to be output with the information we need
    input_string <- basename(subdir)
    date <- gsub(".*_LST____(\\d{8}).*", "\\1", input_string)
    instance_id <- extract_instance_id(input_string)
    
    ########## LST Extraction ##########

    # Save LST raster as tif in the LST directory
    output_file_lst <- file.path(lst_dir, paste0(paste("/LST", date, instance_id, sep = "_"), ".tif"))
    output_file_ndvi <- file.path(ndvi_dir, paste0(paste("/NDVI", date, instance_id, sep = "_"), ".tif"))

    if (!file.exists(output_file_lst) || !file.exists(output_file_ndvi)) {
      # Read coordinates data
      longitude = rast(file.path(subdir, "geodetic_in.nc"), subds = "longitude_in")
      latitude = rast(file.path(subdir, "geodetic_in.nc"), subds = "latitude_in")

      # Read LST data
      lst = rast(file.path(subdir, "/LST_in.nc"), subds = "LST")

      # Create empty raster
      e = ext(cbind(values(longitude), values(latitude)))
      r = rast(e, ncol = ncol(longitude), nrow = nrow(longitude))

      # Assign LST values to coordinates
      lst_array = cbind(values(longitude), values(latitude), values(lst))

      # Save max LST value
      max_lst <- max(lst_array[,3], na.rm = TRUE)

      # Set all NA values to 9999, so that we can later distinguish 
      # which pixels were already NA before rasterizing
      lst_array[,3][is.na(lst_array[,3])] <- 9999

      # Create raster
      lst_raster = terra::rasterize(lst_array[,1:2], r, lst_array[,3], fun = 'first') # We are using the function "first" so that we get the matching cloud mask information later on
      crs(lst_raster) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

      # Apply focal function to interpolate missing values caused by rasterization
      lst_raster <- terra::focal(
        lst_raster, w=matrix(
          1, nrow=3, ncol=3), fun=mean, 
          na.policy = "only", na.rm = TRUE
      )

      # Set all values where a no data pixel with the value of 9999 was included 
      # in the focal function to NA
      lst_raster[lst_raster > max_lst & lst_raster <= 9999] <- NA

      # Crop raster to the extent of the shapefile
      if (!is.null(aoi_filepath)){
        shapefile <- vect(aoi_filepath)
        crs(shapefile) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
        cropped_lst_raster <- crop(lst_raster, ext(shapefile))
        cropped_lst_raster <- mask(cropped_lst_raster, shapefile)
      }
      
      ########## NDVI Extraction ##########

      # Read NDVI data
      ndvi = rast(file.path(subdir, "LST_ancillary_ds.nc"), subds = "NDVI")

      # Assign NDVI values to coordinates
      ndvi_array = cbind(values(longitude), values(latitude), values(ndvi))

      # Set all NA values to 9999, so that we can later distinguish 
      # which pixels were already NA before rasterizing
      ndvi_array[,3][is.na(ndvi_array[,3])] <- 9999

      # Create raster
      ndvi_raster = terra::rasterize(ndvi_array[,1:2], r, ndvi_array[,3], fun = 'first') # We are using the function "first" so that we get the matching cloud mask information later on
      crs(ndvi_raster) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

      # Apply focal function to interpolate missing values caused by rasterization
      ndvi_raster <- terra::focal(
        ndvi_raster, w=matrix(1, nrow=3, ncol=3), 
        fun=mean, na.policy="only", na.rm = TRUE
      )

      # Set all values where a no data pixel with the value of 9999 was included 
      # in the focal function, all values > 1 and all values < -1 to NA
      ndvi_raster[(ndvi_raster > 1) | (ndvi_raster < -1)] <- NA #insert 'or' criteria

      # Crop raster to the extent of the shapefile
      if (!is.null(aoi_filepath)){
        cropped_ndvi_raster <- crop(ndvi_raster, ext(shapefile))
        cropped_ndvi_raster <- mask(cropped_ndvi_raster, shapefile)
      }
      
      ########## Cloud Cover Information Extraction I ##########

      # Extract the information about cloud cover from the "confidence_in" 
      # variable of the "flags_in.nc" file to generate cloud mask
      confidence_in = rast(file.path(subdir, "flags_in.nc"), subds = "confidence_in")

      # Assign confidence values to coordinates
      confidence_in_array = cbind(
        values(longitude), 
        values(latitude), 
        values(confidence_in)
      )

      # Create raster
      confidence_in_raster = terra::rasterize(
        confidence_in_array[,1:2], 
        r, confidence_in_array[,3], 
        fun = "first") # We are using the function "first" so that we get the matching LST information
      crs(confidence_in_raster) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
      
      # Apply focal function to fill missing values 
      # (we are using the modal function, since it is an categorical variable)
      confidence_in_raster <- raster::focal(
        confidence_in_raster, w = matrix(1, nrow = 3, ncol = 3), 
        fun = modal, Naonly = TRUE, na.rm = TRUE)
      # confidence_in_raster <- terra::focal(confidence_in_raster, w = matrix(1, nrow = 3, ncol = 3), fun = modal, na.policy = "only", na.rm = TRUE)

      # Reclassify all pixel: Pixel with values >=16384 && < 32678 are influenced 
      # by clouds and set to 1
      reclass_mat_1 <- matrix(
        c(0, 16384, 0, 16384, 32768, 1, 32768, Inf, 0), ncol = 3, byrow = TRUE
      )
      confidence_in_raster <- classify(confidence_in_raster, reclass_mat_1)

      # Crop raster to the extent of the shapefile
      if (!is.null(aoi_filepath)){
        cropped_confidence_in_raster <- crop(confidence_in_raster, ext(shapefile))
        cropped_confidence_in_raster <- mask(cropped_confidence_in_raster, shapefile)
      }
      
      ########## Cloud Cover Information Extraction II ##########

      # Extract the information about cloud cover from the "bayes_in" variable 
      # of the "flags_in.nc" file to generate cloud mask
      bayes_in = rast(file.path(subdir, "flags_in.nc"), subds = "bayes_in")

      # Assign confidence values to coordinates
      bayes_in_array = cbind(values(longitude), values(latitude), values(bayes_in))

      # Create raster
      bayes_in_raster = terra::rasterize(bayes_in_array[,1:2], r, bayes_in_array[,3], fun = "first")
      crs(bayes_in_raster) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

      # Apply focal function
      bayes_in_raster <- raster::focal(
        bayes_in_raster, w = matrix(1, nrow = 3, ncol = 3),
        fun = modal, Naonly = TRUE, na.rm = TRUE
      )
      # bayes_in_raster <- terra::focal(bayes_in_raster, w=matrix(1, nrow=3, ncol=3), fun=modal, na.policy="only", na.rm=TRUE)

      # Reclassify all pixel: Pixel with values == 2 are influenced by clouds and set to 1
      reclass_mat_2 <- matrix(c(0, 0, 2, 1), ncol = 2, byrow = TRUE)
      bayes_in_raster <- classify(bayes_in_raster, reclass_mat_2)

      # Crop raster to the extent of the shapefile
      if (!is.null(aoi_filepath)){
        cropped_bayes_in_raster <- crop(bayes_in_raster, ext(shapefile))
        cropped_bayes_in_raster <- mask(cropped_bayes_in_raster, shapefile)
      }
      
      ########## Cloud Masking ##########
      
      # Cloud masking
      cloud_mask <- (cropped_confidence_in_raster == 1 | cropped_bayes_in_raster == 1)

      cropped_lst_raster_filtered <- cropped_lst_raster
      cropped_lst_raster_filtered[cloud_mask] <- NA

      cropped_ndvi_raster_filtered <- cropped_ndvi_raster
      cropped_ndvi_raster_filtered[cloud_mask] <- NA
      
      ########## Writing LST and NDVI Raster ##########

      # If all values of the cropped and filtered LST or NDVI raster are NA, 
      # we will not proceed, because the layers have no overlap with our area of 
      # interest and therefore no information for us
      
      if (!all(is.na(values(cropped_lst_raster_filtered))) || !all(is.na(values(cropped_ndvi_raster_filtered)))) {
        
        # If they do contain information, save LST and NDVI raster as tif in the 
        # corresponding directories
        writeRaster(
          cropped_lst_raster_filtered, output_file_lst, 
          filetype = "GTiff", overwrite = TRUE
        )
        writeRaster(
          cropped_ndvi_raster_filtered, output_file_ndvi, 
          filetype = "GTiff", overwrite = TRUE
        )
        1 # Return 1 if the condition is met
      }
      else {
        print(0) # Return 0 if the condition is not met
      }
    }
  }
  count_new_tifs <- sum(new_tifs != 0) # Sum all the 1s returned by the foreach loop
  stopCluster(cl)
  cat("Finished processing sentinel-3 S3_SL_2_LST data!", "\n")
  cat("New LST and NDVI tifs created:", count_new_tifs, "\n")
}

# ---------------------------------------------------------------------------- #
#' Function, that processes sentinel-3 SY_2_SYN products.
#'
#' First, the function extracts the downloaded zip SY_2_SYN products, if necessary. 
#' Then it processes each product to obtain the desired S3 and S5 band reflectances 
#' together with the associated cloud cover information. It then crops the layers 
#' to the given shapefile and mask out all the pixels covered by clouds. For each 
#' Sentinel 3 product an NDII tif is stored in the corresponding folder, if it 
#' contains information within the shapefile. The name of the tif files starts 
#' with 'NDII' followed by the sensing starting time and the instance id of the 
#' corresponding sentinel product.
#'
#' @param SYN_dirpath A character, path to directory containing zipped or 
#'                    unzipped sentinel-3 SY_2_SYN products.
#' @param aoi_filepath (optional) A character, path to the shapefile 
#'                     to which the data is to be clipped.
#' @param output_dir A character, path to a folder where 
#'                   the processed products will be stored.
#' 
#' @return All processed NDII tifs in the output directory.
#'
#' @export
# ---------------------------------------------------------------------------- #

process_S3_SY_2_SYN_products <- function(SYN_dirpath, aoi_filepath, output_dir){
  
  # Set the number of cores to be used for parallel execution
  num_cores <- detectCores()/2
  
  # Initialize parallel backend
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  ########## Data Preparation ##########
  
  # Set path to the directory containing zipped or unzipped S3_SL_2_LST files
  directory_path <- SYN_dirpath
  
  # Path to where the zip files will be extracted to 
  output_dir <- output_dir
  
  # Define output directory paths
  ndii_dir <- file.path(output_dir, "NDII")
  
  # Check if directories exist and create them if necessary
  create_directory(output_dir)
  create_directory(ndii_dir)
  
  # List all zip files in the directory containing the sentinel-3 SY_2_SYN products
  zip_files <- list.files(directory_path, pattern = ".zip$", full.names = TRUE)
  
  # If the directory contains any zip files..
  if (length(zip_files) > 0) {
    # .. check if all the zip files start with the correct pattern
    if(!all(grepl("^S3[AB]_SY_2_SYN", basename(zip_files)))) {
      stop("One or more of the zip files are not sentinel-3 SY_2_SYN products!")
    }
    cat("Extracting zip files before processing all files...\n")
    
    # .. and use foreach to extract each zip file parallel
    foreach(zip_file = zip_files) %dopar% {

      # Extract name of the data folder
      intermediate_path <- sub("\\.zip$", x = basename(zip_file), replacement = ".SEN3")
      
      # Extract the zip file
      unzip(
        zip_file, 
        exdir = directory_path, 
        files = c(
          paste0(intermediate_path, "/geolocation.nc"),
          paste0(intermediate_path, "/Syn_S3N_reflectance.nc"),
          paste0(intermediate_path, "/Syn_S5N_reflectance.nc"),
          paste0(intermediate_path, "/flags.nc")
        )
      )
  
      # Delete the original zip file
      file.remove(zip_file)

    }
    cat("Finished extracting zip files!", "\n")
  } else {
    cat("No zip files to extract!", "\n")
  }
  
  # The passed directory now contains only extracted zip files with the necessary .nc files
  
  # Create list of extracted directories we want to iterate over
  dir_list <- list.dirs(path = directory_path, recursive = FALSE) # List all files in the processed directory
  
  # -------------------------------------------------------------------------- #
  #' Extract the instance ID from a given string.
  #'
  #' This function extracts the instance ID from a string based on specific patterns.
  #'
  #' @param string The input string from which the instance ID will be extracted.
  #' @return The instance ID extracted from the string, or NA if no match is found.
  #' @examples
  #' string1 <- "S3A_SL_2_LST____20160906T065849_20160906T083948_20180929T155649_6059_008_220______LR1_R_NT_003.SEN3"
  #' string2 <- "S3A_SL_2_LST____20230130T073307_20230130T073607_20230131T164133_0180_095_049_2880_PS1_O_NT_004"
  #' instance_id1 <- extract_instance_id(string1)
  #' instance_id2 <- extract_instance_id(string2)
  #' print(instance_id1)  # "6059_008_220_____"
  #' print(instance_id2)  # "0180_095_049_2880"
  #' 
  #' @export
  # -------------------------------------------------------------------------- #
  
  extract_instance_id <- function(string) {
    # Pattern for instrument data products disseminated in "stripes"
    pattern_stripes <- "_(\\d{4}_\\d{3}_\\d{3}_____)"
    
    # Pattern for instrument data products disseminated in "frames"
    pattern_frames <- "_(\\d{4}_\\d{3}_\\d{3}_\\d{4})"
    
    # Check if the string matches the pattern for "stripes"
    if (grepl(pattern_stripes, string)) {
      instance_id <- regmatches(string, regexpr(pattern_stripes, string))
      return(sub("_", "", instance_id))
    }
    
    # Check if the string matches the pattern for "frames"
    if (grepl(pattern_frames, string)) {
      instance_id <- regmatches(string, regexpr(pattern_frames, string))
      return(sub("_", "", instance_id))
    }
    
    # Return NA if no match is found
    return(NA)
  }
  
  count_new_tifs <- 0
  cat("Processing files...", "\n")
  
  # Loop through extracted zip files and process data to get a NDII layers
  new_tifs <- foreach(subdir = dir_list, .packages = c("terra", "raster"), .combine = c) %dopar% {
    
    # Generate name for the TIF to be output with the information we need
    input_string <- basename(subdir)
    date <- gsub(".*_SYN____(\\d{8}).*", "\\1", input_string)
    instance_id <- extract_instance_id(input_string)
    
    ########## NDII Calculation ##########
    
    # Save LST raster as tif in the LST directory
    output_file_ndii <- file.path(ndii_dir, paste0(paste("NDII", date, instance_id, sep = "_"), ".tif"))
    
    # Check if the products to be generated already exist
    if (!file.exists(output_file_ndii)) {
      # Read coordinates data
      longitude = rast(file.path(subdir, "geolocation.nc"), subds = "lon")
      latitude = rast(file.path(subdir, "geolocation.nc"), subds = "lat")
      
      # Read S3 and S5 band
      s3_band = rast(file.path(subdir, "Syn_S3N_reflectance.nc"), subds = "SDR_S3N")
      s5_band = rast(file.path(subdir, "Syn_S5N_reflectance.nc"), subds = "SDR_S5N")
      
      # Compute NDII values
      ndii_band <- (s3_band - s5_band)/(s3_band + s5_band)
      
      # Create empty raster
      e = ext(cbind(values(longitude), values(latitude)))
      r = rast(e, ncol = ncol(longitude), nrow = nrow(longitude))
      
      # Assign NDII values to coordinates
      ndii_array = cbind(values(longitude), values(latitude), values(ndii_band))
      
      # Save max NDII value
      max_ndii <- max(ndii_array[,3], na.rm = TRUE)
      
      # Set all NA values to 9999, so that we can later distinguish which pixels 
      # were already NA before rasterizing
      ndii_array[,3][is.na(ndii_array[,3])] <- 9999
      
      # Create raster
      ndii_raster = terra::rasterize(ndii_array[,1:2], r, ndii_array[,3], fun = 'first') # We are using the function "first" so that we get the matching cloud mask information later on
      crs(ndii_raster) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
      
      # Apply focal function to interpolate missing values caused by rasterization
      ndii_raster <- terra::focal(ndii_raster, w=matrix(1, nrow=3, ncol=3), 
                                  fun=mean, na.policy="only", na.rm = TRUE) 
      
      # Set all values where a no data pixel with the value of 9999 was included 
      # in the focal function to NA
      ndii_raster[ndii_raster > max_ndii & ndii_raster <= 9999] <- NA
      
      # Crop raster to the extent of the shapefile
      if (!is.null(aoi_filepath)){
        shapefile <- vect(aoi_filepath)
        crs(shapefile) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
        cropped_ndii_raster <- crop(ndii_raster, ext(shapefile))
        cropped_ndii_raster <- mask(cropped_ndii_raster, shapefile)
      }
      
      ########### Cloud Cover Information Extraction ##########
      
      # Read cloud flags
      cloud_flags = rast(file.path(subdir, "flags.nc"), subds = "CLOUD_flags")
      
      # Assign cloud flags values to coordinates
      cloud_flags_array = cbind(values(longitude), values(latitude), values(cloud_flags))
      
      # Create raster
      cloud_flags_raster = terra::rasterize(cloud_flags_array[,1:2], r, 
                                            cloud_flags_array[,3], fun = "first") # We are using the function "first" so that we get the matching LST information
      crs(cloud_flags_raster) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
      
      # Apply focal function to fill missing values 
      # (we are using the modal function, since it is an categorical variable)
      cloud_flags_raster <- raster::focal(cloud_flags_raster, w = matrix(1, nrow = 3, ncol = 3), 
                                          fun = modal, NAonly = TRUE, na.rm = TRUE)
      # cloud_flags_raster <- terra::focal(cloud_flags_raster, w=matrix(1, nrow=3, ncol=3), fun=modal, na.policy="only", na.rm=TRUE)
      
      # Crop raster to the extent of the shapefile
      if (!is.null(aoi_filepath)){
        cropped_cloud_flags_raster <- crop(cloud_flags_raster, ext(shapefile))
        cropped_cloud_flags_raster <- mask(cropped_cloud_flags_raster, shapefile)
      }
      
      ############ Cloud Masking ##########
      
      # Cloud masking
      cloud_mask <- (cropped_cloud_flags_raster > 1)
      cropped_ndii_raster_filtered <- cropped_ndii_raster
      cropped_ndii_raster_filtered[cloud_mask] <- NA
      
      # If all values of the cropped and filtered NDII raster are NA, we will not 
      # proceed, because the layer has no overlap with our area of interest and 
      # therefore no information for us
      
      if (!all(is.na(values(cropped_ndii_raster_filtered)))) {
        # If they do contain information, save LST and NDVI raster as tif in the 
        # corresponding directories
        terra::writeRaster(cropped_ndii_raster_filtered, output_file_ndii, 
                           filetype = "GTiff", overwrite = TRUE)
        1 # Return 1 if the condition is met
      } else {
        0 # Return 0 if the condition is not met
      }
    }
  }
  count_new_tifs <- sum(new_tifs != 0) # Sum all the 1s returned by the foreach loop
  stopCluster(cl)
  cat("Finished processing sentinel-3 S3_SY_2_SYN data!", "\n")
  cat("New NDII tifs created:", count_new_tifs, "\n")
}