#==============================================================================#
#' Script with functions for pre-processing the Sentinel-3 data and            #
#' calculation of monthly anomalies. Additionally, it can be used to           #
#' re-calculate the MODIS and Sentinel-3 baselines and create new              #
#' monthly anomalies.                                                          #
#'                                                                             #
#' This code includes functions to pre-process the raw Sentinel-3 data and     #
#' create monthly composites of NDVI, NDII, and LST. Additionally, the         #
#' composites are calibrated to match the MODIS composites. Then, the script   #
#' can be used to compute monthly anomalies of the three variables.            #
#' Optionally, there is a function to re-calculate the baselines every few     #
#' years to keep them up-to-date.If the baselines and anomalies should be      #
#' re-calculated, the user needs to comment out the function                   #
#' "S3_calculateAnomalies" and uncomment the function                          # 
#' "MODIS_S3_calculateBaselinesAnomalies".                                     #
#'                                                                             #
#' The user needs to specify:                                                  #
#'                                                                             #
#' 1. Pre-Processing Parameters:                                               #
#' @param aoi_filepath The path to the shapefile of Kenya.                     #
#' @param LST_dir      The path to the raw Sentinel-3 SL_2_LST products.       #
#' @param SYN_dir      The path to the raw Sentinel-3 SY_2_SYN products.       #
#' @param S3_preProcessed_dir The path to the directory containing the         #
#'                            pre-processed S3 LST and SYN data.               #
#'                                                                             #
#' 2. Composite creation and Calibration Parameters:                           #
#' @param gdal_dir    The path to the GDAL installation directory.             #
#' @param date_range  The date range for which the                             #
#'                    monthly composites are created.                          #
#' @param S3_preprocessed_dir      The directory containing the                #
#'                                 pre-processes S3 folders                    #
#' @param S3_monthlyComposites_dir The directory containing the S3 monthly     # 
#'                                 composites (NDVI, NDII, LST).               #                                             
#' @param S3_monthlyCalibratedComposites_dir The directory containing the S3   #
#'                                           monthly calibrated composites.    #
#'                                                                             #
#' 3. S3 Monthly Anomalies Calculation Parameters:                             #
#' @param temp_dir The path to the temporary directory.                        #
#' @param MODIS_S3_Baselines_dir The directory containing the combined MODIS   #
#'                               and S3 baselines (mean/ STD).                 #
#' @param MODIS_S3_monthlyAnomalies_dir The directory containing the monthly   #
#'                                      NDVI, NDII, and LST anomalies.         #                 
#'                                                                             #
#' 4. Re-Calculation of MODIS/ S3 Baselines and Anomalies Parameters:          #
#' @param MODIS_monthlyComposites_dir The directory containing the MODIS       #
#'                                    monthly composites (NDVI, NDII, LST).    #
#'                                                                             #
#' 5. Calculation of 3-Monthly SPI from TAMSAT Rainfall Data:                  #
#' @param TAMSAT_InputFile     The TAMSAT Monthly Rainfall data of the AOI     # 
#' @param SPI_OutputDir        The output directory for the SPI results        #
#'                                                                             #
#==============================================================================#


#==============================================================================#
##################  Install and load required packages  ########################

packages <- c(
  "raster", "terra", "doParallel", "foreach", "tools", 
  "furrr", "R.utils", "Rcpp", "tcltk", "ncdf4",
  "SPEI", "lubridate", "ClusterR", "snow", "parallel"
)
for (i in packages) {
  if(!require(i, character.only = TRUE)) {
    install.packages(i)
    library(i, character.only = TRUE)
  }
}


#==============================================================================#
##  Custom functions for using the correct "choose()" function depends on OS  ##

chooseDirCustom <- function(caption = "Select Directory") {
  if (exists("choose.dir", where = "package:utils")) {
    # Use utils::choose.dir
    selectedDir <- utils::choose.dir(caption=caption)
  } else {
    # Fallback to tcltk::tk_chooseDirectory
    if (!requireNamespace("tcltk", quietly = TRUE)) {
      stop("The 'tcltk' package is needed but not available.")                 
    }
    selectedDir <- tcltk::tk_choose.dir(caption=caption)
  }
  
  return(selectedDir)
}

chooseFilesCustom <- function(caption = "Select Files") {
  if (exists("choose.files", where = "package:utils")) {
    # Use utils::choose.files
    selectedFiles <- utils::choose.files(caption=caption)
  } else {
    # Fallback to tcltk::tk_getOpenFile
    if (!requireNamespace("tcltk", quietly = TRUE)) {
      stop("The 'tcltk' package is needed but not available.")
    }
    selectedFiles <- tcltk::tk_choose.files(caption=caption)
  }
  
  return(selectedFiles)
}

#==============================================================================#
##########################  Step 1: Pre-processing  ############################
# ---------------------------------------------------------------------------- #
#' Pre-Process Sentinel-3 SL_2_LST and SY_2_SYN products
#' 
#' This function is used to extract the downloaded SL_2_LST/ SY_2_SYN zip files 
#' to obtain the desired LST, NDVI, and NDII layers, crop the layers to the 
#' extent of Kenya and mask out all cloud covered pixels. The output is then 
#' stored in the respective NDVI, LST, and NDII folder.
#' 
#' @param aoi         The path to the shapefile of Kenya
#'                    to which the data is clipped.
#' @param LST_dirpath The directory containing the zipped 
#'                    or unzipped S3 SL_2_LST products.
#' @param SYN_dirpath The directory containing zipped or 
#'                    unzipped sentinel-3 SY_2_SYN products.
#' @param output_dir  The path to a folder where the processed 
#'                    products will be stored.
# ---------------------------------------------------------------------------- #

### Load the step 1 script
step1_script <- chooseFilesCustom("Choose R-Script 'S3_PreProcessing.R'")
source(step1_script)

### SET DIRECTORIES
aoi <- chooseFilesCustom("Select shapefile with AOI of Kenya")
LST_dir <- chooseDirCustom("Select Folder with raw Sentinel-3 SL_2_LST products")
SYN_dir <- chooseDirCustom("Select Folder with raw Sentinel-3 SY_2_SYN products")
S3_preProcessed_dir <- chooseDirCustom(
  "Select Folder to store the pre-processed Sentinel-3 files"
)

#### RUN FUNCTIONS
process_S3_SL_2_LST_products(LST_dir, aoi, S3_preProcessed_dir)
process_S3_SY_2_SYN_products(SYN_dir, aoi, S3_preProcessed_dir)

# Clean up the environment (keep directories if needed for the next step)
to_keep1 <- c(
  "packages", "chooseDirCustom", "chooseFilesCustom",
  "step1_script",
  "aoi", "LST_dir", "SYN_dir", "S3_preProcessed_dir"
)
# Items to be removed
to_remove1 <- setdiff(ls(all.names = TRUE), to_keep1)

# Optionally print items to be removed for verification
print(paste("Removing:", to_remove1))

# Proceed with removal if verified
rm(list = to_remove1)


#==============================================================================#
#################  Step 2: Composite Creation and Calibration  #################
# ---------------------------------------------------------------------------- #
#' Create Monthly NDVI, NDII, and LST Composites
#' 
#' These functions creates monthly composites of the NDVI, 
#' NDII, and LST from three lists of input files.
#' 
#' @param gdal_path         The path to the GDAL installation directory.
#' @param preprocessed_dir  The directory containing the pre-processed
#'                          folders of LST, NVDI, and NDII
#' @param date_range        A character vector with the start and end dates 
#'                          of the desired range in the format "YYYY-MM-DD".
#' @param output_dir        The output directory where the monthly
#'                          composite files will be saved.
# ---------------------------------------------------------------------------- #

### Load the composites and calibration script
step2_script <- chooseFilesCustom("Choose R-Script 'S3_Monthly_Composites_MODIS_Calibration.R'")
source(step2_script)  # Re-run this line when adjust function during execution

### SET VARIABLES and DIRECTORIES

# 'date_range' variable
if (.Platform$OS.type == "windows") {
  date_range <- winDialogString(
    "Enter the date range for which the monthly composites are created. Use format: YYYY-MM-DD, YYYY-MM-DD", ""
  )
} else {
  date_range <- readline(
    prompt = "Enter the date range for which the monthly composites are created. Use format: YYYY-MM-DD, YYYY-MM-DD: ")
}

# 'QGIS' folder
if (.Platform$OS.type == "windows") {
  gdal_dir <- chooseDirCustom("Select Folder of the QGIS bin.")  # For Windows
} else {
  gdal_dir <- "/usr/bin"  # Default for Linux and other non-Windows
}

# 'S3_monthlyComposites_dir' folder
S3_monthlyComposites_dir <- chooseDirCustom(
  "Select Folder to store the monthly S3 NDVI, NDII, and LST composites"
)

#### Execute 'MonthlyComposites' FUNCTIONS
createMonthlyComposites_ndvi_lst(gdal_dir, S3_preProcessed_dir, date_range, S3_monthlyComposites_dir)
createMonthlyComposites_ndii(gdal_dir, S3_preProcessed_dir, date_range, S3_monthlyComposites_dir)

# ---------------------------------------------------------------------------- #
#' Calibrate Monthly NDVI, NDII, and LST Composites
#' 
#' This function calibrates the monthly composites of the NDVI, NDII, and LST 
#' from three lists of input files.
#' 
#' @param composites_dir       The directory containing the monthly composites
#'                             of LST, NDII and NDVI composites.
#' @param calibratedOutput_dir The output directory where the calibrated 
#'                             composite files will be saved.
# ---------------------------------------------------------------------------- #

#### SET DIRECTORIES
S3_calibratedMonthlyComposites_dir <- chooseDirCustom(
  "Select Folder to store the calibrated monthly composites (NDVI, NDII, and LST)"
)

#### RUN 'Calibration' FUNCTIONS
calibrateMonthlyComposites_ndvi_lst(S3_monthlyComposites_dir, S3_calibratedMonthlyComposites_dir)
calibrateMonthlyComposites_ndii(S3_monthlyComposites_dir, S3_calibratedMonthlyComposites_dir)

# Clean up the environment (keep directories if needed for the next step)
to_keep2 <- c(
  "packages", "chooseDirCustom", "chooseFilesCustom",
  "step1_script", "step2_script",
  "aoi", "LST_dir", "SYN_dir", "S3_preProcessed_dir",
  "gdal_dir", "date_range", "S3_monthlyComposites_dir", "S3_calibratedMonthlyComposites_dir"
)

# Items to be removed
to_remove2 <- setdiff(ls(all.names = TRUE), to_keep2)

# Optionally print items to be removed for verification
print(paste("Removing:", to_remove2))

# Proceed with removal if verified
rm(list = to_remove2)


#==============================================================================#
###################  Step 3: Monthly Anomalies Calculation #####################
# ---------------------------------------------------------------------------- #
#' Calculate monthly Anomalies for NDVI, NDII, and LST
#' 
#' This function calculates baseline mean and std to compute 
#' the monthly anomalies of ndvi, ndii, and lst and saves 
#' the results in the defined output folder.
#' 
#' @param baselines_dir The directory contains the baseline mean and std files
#' @param S3_dir The directory containing the Sentinel-3 monthly composites.
#' @param out_dir The directory where the output files will be saved.
# ---------------------------------------------------------------------------- #

# Load the anomalies calculation script
step3_script <- chooseFilesCustom("Choose R-Script 'S3_AnomaliesCalculation_no_Baseline_Recalculation.R'")
source(step3_script)

# Set location of 'r_temp' folder (for Windows only) if your "C:" Drive is full 
temp_dir <- chooseDirCustom(caption = "Select temporary directory.")

#### SET DIRECTORIES
MODIS_S3_Baselines_dir <- chooseDirCustom("Select Folder containing the mean and std baseline files.")
MODIS_S3_monthlyAnomalies_dir <- chooseDirCustom("Select Folder with old (to update) monthly anomalies.")
MODIS_S3_monthlyAnomalies_Update_dir <- chooseDirCustom("Select Folder to store updated monthly anomalies")

#### RUN FUNCTION
S3_calculateAnomalies(
  baselines_dir = MODIS_S3_Baselines_dir, S3_dir = S3_calibratedMonthlyComposites_dir, old_anomalies = MODIS_S3_monthlyAnomalies_dir, 
  out_dir = MODIS_S3_monthlyAnomalies_Update_dir
)

# Clean up the environment (keep directories if needed for the next step)
to_keep3 <- c(
  "packages", "chooseDirCustom", "chooseFilesCustom",
  "step1_script", "step2_script", "step3_script",
  "aoi", "LST_dir", "SYN_dir", "S3_preProcessed_dir",
  "gdal_dir", "date_range", "S3_monthlyComposites_dir","S3_calibratedMonthlyComposites_dir", 
  "MODIS_S3_Baselines_dir", "MODIS_S3_monthlyAnomalies_dir"
)

# Items to be removed
to_remove3 <- setdiff(ls(all.names = TRUE), to_keep3)

# Optionally print items to be removed for verification
print(paste("Removing:", to_remove3))

# Proceed with removal if verified
rm(list = to_remove3)


#==============================================================================#
#########  Step: 4. Re-Calculation of Baselines and Monthly Anomalies  #########
# ---------------------------------------------------------------------------- #
#' Calculate Baseline Mean/ STD and Anomalies for NDVI, NDII, and LST
#' 
#' This function calculates the baseline mean and std for a combined MODIS and 
#' S3 raster stacks of NDVI, NDII, and LST. Then it uses these baselines to 
#' compute the monthly anomalies of the three variables.
#' 
#' @param MODIS_dir The directory containing the MODIS monthly composites.
#' @param S3_dir The directory containing the Sentinel-3 monthly composites.
#' @param out_dir The directory where the output files will be saved.
# ---------------------------------------------------------------------------- #

# Load the baseline recalculation script
step4_script <- chooseFilesCustom(
  "Choose R-Script 'MODIS_S3_Anomalies_with_Baseline_Recalculation.R'"
)
source(step4_script)

### SET DIRECTORIES
MODIS_monthlyComposites_dir <- chooseDirCustom(
  "Select folder containing the MODIS monthly NDVI, NDII, and LST composites."
)

### Execute Function
MODIS_S3_calculateBaselinesAnomalies(
  MODIS_monthlyComposites_dir, S3_calibratedMonthlyComposites_dir, 
  MODIS_S3_monthlyAnomalies_dir, MODIS_S3_Baselines_dir
)

# Clean up the environment (keep directories if needed for the next step)
to_keep4 <- c(
  "packages", "chooseDirCustom", "chooseFilesCustom",
  "step1_script", "step2_script", "step3_script", "step4_script",
  "aoi", "LST_dir", "SYN_dir", "S3_preProcessed_dir",
  "gdal_dir", "date_range", "S3_monthlyComposites_dir","S3_calibratedMonthlyComposites_dir", 
  "MODIS_S3_Baselines_dir", "MODIS_S3_monthlyAnomalies_dir",
  "MODIS_monthlyComposites_dir"
)

# Items to be removed
to_remove4 <- setdiff(ls(all.names = TRUE), to_keep4)

# Optionally print items to be removed for verification
print(paste("Removing:", to_remove4))

# Proceed with removal if verified
rm(list = to_remove4)


#==============================================================================#
##########  Step: 5. SPI Calculation for 3-monthly using TAMSAT data  ##########
# ---------------------------------------------------------------------------- #
#' Calculation of the Standardized Precipitation Index (SPI) 
#' from TAMSAT rainfall estimates.    
#' 
#' @param TAMSAT_InputFile The monthly TAMSAT rainfall data.
#' @param SPI_OutputDir The ouput directory where the SPI results will be saved
# ---------------------------------------------------------------------------- #

# Load the spi calculation script
step5_script <- chooseFilesCustom("Choose R-Script 'TAMSAT_SPI3_Calculation_Function.R'")
source(step5_script)

# Define variables for the function
TAMSAT_InputFile <- chooseFilesCustom("Select TAMSAT Monthly Rainfall File")
SPI_OutputDir <- chooseDirCustom(
  "Select output folder where the SPI results will be saved at"
)

### Execute Function
Calculate_SPI_from_TAMSAT(TAMSAT_InputFile, SPI_OutputDir)

# Clean up the environment (keep directories if needed for the next step)
to_keep5 <- c(
  "packages", "chooseDirCustom", "chooseFilesCustom",
  "step1_script", "step2_script", "step3_script", "step4_script",
  "aoi", "LST_dir", "SYN_dir", "S3_preProcessed_dir",
  "gdal_dir", "date_range", "S3_monthlyComposites_dir","S3_calibratedMonthlyComposites_dir", 
  "MODIS_S3_Baselines_dir", "MODIS_S3_monthlyAnomalies_dir",
  "MODIS_monthlyComposites_dir", 
  "TAMSAT_InputFile", "SPI_OutputDir"
)

# Items to be removed
to_remove5 <- setdiff(ls(all.names = TRUE), to_keep5)

# Optionally print items to be removed for verification
print(paste("Removing:", to_remove5))

# Proceed with removal if verified
rm(list = to_remove5)
