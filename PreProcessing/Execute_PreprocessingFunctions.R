#==============================================================================#
#' Script with functions for pre-processing the Sentinel-3 data and            #
#' calculation of monthly anomalies. Additionally, it can be used to re-       #
#' calculate the MODIS and Sentinel-3 baselines and create new monthly         #
#' anomalies.                                                                  #
#'                                                                             #
#' This code includes functions to pre-process the raw Sentinel-3 data and     #
#' create monthly composites of NDVI, NDII, and LST. Additionally, the         #
#' composites are calibrated to match the MODIS composites. Then, the script   #
#' can be used to compute monthly anomalies of the three variables.            #
#' Optionally, there is a function to re-calculate the baselines every few     #
#' years to keep them up-to-date.If the baselines and anomalies should be re-  #
#' calculated, the user needs to comment out the function                      #
#' "S3_calculateAnomalies" and uncomment the function                          # 
#' "MODIS_S3_calculateBaselinesAnomalies".                                     #
#'                                                                             #
#' The user needs to specify:                                                  #
#'                                                                             #
#' 1. Pre-Processing Parameters:                                               #   
#' @param LST_dir The path to the raw Sentinel-3 SL_2_LST products.            #
#' @param SYN_dir The path to the raw Sentinel-3 SY_2_SYN products.            #
#' @param aoi The path to the shapefile of Kenya.                              #
#' @param S3_preProcessed_dir The path to the directory containing the pre-    #
#'                            processed S3 LST and SYN data.                   #
#'                                                                             #
#' 2. Composite creation and Calibration Parameters:                           #
#' @param gdal_dir The path to the GDAL installation directory.                #
#' @param S3_ndvi_dir The directory containing the pre-processes S3 NDII files.#
#' @param S3_ndii_dir The directory containing the pre-processes S3 NDVI files.#
#' @param S3_lst_dir The directory containing the pre-processes S3 LST files.  #
#' @param date_range The date range for which the monthly composites are       #
#'                   created.                                                  #
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
#==============================================================================#


# Install and load required packages
packages <- c("raster", "terra", "doParallel", "foreach", "tools", "furrr", "R.utils", "Rcpp")

for (i in 1:length(packages)) {
  if(packages[i] %in% rownames(installed.packages()) == FALSE) {
    install.packages(packages[i])}
}

lapply(packages, library, character.only = T)

# Set working directory
shell.exec("explorer.exe")
wd <- choose.dir(caption = "Select working directory containing all scripts and data to run this code")
setwd(wd)


# Load Scripts:
# List of scripts to access functions
S3_PreProcessing <- choose.files(getwd(), caption = "Choose R-Script 'S3_PreProcessing.R'")
S3_Monthly_Composites_MODIS_Calibration <- choose.files(getwd(), caption = "Choose R-Script 'S3_Monthly_Composites_MODIS_Calibration.R'")
S3_AnomaliesCalculation_no_Baseline_Recalculation <- choose.files(getwd(), caption = "S3_AnomaliesCalculation_no_Baseline_Recalculation.R'")
MODIS_S3_Anomalies_with_Baseline_Recalculation <- choose.files(getwd(), caption = "Choose R-Script 'MODIS_S3_Anomalies_with_Baseline_Recalculation.R'")

script_list <- c(S3_PreProcessing, S3_Monthly_Composites_MODIS_Calibration, 
                 S3_AnomaliesCalculation_no_Baseline_Recalculation, MODIS_S3_Anomalies_with_Baseline_Recalculation)

lapply(script_list, source)


########## Step 1: Pre-processing ##########

#' Pre-Process Sentinel-3 SL_2_LST and SY_2_SYN products
#' 
#' This function is used to extract the downloaded SL_2_LST/ SY_2_SYN zip files 
#' to obtain the desired LST, NDVI, and NDII layers, crop the layers to the extent 
#' of Kenya and mask out all cloud covered pixels. The output is then stored in
#' the respective NDVI, LST, and NDII folder.
#' 
#' @param LST_directory_path The directory containing the zipped or unzipped S3
#'                           SL_2_LST products.
#' @param SYN_directory_path The directory containing zipped or unzipped 
#'                            sentinel-3 SY_2_SYN products.
#' @param aoi The path to the shapefile of Kenya to which the data is clipped.
#' @param output_dir The path to a folder where the processed products will be stored.

#### SET DIRECTORIES
LST_dir <- choose.dir(getwd(), caption = "Select Folder with raw Sentinel-3 SL_2_LST products.")
SYN_dir <- choose.dir(getwd(), caption = "Select Folder with raw Sentinel-3 SY_2_SYN products.")
aoi <- choose.files(getwd(), caption = "Select shapefile with AOI of Kenya.")
S3_preProcessed_dir <- choose.dir(getwd(), caption = "Select Folder to store the pre-processed Sentinel-3 files.")

#### RUN FUNCTIONS
process_S3_SL_2_LST_products(LST_directory_path = LST_dir, aoi = aoi, output_dir = S3_preProcessed_dir)

process_S3_SY_2_SYN_products(SYN_directory_path = SYN_dir, aoi = aoi, output_dir = S3_preProcessed_dir)




########## Step 2: Composite Creation and Calibration ##########

#' Create Monthly NDVI, NDII, and LST Composites
#' 
#' These functions creates monthly composites of the NDVI, NDII, and LST from three 
#' lists of input files.
#' 
#' @param gdal_path The path to the GDAL installation directory.
#' @param lst_dir The directory containing the LST files.
#' @param ndvi_dir The directory containing the NDVI files.
#' @param ndii_dir The directory containing the NDII files.
#' @param date_range  A character vector containing the start and end dates of the 
#'                    desired range in the format "YYYY-MM-DD".
#' @param output_dir The output directory where the composite files will be saved.
#' 
#' @example 
#' createMonthlyComposites_ndvi_lst(gdal_path, lst_dir, ndvi_dir, date_range, output_dir)
#' createMonthlyComposites_ndii(gdal_path, ndii_dir, date_range, output_dir)

### SET DIRECTORIES
gdal_dir <- choose.dir(getwd(), caption = "Select Folder of the QGIS bin.")

S3_ndvi_dir <- choose.dir(getwd(), caption = "Select Folder containing the pre-processed S3 NDVI files created in Step 1.")
S3_ndii_dir <- choose.dir(getwd(), caption = "Select Folder containing the pre-processed S3 NDII files created in Step 1.")
S3_lst_dir <- choose.dir(getwd(), caption = "Select Folder containing the pre-processed S3 LST files created in Step 1.")

S3_monthlyComposites_dir <- choose.dir(getwd(), caption = "Select Folder to store the monthly S3 NDVI, NDII, and LST composites")


### SET VARIABLES
date_range <- winDialogString("Enter the date range for which the monthly composites are created. Use format: YYYY-MM-DD, YYYY-MM-DD", "")

#### RUN FUNCTIONS
createMonthlyComposites_ndvi_lst(gdal_dir = gdal_dir, lst_dir = S3_lst_dir, ndvi_dir = S3_ndvi_dir, 
                                 date_range = date_range, output_dir = S3_monthlyComposites_dir)

createMonthlyComposites_ndii(gdal_dir = gdal_dir, ndii_dir = S3_ndii_dir, date_range = date_range, 
                             output_dir = S3_monthlyComposites_dir)


#' Calibrate Monthly NDVI, NDII, and LST Composites
#' 
#' This function calibrates the monthly composites of the NDVI, NDII, and LST 
#' from three lists of input files.
#' 
#' @param composites_dir The directory containing the monthly LST and NDVI composites.
#' @param calibratedOutput_dir The output directory where the calibrated composite 
#'                             files will be saved.
#'
#' @example 
#' calibrateMonthlyComposites_ndvi_lst(composites_dir, calibratedOutput_dir)
#' calibrateMonthlyComposites_ndii(composites_dir, calibratedOutput_dir)

#### SET DIRECTORIES
S3_calibratedMonthlyComposites <- choose.dir(getwd(), caption = 
"Select Folder to store the calibrated monthly composites (NDVI, NDII, and LST).")

#### RUN FUNCTIONS
calibrateMonthlyComposites_ndvi_lst(composites_dir = S3_monthlyComposites_dir, 
                                    calibratedOutput_dir = S3_calibratedMonthlyComposites_dir)
calibrateMonthlyComposites_ndii(composites_dir = S3_monthlyComposites_dir, 
                                calibratedOutput_dir = S3_calibratedMonthlyComposites_dir)



########## Step 3: Monthly Anomalies Calculation ##########

#' Calculate monthly Anomalies for NDVI, NDII, and LST
#' 
#' This function calculates baseline mean and std to compute the monthly 
#' anomalies of ndvi, ndii, and lst and saves the results
#' in the defined output folder.
#' @param baselines_dir The directory containing the baseline mean and std files.
#' @param S3_dir The directory containing the Sentinel-3 monthly composites.
#' @param out_dir The directory where the output files will be saved.
#' 
#' @examples
#' S3_calculateAnomalies(baselines_dir, S3_dir, out_dir)

#### SET DIRECTORIES
MODIS_S3_Baselines_dir <- choose.dir(getwd(), caption = "Select Folder containing the mean and std baseline files.")
MODIS_S3_monthlyAnomalies_dir <- choose.dir(getwd(), caption = "Select Folder to store the monthly anomalies.")

temp_dir <- choose.dir(getwd(), caption = "Select temporary directory.")

#### RUN FUNCTION
S3_calculateAnomalies(baselines_dir = MODIS_S3_Baselines_dir, 
                      S3_dir = S3_calibratedMonthlyComposites_dir, 
                      out_dir = MODIS_S3_monthlyAnomalies_dir)



########## Optional: 4. Re-Calculation of Baselines and Monthly Anomalies ##########

#' Calculate Baseline Mean/ STD and Anomalies for NDVI, NDII, and LST
#' 
#' This function calculates the baseline mean and std for a combined MODIS and 
#' S3 raster stacks of NDVI, NDII, and LST. Then it uses these baselines to 
#' compute the monthly anomalies of the three variables.
#' 
#' @param MODIS_dir The directory containing the MODIS monthly composites.
#' @param S3_dir The directory containing the Sentinel-3 monthly composites.
#' @param out_dir The directory where the output files will be saved.
#' 
#' @examples
#' calculate_baselines_anomalies(MODIS_monthlyComposites_path, 
#' S3_monthlyCalibratedComposites_path, out_dir)

### SET DIRECTORIES
MODIS_monthlyComposites_dir <- choose.dir(getwd(), caption = "Select folder containing the MODIS monthly 
                                                              NDVI, NDII, and LST composites.")

# MODIS_S3_calculateBaselinesAnomalies(MODIS_dir = MODIS_monthlyComposites_dir, 
#                                      S3_dir = S3_calibratedMonthlyComposites_dir, 
#                                      out_dir = MODIS_S3_Baselines_dir) 