#==============================================================================#
#' Script with Functions for Resampling the Model data, Collecting Training    #
#' Data and Executing the Drought Model                                        #
#'                                                                             #
#' This code is used to resample the model input data (NDVI, NDII, and LST     #
#' Anomalies + SPI3 Precipitation) to the same spatial resolution (~0.01       #
#' degree of latitude or ~1.11km) and mask the data with a Copernicus Land     #
#' Cover Mask. Drought and non-drought years are identified with FAO crop      #
#' yield data. In the third step, a total of 100,000 training samples per      #
#' index are collected for specified range of years and months. These training #
#' samples of the different variables are then tested for autocorrelation.     #
#' Afterwards, a logistic regression model (GLM) is build and trained based on #
#' the training data. Then this model is used to predicted the drought         #
#' probability for a specified number of years and months.                     #
#'                                                                             #
#' Once the modelling part is finished, the percentage of irrigation is        #
#' calculated from the farming systems classification results. The resulting   #
#' map is then used as input alongside with other data (GDP, population density#
#' livestock density, relative wealth index) to calculate the drought          #
#' vulnerability index. In the last step the drought probability maps alongside#
#' the results of the DVI calculation are used to determine the drought risk   #
#' for the predicted months and years as well as for the growing season.       #
#'                                                                             #
#' If necessary, the user can exclude variables from the statistical           #
#' correlation test, the modeling, and prediction. This is done by choosing the#
#' variable to be excluded in the dialog box that is popping up when executing #
#' the desired function.                                                       #
#'                                                                             #
#' 1. Resampling Parameters:                                                   #
#' @param anomalies_dir The path to directory containing the NDVI, NDII, and   #
#'                      LST Anomalies.                                         #           
#' @param spi3_dir The path to the directory containing the TAMSAT SPI3        #
#'                 Precipitation data.                                         #
#' @param lc_file The Copernicus Land Cover file (land cover mask).            #
#' @param output_dir The directory to store the resampled and masked anomalies #
#'                   and SPI3 data.                                            #
#'                                                                             #
#' 2. Identify drought years with yield data:                                  #
#' FAO crop yield data needs to be downloaded as .csv file from:               #
#' https://www.fao.org/faostat/en/#data/QCL                                    #
#'    -"country": select your country of interest                              #
#'    -"elements": select "yield"                                              #
#'    -"items": select four to five main food crops of the country             #
#'    -"years": select years 1990 to today                                     #                                          
#' @param crop_data_path Path to downloaded FAO crop yield data csv-file.      #
#' @param drought_years_fao Path to save the output file containing the        #
#'                          information about drought and non-drought years.   #
#'                                                                             # 
#' 3. Training Parameters:                                                     #
#' @param start_date The start date of the years to collect the training data. #
#' @param end_date The end date for the collection of the training data.       #
#' @param season A list of months to collect the training data in              #
#'               (growing season).                                             #
#'                                                                             #
#' @param input_dir The path to the model data (NDVI/ NDII/ LST Anomalies and  #
#'                  SPI3 Precipitation Data) and (training data).              #
#'                                                                             #
#'                                                                             #
#' 4. - 6. Correlation Test, Model Building, Predicting                        #
#' @param predicted_years The years for which the drought probability will be  #
#'                        predicted.                                           #
#' @param predicted_months The months for which the drought probability will   #
#'                         be predicted.                                       #
#'                                                                             #
#' @param model_dir The path to the folder where the drought probability maps  #
#'                  will be stored.                                            #  
#'                                                                             #
#' 7. Percentage of Irrigation                                                 #
#' @param farmingSystems_file The path to the file with the results of the     #
#'                            farming systems classification.                  #
#' @param aoi The path to the shapefile of Kenya.                              #
#' @param output_dir The path to the directory where the map of the percentage #
#'                   of irrigation, the vulnerability map and the drought risk #
#'                   results are stored.                                       #
#'                                                                             #
#' 8. Drought Vulnerability Index (DVI)                                        #
#' @param dviInput_dir The path to the directory containing the input data for #
#'                     the DVI calculation (GDP, population and livestock      #
#'                     density, percentage of irrigation, and relative wealth  #
#'                     index (RWI)).                                           #
#' @param lc_dir The path to the directory containing the Copernicus Land Cover# 
#'               classification as well as different land cover masks.         #
#' @param normData_dir The path to the directory where the normalized input    #
#'                     data will be stored.                                    #
#'                                                                             #
#' 9. Drought Risk                                                             #
#' @param dvi_file The path to the drought vulnerability index file.           #
#==============================================================================#


#==============================================================================#
################# Install and load all required packages #######################
packages = c("furrr", "terra", "R.utils", "Rcpp", "corrplot", "pscl", "R.utils",
             "utils", "strucchange", "segmented", "data.table", "dplyr")

for (i in 1:length(packages)) {
  if(packages[i] %in% rownames(installed.packages()) == FALSE) {
    install.packages(packages[i])}
}

lapply(packages, library, character.only = T)

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
################## Step 1: Resample and Mask Data ##############################
# ---------------------------------------------------------------------------- #
#' Resample and Mask Data to be ready for model input
#' 
#' This function is used to resamples the model input data (NDVI, NDII, and LST 
#' Anomalies + SPI3 Precipitation data) to the same spatial resolution 
#' (~ 0.01 degree of latitude) and masks the resampled data with a Copernicus 
#' Land Cover mask.
#' 
#' @param anomalies_dir The path to directory containing the NDVI, NDII, and LST
#'                      Anomalies.       
#' @param spi3_dir The path to the directory containing the Tamsat SPI3 Precipitation.
#' @param lc_file The Copernicus Land Cover file (land cover mask).  
#' @param input_dir The directory to store the resampled and masked anomalies 
#'                   and SPI3 data.      
# ---------------------------------------------------------------------------- #

### Load the step 1 script
step1_script <- chooseFilesCustom("Choose R-Script 'Resampling_Functions.R'")
source(step1_script)

#### SET DIRECTORIES
anomalies_dir <- chooseDirCustom(caption = "Select Folder containing the NDVI, NDII, and LST Anomalies")
spi3_dir <- chooseDirCustom(caption = "Select Folder containing the SPI3 Precipitation Data")
lc_file <- chooseFilesCustom(caption = "Select Copernicus Land Cover File (Land Cover Mask)")
output_dir <- chooseDirCustom(caption = "Select Folder to store the Model Input Data")

#### RUN FUNCTION
resample_LCmask(anomalies_dir = anomalies_dir, spi3_dir = spi3_dir, 
                            lc_file = lc_file, output_dir = output_dir)

#==============================================================================#
###############  Step 2: Identification of Drought years  ######################
# ---------------------------------------------------------------------------- #
#' Identify drought years with yearly FAO yield data of main Food Crops:
#'  
#' This function is used to identify drought years with yield data from the FAO.
#' It outputs a table indicating drought (1) or no drought (0) for every year
#' being used as an input for the training function in Step 3. 
#' 
#' !!! TO DO !!!
#' Since the function identifies years with low yield as drought years and does 
#' not consider flooding events or pests and diseases, please do a quick internet 
#' research on drought events to double-check the classification after running 
#' the function. Corrections can be made by changing the numbers 0 (no drought) 
#' and 1 (drought) in the csv-file.
#' 
#' @param crop_data_path Path to downloaded FAO crop yield data csv-file.      
#' @param drought_years_fao Path to save the output file containing the        
#'                          information about drought and non-drought years.  
# ---------------------------------------------------------------------------- #

### Load the step 2 script
step2_script <- chooseFilesCustom("Choose R-Script 'Crop_Yield_Drought_Function.R'")
source(step2_script)


#### SET DIRECTORIES
crop_data_path <- chooseFilesCustom(caption = "Select FAO crop yield file")
drought_years_fao <- chooseDirCustom(caption = "Select path to save the output file 
                                containing the information about drought and 
                                non-drought years")


#### RUN FUNCTION
get_drought_years_from_yield_data(crop_data_path = crop_data_path, 
                                  drought_years_fao = drought_years_fao)

#==============================================================================#
##########################  Step 3: Training Data ##############################
# ---------------------------------------------------------------------------- #  
#' 3. Sample Training Data based on identified drought years:
#' 
#' This function is used to sample the training data from the NDVI, NDII, LST 
#' Anomalies and the SPI3 precipitation data for drought and non-drought years
#' within a specified period of months (growing season) and years.
#' 
#' @param start_date The start date of the years to collect data in.
#' @param end_date The end date of the years to collect data in.
#' @param season A list of the months to collect data in (growing season).
#' @param input_dir The path to the model data (NDVI, NDII, LST Anomalies and SPI3).
#' @param drought_years_fao The FAO excel file with information about the drought and 
#'                          non-drought years.
# ---------------------------------------------------------------------------- #

### Load the step 3 script
step3_script <- chooseFilesCustom("Choose R-Script 'SE_GLM_Training_Kenya_Functions.R'")
source(step3_script)

#### Directories (not all defined in step 1 and 2 already)
raster_dir <- chooseDirCustom(caption = "Select Folder with Model Input Data")
drought_years_fao_file <- chooseFilesCustom("Choose File with Drought Years 
                                            from Previos Step")

#### SET VARIABLES
start_date <- winDialogString("Enter starting date for training period. 
                              Use format: YYYY-MM-DD","") #2001-01-01 (fixed)
end_date <- winDialogString("Enter end date for training period. 
                            Use format: YYYY-MM-DD","") #2021-12-31
season <- winDialogString("Enter the months for training data collection. 
                          Use format: MM, MM, ...", "") 
#season <- c("04", "05", "06", "07", "11", "12")
season <- strsplit(season, ",\\s*")[[1]]

#### RUN FUNCTION
sample_training_data(start_date = start_date, end_date = end_date, season = season, 
                     input_dir = raster_dir, drought_years_fao = drought_years_fao_file)

#==============================================================================#
############################  Step 4: Correlation  #############################
# ---------------------------------------------------------------------------- #
#' 4. Test Statistical Correlation between model variables
#' 
#' This function checks for autocorrelation between the model input variables
#' (NDVI, NDII, and LST Anomalies + SPI3 Precipitation data) in the training data. 
#' 
#' @param input_dir The path to directory with the training data.
# ---------------------------------------------------------------------------- #

### Load the step 4 script
step4_script <- chooseFilesCustom("Choose R-Script 'SE_GLM_Functions.R'")
source(step4_script)

#### Directories
training_data_path = chooseDirCustom("Choose Directory that contains the
                                     Training data for the Model")

#### RUN FUNCTION
test_correlation(input_dir = training_data_path)

#==============================================================================#
##########################  Step 5: Modelling  #################################
# ---------------------------------------------------------------------------- #
#' 5. Build and Evaluate Logistic Model (GLM)
#' 
#' This function builds a general linear model (GLM) and trains it with the 
#' training data. Additionally, the model is statistically evaluated based on
#' the z-value, P, Std. Error, and McFadden's pseudo R^2.
#' 
#' @param input_dir The path to the directory containing the training data.
# ---------------------------------------------------------------------------- #

### Load the step 5 script: Already loaded in Step 4

#### Directory set in step 4.

#### RUN FUNCTION
model_training_evaluation(input_dir = training_data_path)

#==============================================================================#
##########################  Step 6: Drought Hazard  ############################
# ---------------------------------------------------------------------------- #
#' 6. Create Drought Probability Maps with the developed model
#' 
#' This function is used to apply the model to predict the drought probability
#' for defined years and months. The output are drought probability maps for the
#' specified months.
#' 
#' @param start_date The start date of the period in which the training data 
#'                   was collected.
#' @param end_date The end date of the training data collection period.
#' @param predicted_years The years for which the drought probability will be predicted.
#' @param predicted_months The months for which the drought probability will be predicted.
#' @param input_dir The path to directory with the model data (NDVI, NDII, and
#'                  LST Anomalies + SPI3 Precipitation Data) and training data.
#' @param model_dir The path to the folder where the drought probability maps
#'                   will be stored.
# ---------------------------------------------------------------------------- #

### Load the step 6 script: Already loaded in Step 4

#### SET DIRECTORIES
model_dir <- chooseDirCustom(caption = "Select Folder to store the Drought Probability Maps")
input_data_dir <- chooseDirCustom("Select Folder with Model Input Data")

#### SET VARIABLES
start_date <- winDialogString("Enter starting date for Application period. 
                              Use format: YYYY-MM-DD","") #2001-01-01 -> has to be starting date of dataset
end_date <- winDialogString("Enter end date for training period. 
                            Use format: YYYY-MM-DD","") #2023-12-31
predicted_months <- winDialogString("Enter the months for predicting the drought probability. Use format: MM, MM, ...", "")
predicted_months <- strsplit(predicted_months, ",\\s*")[[1]]
predicted_years <- winDialogString("Enter the years for predicting the drought probability. Use format: YYYY, YYYY, ...", "")
predicted_years <- strsplit(predicted_years, ",\\s*")[[1]]

# predicted_months <- c("01","02","03","04","05","06", "07", "08", "09", "10","11","12")
# predicted_years <- c("2018", "2022")

#### RUN FUNCTION
drought_probability_maps(start_date = start_date, end_date = end_date, predicted_years = predicted_years, 
                         predicted_months = predicted_months, input_dir = input_data_dir, output_dir = model_dir)

#==============================================================================#
####################  Step 7: Process Irrigation Data  #########################
# ---------------------------------------------------------------------------- #
#' 7. Calculate Percentage of Irrigation
#' 
#' This function is used to calculate the percentage of irrigated areas based on
#' the farming systems classification by computing the sum of pixels that are 
#' irrigated within a moving window of the size 100 x 100 pixels and dividing this
#' value by the total number of pixels within this window (10,000). Additionally,
#' the raster is resampled to the spatial resolution of the drought model output
#' (0.0125 degree or ~1km), clipped to the extent of Kenya and saved.
#' 
#' @param farmingSystems_file The path to the file with the results of the 
#'                            farming systems classification.
#' @param model_dir The directory with the results of the drought model.
#' @param aoi The path to the shapefile of Kenya.
#' @param irrig_perc_dir The path to the directory where the map of the      
#'                         percentage of irrigation is stored.  
# ---------------------------------------------------------------------------- #

### Load the step 7 script
step7_script <- chooseFilesCustom("Choose R-Script 'DroughtVulnerabilityIndex_Functions.R'")
source(step7_script)

#### SET DIRECTORIES
farmingSystems_file <- chooseFilesCustom(caption = "Select Farming Systems Classification File")
aoi <- chooseFilesCustom(caption = "Select Shapefile of Kenya")
irrig_perc_dir <- chooseDirCustom(caption = "Select Folder with the Subfolders to store the Percentage 
                                              of Irrigation and the Drought Vulnerability Data")

#### RUN FUNCTION
percentage_irrigation(farmingSystems_file = farmingSystems_file, model_dir = model_dir, 
                      aoi = aoi, output_dir = irrig_perc_dir)

#==============================================================================#
###################  Step 8: Drought Vulnerability  ############################
# ---------------------------------------------------------------------------- #
#' 8. Calculate Drought Vulnerability Index (DVI)
#' 
#' This function calculates the drought vulnerability index from different indicators
#' (GDP, population density, livestock density, percent irrigation, and relative 
#' wealth index(RWI)). Before the DVI can be calculated,  the population density is
#' capped to 300 counts per km^2 which referrers to urban areas (which are excluded).
#' Then it calculates the mean livestock density from the sheep, goat and cattle 
#' density. Additionally, it resamples the vulnerability data, the Copernicus Land 
#' Cover raster, and the land cover mask (cropland, herbs, shrubs) to the spatial 
#' resolution of the drought model output (0.0125 degree or ~1.113 km). After the 
#' resampling step, it uses the Copernicus Land Cover to create a mask for urban 
#' areas and applies it to the vulnerability data. Then, it normalizes the input 
#' data and inverts the values of GDP, percentage of irrigation, and relative wealth 
#' index. Afterwards, the DVI is calculated by computing the sum of the vulnerability 
#' layers and dividing it by the number of input layers. The result is masked with
#' the crop/ herb/ shrub mask and then saved.
#' 
#' @param dviInput_dir The path to the directory containing the input data for the
#'                    DVI calculation (GDP, population and livestock density, 
#'                    percentage of irrigation, and relative wealth index (RWI)).
#' @param lc_dir The path to the directory containing the Copernicus Land Cover 
#'               classification as well as different land cover masks.
#' @param model_dir The path to the directory with the results of the drought model.
#' @param normData_dir The path to the directory where the normalized input data 
#'                     will be stored, which are an intermediate output.
#' @param dvi_dir The path to the directory where the output should be saved.
# ---------------------------------------------------------------------------- #

### Load the step 8 script: Already loaded in Step 7

#### SET DIRECTORIES
gdp <- chooseFilesCustom(caption = "Select the Copernicus Land Cover Data")
popdens <- chooseFilesCustom(caption = "Select the Copernicus Land Cover Data")
rwi <- chooseFilesCustom(caption = "Select the Copernicus Land Cover Data")
livestock <- chooseDirCustom(caption = "Select the Copernicus Land Cover Data")
percprec <- chooseFilesCustom(caption = "Select the Copernicus Land Cover Data")
lc_dir <- chooseFilesCustom(caption = "Select the Copernicus Land Cover Data")
lcmask_dir <- chooseFilesCustom(caption = "Select the Copernicus Land Cover Data")
normData_dir <- chooseDirCustom(caption = "Select Folder to store the normalized data created in this function")
dvi_dir <- chooseDirCustom(caption = "Select Folder where the output should be saved")

#### RUN FUNCTION
drought_vulnerability_index(gdp = gdp, popdens = popdens, rwi = rwi, 
                            livestock = livestock, percprec = percprec,lc_dir = lc_dir, 
                            lcmask_dir = lcmask_dir, model_dir = model_dir, 
                            normData_dir = normData_dir, output_dir = dvi_dir)

#==============================================================================#
##########################  Step 9: Drought Risk ###############################
# ---------------------------------------------------------------------------- #
#' 9. Calculate the Drought Risk 
#' 
#' This function calculates the drought risk for the predicted months based on 
#' the results of the drought model and the Drought Vulnerability Index (DVI). 
#' The drought risk is calculated for every month of the predicted years. 
#' Additionally, the mean drought risk for every year and growing season is computed.                                                                #
#'                                                                             
#' @param model_dir The path to the directory with the results from the drought
#'                  model (drought probability maps).                          
#' @param dvi_file The path to the drought vulnerability index file.           
#' @param output_dir The path to the directory where the output files of the   
#'                   drought risk calculation will be stored. 
# ---------------------------------------------------------------------------- #

### Load the step 9 script
step9_script <- chooseFilesCustom("Choose R-Script 'DroughtRisk_Functions.R'")
source(step9_script)

####SET DIRECTORIES
hazard_dir <- chooseDirCustom(caption = "Select Dir with Hazard Data")
dvi_file <- chooseFilesCustom(caption = "Select File with DVI Data")
output_dir <- chooseDirCustom(caption = "Select Folder to save the drought risk files")

#### RUN FUNCTION
drought_risk_calc(hazard_dir = hazard_dir, dvi_file = dvi_file, output_dir = output_dir)
