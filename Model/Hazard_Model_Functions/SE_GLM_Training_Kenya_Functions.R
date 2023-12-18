#==============================================================================#
#' Script for creating Training Data for the Drought Model                     #
#'                                                                             #
#' This code is used for collecting training data for the drought model.       #
#' It allows to specify a range of years and months (growing season) from which# 
#' the training data is extracted. The training data is sampled from NDVI,     #
#' NDII, and LST Anomalies as well as SPI3 precipitation data. Years           #
#' of drought/ non-drought are defined in an excel file from the FAO. For the  #
#' training data, a total of 100,000 observations are extracted for each index.#
#' This is done for drought years and non-drought years.                       #
#'                                                                             #
#' The user needs to specify:                                                  #
#' @param input_path The path to directory with the input data (NDVI/ NDII/    #
#'                   LST Anomalies, SPI3 Precipitation).                       #
#' @param drought_years_fao The path to the FAO excel file with the drought/   #
#'                          non-drought information.                           #
#' @param output_path The path to the directory to store the training data.    #
#'                                                                             #
#' @param start_date The date to start the collection of the training data.    #
#' @param end_date The date to end the collection of training data.            #
#' @param season The range of months in which to collect the training data     #
#'               (usually the growing season).                                 #
#==============================================================================#



packages = c("terra", "data.table", "openxlsx", "ncdf4")

for(i in 1:length(packages)) {
  if(packages[i] %in% rownames(installed.packages()) == FALSE) {
    install.packages(packages[i])}
}

lapply(packages, library, character.only = T)



########## Variables and Data Import #########

### PATHS:

input_dir <- "E:/Maxi/06_Drought_Model/02_inputData_Sarah/05_inputData_copernicusLC_cropHerbShrub/"
drought_years_fao <-  "E:/Maxi/06_Drought_Model/02_inputData_Sarah/Drought_Years_FAO_Residuals_Analysis.xlsx"

### VARIABLES: if you do not want to use the dialog boxes use the commented out lines

start_date <- winDialogString("Enter starting date for training period. Use format: YYYY-MM-DD","")
# start_date <- "2001-01-01"
end_date <- winDialogString("Enter end date for training period. Use format: YYYY-MM-DD","")
# end_date <- "2021-12-01"
season <- winDialogString("Enter the months for training data collection. Use format: MM, MM, ...", "")
season <- strsplit(season, ",\\s*")[[1]]
#season <- c("04", "05", "06", "07", "11", "12")



########## Functions ##########

#' Extract Years and Months
#' 
#' This function extracts the defined years and months (growing season) of data
#' from the model data (NDVI/ NDII/ LST Anomalies and SPI3) and the FAO excel file.
#' 
#' @param start_date The start date of the period in which to collect training data.
#' @param end_date The end date of the training data collection period.
#' @param season The range of months in which to collect data (growing season).
#' @param drought_years_fao The FAO excel file with information about drought and non-
#'                          drought years.
#' @param model_data The data to collect the training data from (NDVI/ NDII/ LST
#'                   Anomalies and SPI3 Precipitation).
#'                   
#' @return A list containing a list with the dates in which to collect the 
#'         training data, the filtered FAO data and the filtered model data.

extract_years_months <- function(start_date, end_date, season, sem_data, model_data){
  
  # Extract start and end year
  start_year <- as.integer(format(as.Date(start_date), "%Y"))
  end_year <- as.integer(format(as.Date(end_date), "%Y"))
  
  # Sequence of dates with a step of one month
  dates = seq(as.Date(start_date), as.Date(end_date), by = "1 month")
 
  # Create a subset of the model data for the defined month (growing season)
  for(i in 1:length(model_data)) {
    model_data[[i]] = subset(model_data[[i]], which(format.Date(dates, "%m") %in% 
                                                      season))
  }
  
  # Filter the sequence of dates to the defined months
  dates = dates[format.Date(dates, "%m") %in% season]
  
  # Extract model data for the defined years 
  for(i in 1:length(model_data)) {
    model_data[[i]] = subset(model_data[[i]], which(format.Date(dates, "%Y") %in% 
                                                      c(start_year:end_year)))
  }
  
  # Subset the FAO data to the defined years
  sem_data = subset(sem_data, sem_data$Year %in% c(start_year:end_year))
  
  return(list(dates, sem_data, model_data))
}


#' Extract Drought/ Non-Drought Index
#' 
#' This function is used to extract the data in the drought and non-drought years 
#' from the model data, which is later used to collect the training data from. 
#' The drought and non-drought years are defined in the FAO file.
#' 
#' @param sem_data The FAO excel file containing the information about drought 
#'                 and non-drought years.
#' @param model_data The data to collect the training samples from (NDVI/ NDII/
#'                   LST Anomalies and SPI3 Precipitation).
#' @param dates A list containing all months of the years in which to collect the data.

extract_drought_nonDrought <- function(sem_data, model_data, dates){
  
  # Define drought and non-drought years base on FAO data and remove undefined years
  drought_years = na.omit(sem_data$Year[sem_data$Drought == 1])
  nondrought_years = na.omit(sem_data$Year[sem_data$Drought == 0])
  
  # Create two lists: one containing the model data for the drought years for the 
  # defined period and one containing the data for the non-drought years
  drought_data = list()
  nondrought_data = list()
  
  for(i in 1:length(model_data)) {
    # Get pixel values in drought years and non-drought years 
    # -> matrix with each column representing values from 1 layer
    drought_data[[i]] = terra::values(terra::subset(model_data[[i]], 
                                                    which(format.Date(dates, "%Y") %in% 
                                                            as.character(drought_years))))
    nondrought_data[[i]] = terra::values(terra::subset(model_data[[i]], 
                                                       which(format.Date(dates, "%Y") %in% 
                                                               as.character(nondrought_years))))
    
    # Convert matrix to vector by column
    drought_data[[i]] = as.vector(drought_data[[i]])
    nondrought_data[[i]] = as.vector(nondrought_data[[i]])
  }
  
  # Convert lists with data to data frame 
  # (one column with all values of drought/non-drought per index)
  drought_data_df = do.call(cbind.data.frame, drought_data)
  nondrought_data_df = do.call(cbind.data.frame, nondrought_data)
  
  # Define and set names of the data frame
  names = c("LST", "NDII", "NDVI", "SPI3")
  names(drought_data_df) = names
  names(nondrought_data_df) = names
  
  # Convert the data frames to data tables for a faster process
  setDT(drought_data_df)
  setDT(nondrought_data_df)
  
  # Remove all rows with at least one NA 
  drought_data_df = na.omit(drought_data_df)
  nondrought_data_df = na.omit(nondrought_data_df)
  
  
  return(list(drought_data_df, nondrought_data_df))
}


#' Sample Training Data
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
#'                 non-drought years.

sample_training_data <- function(start_date, end_date, season, input_dir, drought_years_fao){
  
  ########## Data Import ##########
  
  ### FAO DATA: without sugar cane
  no_sugar_cane = read.xlsx(drought_years_fao, sheet = "No Sugar cane", startRow = 2)

  # remove not needed columns: Area, Value, Mean and Stddev
  sem_data = no_sugar_cane[c("Year", "Drought")]

  ### MODEL INPUT DATA: NDVI, NDII, LST Anomalies and SPI3 Precipitation
  model_files = list.files(input_dir, pattern = ".nc$", recursive = T, full.names = T)

  model_data <- list()

  for (i in 1:length(model_files)) {
    model_data[i] <- rast(model_files[i])
  }
  
  ########## Extract Years and Months ##########
  
  extracted_data <- extract_years_months(start_date, end_date, season, sem_data, model_data)
  
  dates <- extracted_data[[1]]
  sem_data <- extracted_data[[2]]
  model_data <- extracted_data[[3]]
  
  ########## Extract Drought/ Non-Drought Index ##########
  extracted_drought_nonDrought <- extract_drought_nonDrought(sem_data, model_data, dates)
  
  drought_data_df <- extracted_drought_nonDrought[[1]]
  nondrought_data_df <- extracted_drought_nonDrought[[2]]
  
  ########## Training Data Sampling ##########
  
  # Extract random sample of 100k observations per index (each row equals 1 pixel)
  drought_data_df <- drought_data_df[sample(.N, 100000)]
  nondrought_data_df <- nondrought_data_df[sample(.N, 100000)]
  
  ########## Post-Processing ##########
  
  # Convert data tables back to data frame
  drought_data_df = as.data.frame(drought_data_df)
  nondrought_data_df = as.data.frame(nondrought_data_df)
  
  # Check if there are NA values left in data frames and add column with drought
  # and non-drought indication (drought = 1, non-drought = 0)
  has_na_drought <- any(is.na(drought_data_df))
  if(has_na_drought) {
    print("The data frame contains NA values.") } else {
      for(i in 1:nrow(drought_data_df)) {
        drought_data_df$Drought[i] = 1 
      }
    }
  
  has_na_nondrought <- any(is.na(nondrought_data_df))
  if(has_na_nondrought) {
    print("The data frame contains NA values.") } else {
      for(i in 1:nrow(nondrought_data_df)) {
        nondrought_data_df$Drought[i] = 0 
      }
    } 
  
  # Create data frame with all training data (combine drought and non-drought data frame)
  training_data = rbind(drought_data_df, nondrought_data_df)
  
  ########## Write Results ##########
  output_name <- "Training_Data_Drought_Model_Kenya_TEST.csv"
  
  write.table(training_data, paste0(input_dir, output_name), 
              col.names = T, row.names = F, sep = ",")
}
sample_training_data(start_date, end_date, season, input_dir, drought_years_fao)
