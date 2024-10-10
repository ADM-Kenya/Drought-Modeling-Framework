#==============================================================================#
#' Script to run the Drought Model (GLM SEM) for Kenya                         #
#'
#' This code is used to create drought probability maps for Kenya based on a   #
#' logistic regression model.                                                  #
#' First, the statistical correlation between the input variables of the model #
#' (NDVI, NDII, LST Anomalies + SPI3 Precipitation Data) is tested. Then the   #
#' model is build and statistically evaluated. The output are drought          #
#' probability maps for the defined months and years.                          #
#'                                                                             #
#' The user needs to specify:                                                  #
#' @param input_dir The path to directory with the model data (NDVI, NDII, and #
#'                  LST Anomalies + SPI3 Precipitation Data) and training data.#
#' @param output_dir The path to the folder where the drought probability maps #
#'                   will be stored.                                           #
#' @param start_date The start date of the period in which the training data   #
#'                   was collected.                                            #
#' @param end_date The end date of the training data collection period.        #
#' @param predicted_years The years for which the drought probability will be  #
#'                        predicted.                                           #
#' @param predicted_months The months for which the drought probability will   #
#'                         be predicted.                                       #
#==============================================================================#



# Install and load required packages
#packages <-  c("corrplot", "terra", "pscl")

#for(i in 1:length(packages)){
#  if(packages[i] %in% rownames(installed.packages()) == FALSE) {
#    install.packages(packages[i])}
#}

#lapply(packages, library, character.only = T)



########## Paths and Variables ##########

### PATHS:
#input_dir <- "E:/Maxi/06_Drought_Model/02_inputData_Sarah/05_inputData_copernicusLC_cropHerbShrub"
#output_dir <- "E:/Maxi/06_Drought_Model/TEST/"


### VARIABLES: if you do not want to use the dialog boxes use the commented out lines
#start_date <- winDialogString("Enter starting date for training period. Use format: YYYY-MM-DD","")
# start_date <- "2001-01-01"
#end_date <- winDialogString("Enter end date for training period. Use format: YYYY-MM-DD","")
# end_date <- "2022-12-01"

#predicted_months <- winDialogString("Enter the months for predicting the drought probability. Use format: MM, MM, ...", "")
#predicted_months <- strsplit(predicted_months, ",\\s*")[[1]]
# predicted_months <- c("01","02","03","04","05","06", "07", "08", "09", "10","11","12")
#predicted_years <- winDialogString("Enter the years for predicting the drought probability. Use format: YYYY, YYYY, ...", "")
#predicted_years <- strsplit(predicted_years, ",\\s*")[[1]]
# predicted_years <- c("2018", "2022")



########## Functions ##########

#' #' Remove Layer from Raster Stack
#' #'
#' 
#' remove_layer <- function(stack, layer_name) {
#'   return(subset(stack, select = which(names(stack) != layer_name)))
#' }


#' Test Statistical Correlation
#' 
#' This function checks for autocorrelation between the model input variables
#' (NDVI, NDII, and LST Anomalies + SPI3 Precipitation data). 
#' 
#' @param input_dir The path to directory with the training data.

test_correlation <- function(input_dir){
  
  ########## Import Training Data ##########
  
  drought_data <- read.table(list.files(path = input_dir, pattern = ".csv$", recursive = T, 
                                        full.names = T), header = T, sep = ",")
  
  ########## Select Variables ##########
  
  variable_choice <- menu(c("Include all variables" ,"Exclude NDVI", "Exclude NDII", 
                            "Exclude LST", "Exclude SPI3"), 
                          "Choose Variables of Correlation Test", graphics = TRUE)
  
  if (variable_choice == 2){
    drought_data <- drought_data[, -which(names(drought_data) == "NDVI"), drop = FALSE]
    variable_choice <- 1
  }
  if (variable_choice == 3){
    drought_data <- drought_data[, -which(names(drought_data) == "NDII"), drop = FALSE]
    variable_choice <- 1
  }
  if (variable_choice == 4){
    drought_data <- drought_data[, -which(names(drought_data) == "LST"), drop = FALSE]
    variable_choice <- 1
  }
  if (variable_choice == 5){
    drought_data <- drought_data[, -which(names(drought_data) == "SPI3"), drop = FALSE]
    variable_choice <- 1
  }
  
  if (variable_choice ==1){
  
    ########## Test for Autocorrelation ##########
  
    # Explanation: pairwise correlation = absolute values of correlation coefficients
    autocorr_df <- drought_data[, c(1: (ncol(drought_data)-1))]
    cor_results <- round(cor(autocorr_df, use = "complete.obs"), 2)
    
    # Define if Corrplot should be displayed
    corrplot_choice <- menu(c("Yes" ,"No"), "Show Correlation Plot", 
                            graphics = TRUE)
    
    if (corrplot_choice == 1){
      
      # Visualize correlation results (should be <0.7)
      corrplot(cor_results, method = "color", type = "lower", order = "hclust", 
               tl.col = "black", tl.srt = 45, addCoef.col = "black", diag = F)
      mtext("Correlation Plot (Training Data)", side = 3, line = 2.5, 
            at =  2.25, cex = 1.5, col = "black")
      corrplot_choice <- 2
    }
    
    if (corrplot_choice == 2){
      
      # Calculate the determinate of the correlation matrix
      # Explanation: Explanation: 1 = no colinearity, 0 = high colinearity  (0.43)
      det <- det(cor_results) 
  
      print(paste('Determinant Correlation Matrix (1 = no colinearity, 0 = high colinearity): ', det))
  
      # Compute the condition number of the correlation matrix
      # Explanation: overall sum of multi-colinearity = highest condition index (should be <30)
      kappa <- kappa(cor_results)
  
      print(paste('Condition Number of the Correlation Matrix (should be <30): ', kappa))
    }
  }
}


#' Logistic Model (GLM)
#' 
#' This functions builds a general linear model (GLM) and trains it with the 
#' training data. Additionally, the model is statistically evaluated based on
#' the z-value, P, Std. Error, and McFadden's pseudo R^2.
#' 
#' @param input_dir The path to the directory containing the training data.

model_training_evaluation <- function(input_dir){
  
    ########## Import Training Data ##########
  
  drought_data <- read.table(list.files(path = input_dir, pattern = ".csv$", recursive = T, 
                                        full.names = T), header = T, sep = ",")
  
  ########## Select Variables ##########
  
  variable_choice <- menu(c("Include all variables" ,"Exclude NDVI", "Exclude NDII", 
                            "Exclude LST", "Exclude SPI3"), 
                            "Choose Variables of Correlation Test", graphics = TRUE)
  
  if (variable_choice == 2){
    drought_data <- drought_data[, -which(names(drought_data) == "NDVI"), drop = FALSE]
    variable_choice <- 1
  }
  if (variable_choice == 3){
    drought_data <- drought_data[, -which(names(drought_data) == "NDII"), drop = FALSE]
    variable_choice <- 1
  }
  if (variable_choice == 4){
    drought_data <- drought_data[, -which(names(drought_data) == "LST"), drop = FALSE]
    variable_choice <- 1
  }
  if (variable_choice == 5){
    drought_data <- drought_data[, -which(names(drought_data) == "SPI3"), drop = FALSE]
    variable_choice <- 1
  }
  
  if (variable_choice ==1){
    
    ########## Build the Model (GLM) ##########
  
    model <- glm(Drought~., family = binomial(link = "logit"), data = drought_data)
  
    ########## Statistical Evaluation ##########
  
    # Calculate the z-value, P, and Std. Error
    sum <- summary(model)
  
    print("Model Summary Statistics: ")
    print(sum)
  
    # Get the McFadden's pseudo R^2 (should be <0.2)
    pseudoR <- pR2(model)
  
    print("McFaddens's Pseudo R^2 (should be <0.2): ")
    print(pseudoR)
  }
}


#' Create Drought Probability Maps
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
#' @param output_dir The path to the folder where the drought probability maps
#'                   will be stored.

drought_probability_maps <- function(start_date, end_date, predicted_years, predicted_months, input_dir, output_dir){
  
  ########## Data Import ##########
  
  # Import model data (NDVI, NDII, and LST Anomalies + SPI3 Precipitation Data)
  model_files <- list.files(path = input_dir, pattern = ".nc$", recursive = T, full.names = T)
  
  model_data <- list()
  for(i in 1:length(model_files)) {
    model_data[[i]] = rast(model_files[i])
  }
  
  # Import Training Data
  drought_data <- read.table(list.files(path = input_dir, pattern = ".csv$", recursive = T, 
                                        full.names = T), header = T, sep = ",")

  ########## Data Preparation ##########
  
  # Extract Months to predict the drought probability
  dates <- seq(as.Date(start_date), as.Date(end_date), by = '1 month')
  
  for(i in 1:length(model_data)){
    model_data[[i]] = subset(model_data[[i]], which(format.Date(dates, "%m") %in% predicted_months))
  }
  
  dates <- dates[format.Date(dates, "%m") %in% predicted_months]
  
  # Extract years to predict the drought probability
  for (i in 1:length(model_data)) {
    model_data[[i]] <- subset(model_data[[i]], which(format.Date(dates, "%Y") %in% predicted_years))
  }
  
  dates <- dates[format.Date(dates, "%Y") %in% predicted_years]
  
  # Create a list of raster stacks with each index as one layer for every month
  predict_stack_list = list()
  
  for(i in 1:nlyr(model_data[[1]])) {
    predict_stack_list[[i]] = c(model_data[[1]][[i]], model_data[[2]][[i]],
                                model_data[[3]][[i]], model_data[[4]][[i]])
    
    names(predict_stack_list[[i]]) = c("LST", "NDII", "NDVI", "SPI3")
  }
  
  ########## Select Variables ##########
  
  variable_choice <- menu(c("Include all variables" ,"Exclude NDVI", "Exclude NDII", 
                            "Exclude LST", "Exclude SPI3"), 
                            "Choose Variables of Correlation Test", graphics = TRUE)
  
  if (variable_choice == 2){
    drought_data <- drought_data[, -which(names(drought_data) == "NDVI"), drop = FALSE]
    for(i in 1:length(predict_stack_list)){
      predict_stack_list[[i]] <- subset(predict_stack_list[[i]], 
                                        subset = which(names(predict_stack_list[[i]]) != "NDVI"))
    }
    variable_choice <- 1
  }
  if (variable_choice == 3){
    drought_data <- drought_data[, -which(names(drought_data) == "NDII"), drop = FALSE]
    for(i in 1:length(predict_stack_list)){
      predict_stack_list[[i]] <- subset(predict_stack_list[[i]], 
                                        subset = which(names(predict_stack_list[[i]]) != "NDII"))
    }
    variable_choice <- 1
  }
  if (variable_choice == 4){
    drought_data <- drought_data[, -which(names(drought_data) == "LST"), drop = FALSE]
    for(i in 1:length(predict_stack_list)){
      predict_stack_list[[i]] <- subset(predict_stack_list[[i]], 
                                        subset = which(names(predict_stack_list[[i]]) != "LST"))
    }
    variable_choice <- 1
  }
  if (variable_choice == 5){
    drought_data <- drought_data[, -which(names(drought_data) == "SPI3"), drop = FALSE]
    for(i in 1:length(predict_stack_list)){
      predict_stack_list[[i]] <- subset(predict_stack_list[[i]], 
                                        subset = which(names(predict_stack_list[[i]]) != "SPI3"))
    }
    variable_choice <- 1
  }
  
  if (variable_choice ==1){
    
    ########## Build Model ##########
  
    model <- glm(Drought~., family = binomial(link = "logit"), data = drought_data)
  
    ########## Create Probability Maps ##########
  
    # List for the created drought probability rasters
    prob_pred = list()
  
    for(i in 1:length(predict_stack_list)) {
      prob_pred[[i]] = terra::predict(predict_stack_list[[i]], model, type = "response")
    }
  
    # Convert the list to a raster stack
    prob_pred = rast(prob_pred)
  
    ########## Write Results ##########
  
    # Set names for output files
    out_name <- "/Drought_Prob_SegmRegr_"
    names_all <- NULL
  
    for(i in 1:length(dates)) {
      names_all[[i]] = paste(output_dir, out_name, dates[i], ".tif", sep = "")  
    }
  
    # Write probability maps as .tif files
    for(i in 1:nlyr(prob_pred)) {
      terra::writeRaster(prob_pred[[i]], filename = names_all[[i]], 
                         filetype = "GTiff", overwrite = TRUE)
    }
  }
}