#==============================================================================#
#' Calculation of the Standardized Precipitation Index (SPI) from TAMSAT       #
#' rainfall estimates.                                                         #
#'                                                                             #
#' Download the TAMSAT Data Precipitation Data from the Website:               #
#' http://www.tamsat.org.uk/index.php/data > "Time and area subsetting tool"   #
#' Tamsat v3.1 Monthly, Regional data (NetCDF), Region: Kenya                  #
#' Start Date: 1/1983                                                          #
#' End Date: 7/2023                                                            #
#                                                                              #
#==============================================================================#


#==============================================================================#
################################ Functions #####################################
#                                                                              #
#' Calculate Standardized Precipitation Index(SPI)                             #
#'                                                                             #
#' This function is used to calculated the SPI from                            #
#' TAMSAT monthly regional precipitation data.                                 #
#'                                                                             #
#==============================================================================#

Calculate_SPI_from_TAMSAT <- function(input_file, output_dir){
  
  ########## Data Import ##########
  precip <- rast(input_file)
  
  ########## Pre-Processing ##########
  # Convert -Inf to NA to avoid issues in calculations
  precip[precip == -Inf] <- NA
  
  # Extract dates from raster time attribute
  dates <- as.Date(terra::time(precip))
  
  ########## SPI Calculation ##########
  # Define a function for SPI calculation
  calcspi_3 <- function(x) {
    # Replace non-positive values with a small positive value
    spi_output <- spi(x, scale=3, distribution="Gamma", na.rm=TRUE)
    if(is.null(spi_output$fitted)) {
      rep(NA, length(x))  # Return NAs if calculation fails
    } else {
      return(as.numeric(spi_output$fitted))
    }
  }
  
  # Setup parallel processing
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, list("spi"))
  clusterEvalQ(cl, library(SPEI))  # ensure SPEI is loaded in each worker
  
  # Apply SPI calculation in parallel
  spiRST_3 <- app(precip, fun=calcspi_3, cores=cl)
  
  # Cleanup
  stopCluster(cl)
  
  ########## Write Output ##########
  out_name <- vector("character", length=nlyr(spiRST_3))
  
  for (i in 1:nlyr(spiRST_3)) {
    date_str <- format(dates[i], "%Y-%m")
    out_name[i] <- file.path(output_dir, paste0(date_str, "_SPI_3.tif"))
  }
  
  # Writing files
  lapply(1:nlyr(spiRST_3), function(i) {
    writeRaster(spiRST_3[[i]], filename=out_name[i], filetype="GTiff", overwrite=TRUE)
  })
  
}
