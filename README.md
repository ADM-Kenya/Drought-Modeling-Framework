# Drought Model

The modeling framework for the drought model consists of multiple R-based modules and functions that work together to provide outputs for drought hazard, drought vulnerability and drought risk. Minimal data input is required by the user. The data for satellite-based indices - the **Normalized Difference Vegetation Index (NDVI), Normalized Difference Infrared Index (NDII), Land Surface Temperature (LST), 3-monthly Standardized Precipitation Index (SPI3)** - that are used as predictors by the model will be downloaded and preprocessed automatically. Overall, the R-based modelling framework covers the download and preprocessing of the Sentinel-3 data, the calibration to MODIS data to ensure the use of a longer time series, training data preparation for the model, model training and application and drought vulnerability and risk assessment. The final output consists of monthly drought hazard and drought risk maps along with a static drought vulnerability layer with a spatial resolution of ~1km.


### **Requirements**
 - System Requirements: Windows 10/11 64 bits, not tested on Linux. 
 - At least 32GB of RAM. 
 - R-Studio, Anaconda, GDAL
 - Dependencies: defined in the environment file


### **Input**
Format needs to geotiff and ESRI shapefile.
 - Study area: Shapefile of the Study Area. 
 - Vulnerability input data: All input layers need to be formatted to GTiffs.
 - Crop yield data: CSV file provided through FAOSTAT.
