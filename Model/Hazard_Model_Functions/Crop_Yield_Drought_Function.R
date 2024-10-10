#===============================================================================#
#' Script to identify drought years with FAO crop yield data                    #
#'                                                                              #
#' This function is used to identify drought years with yield data from the     #
#' Food and Agriculture Organization of the United Nations (FAO). For every     #
#' crop type within each country, the function calculates the linear regression #
#' between yield and time. To take into account yield increases due to          #
#' technical enhancements, breakpoints are calculated for each regression       #
#' with a minimum distance of four years. Residuals of all crop types are       #
#' summed up for each year and each country. If the standard deviation of       #
#' residuals is < -1, the year is considered as a drought year. The final       #
#' output is a table indicating drought (1) or no drought (0) for every year.   #
#'                                                                              #
#' FAO crop yield data can be downloaded as .csv file from:                     #
#' https://www.fao.org/faostat/en/#data/QCL                                     #
#' in "elements", select "yield"                                                #
#' in "items", select four to five main food crops for the country              #
#' in "Years", select years 1990 to the most current                            #
#'                                                                              #
#'                                                                              #
#' The user needs to specify:                                                   #
#' @param crop_data_path The path to the FAO crop yield data.                   #
#' @param drought_years_fao The path to save the final output file.             #
#'==============================================================================#


#install and load required packages
#install.packages("strucchange")
#install.packages("segmented")

#library(strucchange)
#library(segmented)

########## Paths to define #########
#crop_data_path = "D:/Katharina/04_ADM_Kenya/Programming/Total_Model/Training_Data/FAOSTAT_data_Kenya.csv"
#output_file = "D:/Katharina/04_ADM_Kenya/Programming/Total_Model/Training_Data/Kenya_Crop_Yield_Drought_test.csv"

########## Function ###############
get_drought_years_from_yield_data <- function (crop_data_path, drought_years_fao) {

  #Read data
  crop_data = read.table(crop_data_path, header = T, sep = ",")

  #Create data frame with only necessary columns
  crop_data = data.frame(crop_data$Area, crop_data$Year, crop_data$Item, crop_data$Unit, crop_data$Value)

  #Rename columns
  names(crop_data) = c("Country", "Year", "Crop", "Unit", "Yield")

  #Extract countries and crop types
  countries = unique(crop_data$Country)
  crop_types = unique(crop_data$Crop)

  #create empty lists for loop
  res_tot = list()
  res_country = list()

  #Start loop for each country
  for(j in 1:length(countries)){
    #get new dataframe with only crops for current country
    crops_country = subset(crop_data, crop_data$Country == countries[j])
  
    #loop for each crop type within current country
    for(i in 1:length(crop_types)){
      #subset current crop type
      crops_country_type = subset(crops_country, crops_country$Crop == crop_types[i])
    
      #create breaking points for linear segmented regression analysis for current country
      bp = breakpoints(Yield ~ Year, h = 4, data = crops_country_type)
    
      #set up linear model
      reg = lm(Yield ~ Year, data = crops_country_type)
    
      #if no break points for current crop run get residuals from simple linear regression
      ifelse(is.na(breakpoints(bp)$breakpoints),
        
        (res_country[[i]] = residuals(reg)),
        
        #if breakpoints available then run segmented regression -> dependent on number of breakpoints and then get residuals
      
        if(length(breakpoints(bp)$breakpoints) == 1){
          seg = segmented(reg, seg.Z = ~ Year, psi = 2008)
          res_country[[i]] = residuals(seg)
        } else if(length(breakpoints(bp)$breakpoints) == 2){
          seg2 = segmented(reg, seg.Z = ~ Year, psi = c(2002,2015))
          res_country[[i]] = residuals(seg2)
        } else if(length(breakpoints(bp)$breakpoints) == 3){
          seg3 = segmented(reg, seg.Z = ~ Year, psi = c(2002,2008,2015))
          res_country[[i]] = residuals(seg3)
        } else if(length(breakpoints(bp)$breakpoints) == 4){
          seg4 = segmented(reg, seg.Z = ~ Year, psi = c(1994,2002,2008,2015))
          res_country[[i]] = residuals(seg4)
        } else {
          seg5 = segmented(reg, seg.Z = ~ Year, psi = c(1994,1998,2002,2008,2015))
          res_country[[i]] = residuals(seg5)
        })
        
      
    }
    
  
    #create dataframe of residuals
    df = data.frame(matrix(unlist(res_country), nrow = length(res_country[[1]]), byrow = F))
  
    #sum up residuals of all crop types per year
    res_tot[[j]] = rowSums(df)
  }

  #create dataframe of list from loop -> residuals per country and year
  res_df = data.frame(matrix(unlist(res_tot), nrow = length(res_tot[[1]]), byrow = F))

  #changes column names of dataframe
  names(res_df) = countries

  #Drought analysis: If standard deviation of residuals of one year > 1 -> non drought year, if standard deviation of residuals < -1 -> drought year
  #create empty list for loop
  drought = rep(list(c()),length(countries))
  for(i in 1:ncol(res_df)) {
    #calculate standard deviation for residuals of segmented regression analysis for current year
    std_dev = sd(res_df[,i])
  
    #run loop for the residual analysis of each year
    for(j in 1:nrow(res_df)) {
      #If residuals < -1 standard dev -> drought year -> classify as "1"
      if(res_df[j,i] <= -std_dev) {
        drought[[i]][j] = 1
        #If residuals > 1 standard deviation -> non drought year -> classify as "0"
      } else if(res_df[j,i] >= std_dev) {
        drought[[i]][j] = 0
        #If residuals within -1 and 1 standard deviation -> mixed year -> classify as "NA" (cannot be used as training data for the model)
      } else {
        drought[[i]][j] = NA
      }
    }
  }

  #create dataframe of loop with indication of drought (1) or no drought (0) for each year
  drought_df = data.frame(matrix(unlist(drought), nrow = length(drought[[1]]), byrow = F))

  #Change names of data frame
  names(drought_df) = countries

  #add year column to data frame
  drought_df$Year = unique(crop_data$Year)

  #write dataframe to table
  write.table(drought_df, file = paste0(drought_years_fao, "/Kenya_Crop_Yield_Drought.csv"), sep = ";", na = "NA", 
              col.names = TRUE, row.names = FALSE)
}

###################test: execute function
#get_drought_years_from_yield_data(crop_data_path, output_file)
