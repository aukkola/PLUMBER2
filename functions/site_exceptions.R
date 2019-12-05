site_exceptions <- function(site_code, var_data, att_data, qc_val) {
  
  
  if (site_code == "BE-Lon") {
    #Longwave data needs fixing    
    
    
    
    browser()
    
    
    
    
  } else if (site_code == "BE-Vie") {
    #CO2 drift correction, interpolate back from  xxxx
    
    
    
    browser()
    
    
    
  } else if (site_code == "DE-Geb") {
    #Missing CO2
  
    
    browser()
    
    
    
  } else if (site_code == "DK-Fou") {
    #Gapfill missing CO2 using existing site data
    
    
    
    browser()
    
  
  
  } else if (site_code == "DK-Sor") {
    #Correct CO2 by fitting a linear slope to data, and using that to predict CO2
    
    
    browser()
    
    
  } else if (site_code == "RU-Fyo") {
    #Gapfill LWdown from XX onwards using data from previous years
    
    
    browser()
    
    
  } else if (site_code == "RU-Zot") {
    #Gaps in LWdown and CO2, need to gapfill
    
    
    browser()
    
    
    
  } else if (site_code == "SD-Dem") {
    #Gapfill first year of CO2
    
    
    browser()
    
    
    
  } else if (site_code == "US-ARM") {
    #Need to re-gapfill parts of LWdown
    
    
    browser()
    
    
    
    
  } else if (site_code == "US-NR1") {
    #Re-gapfill parts ofLWdown, and missing CO2 
    
    #CO2: fit linear trend to existing data and use to gapfill
    
    
    browser()
    
    
    
    
  } else if (site_code == "US-PFa") {
    
    #Missing CO2 values, gapfill using a linear fit on the available time series
    
    co2 <- var_data$CO2air
    
    #Linear model
    x  <- 1:length(co2)
    lm <- lm(co2 ~ x)
    
    #Linear prediction
    co2_pred <- lm$coefficients[2] * x + lm$coefficients[1]
    
    #Replace missing with predicted and QC with post-processing value
    na_ind <- which(is.na(co2))
    
    var_data$CO2air[na_ind]    <- co2_pred[na_ind]
    var_data$CO2air_qc[na_ind] <- qc_val
    
    #Add information to metadata
    att_data$CO2air$CO2_correction <- "Linear fit"
    
    
    
  } 

  #If not fixes, returns original data
  
  return(list(var_data=var_data, att_data=att_data))
  
  
} #function






