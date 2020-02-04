site_exceptions <- function(site_code, var_data, att_data, qc_val) {
  
    #NB. all indices to fix were found manually
  
    #Variables available
    vars <- names(var_data)
  
  
  if (site_code == "BE-Lon") {
    #Longwave data needs fixing    
  
    ind_start <- 74000
    ind_end   <- 92500
    
    #Check that Tair units are in Kelvin and that RH is available
    if (att_data$Tair$units != "K") stop("Wrong LWdown units")
    if (length(which(vars == "RH")) == 0) stop("RH not available")
    
    lwdown <- SynthesizeLWdown(TairK = var_data$Tair[ind_start:ind_end], 
                               RH = var_data$RH[ind_start:ind_end])
    
    #Add corrected flux to variable data and QC flag
    var_data$LWdown[ind_start:ind_end]    <- lwdown  
    var_data$LWdown_qc[ind_start:ind_end] <- qc_val 
    
    #Add new gapfilled % to attribute data 
    att_data$LWdown$`Gap-filled_%` <- round(length(which(var_data$LWdown_qc > 0)) /
                                      length(var_data$LWdown_qc) * 100, digits=1)
    
    
    
    
  } else if (site_code == "BE-Vie") {
    #CO2 drift correction, interpolate back from time step 60100
    
    #Interpolate back from this point
    ind <- 60100
    
    #Linear prediction
    co2_pred <- linear_pred_co2(co2=var_data$CO2air, start_ind=ind+1, 
                                end_ind=length(var_data$CO2air))
    
 
    #Replace predicted and QC with post-processing value
    var_data$CO2air[1:ind]    <- co2_pred[1:ind]
    var_data$CO2air_qc[1:ind] <- qc_val
    
    #Add information to metadata
    att_data$CO2air$CO2_correction <- "Linear fit"
    
    #Add new gapfilled % to attribute data 
    att_data$CO2air$`Gap-filled_%` <- round(length(which(var_data$CO2air_qc > 0)) /
                                              length(var_data$CO2air_qc) * 100, digits=1)
    
    
    
    
  } else if (site_code %in% c("DE-Geb","DK-Fou","US-PFa")) {
    #Missing CO2 values, gapfill using a linear fit on the available time series
    
    co2 <- var_data$CO2air
    
    #Linear prediction
    co2_pred <- linear_pred_co2(co2=co2, start_ind=1, end_ind=length(co2))
    
    #Replace missing with predicted and QC with post-processing value
    na_ind <- which(is.na(co2))
    
    var_data$CO2air[na_ind]    <- co2_pred[na_ind]
    var_data$CO2air_qc[na_ind] <- qc_val
    
    #Add information to metadata
    att_data$CO2air$CO2_correction <- "Linear fit, gapfilling missing time steps"
    
    #Add new gapfilled % to attribute data 
    att_data$CO2air$`Gap-filled_%` <- round(length(which(var_data$CO2air_qc > 0)) /
                                              length(var_data$CO2air_qc) * 100, digits=1)
    
    
  
  } else if (site_code == "DK-Sor") {
    #Correct CO2 by fitting a linear slope to data, and using that to predict CO2
    
    #Linear prediction
    co2_pred <- linear_pred_co2(co2=var_data$CO2air, start_ind=1, 
                                end_ind=length(var_data$CO2air))
    
    #Replace all values
    var_data$CO2air    <- co2_pred
    var_data$CO2air_qc <- qc_val
    
    #Add information to metadata
    att_data$CO2air$CO2_correction <- "Linear fit, correcting all time steps"
    
    #Add new gapfilled % to attribute data 
    att_data$CO2air$`Gap-filled_%` <- round(length(which(var_data$CO2air_qc > 0)) /
                                              length(var_data$CO2air_qc) * 100, digits=1)
    
    
  } else if (site_code == "RU-Fyo") {
    #Gapfill LWdown from XX onwards using data from previous years
    
    ind_start <- 144000
    ind_end   <- 277000
    
    #Check that Tair units are in Kelvin and that RH is available
    if (att_data$Tair$units != "K") stop("Wrong LWdown units")
    if (length(which(vars == "RH")) == 0) stop("RH not available")
    
    lwdown <- SynthesizeLWdown(TairK = var_data$Tair[ind_start:ind_end], 
                               RH = var_data$RH[ind_start:ind_end])
    
    #Add corrected flux to variable data and QC flag
    var_data$LWdown[ind_start:ind_end]    <- lwdown  
    var_data$LWdown_qc[ind_start:ind_end] <- qc_val 
    
    #Add new gapfilled % to attribute data 
    att_data$LWdown$`Gap-filled_%` <- round(length(which(var_data$LWdown_qc > 0)) /
                                              length(var_data$LWdown_qc) * 100, digits=1)
    
    
    
  } else if (site_code == "RU-Zot") {
    #Gaps in LWdown, RH and CO2, need to gapfill
    
    #-- Gapfill RH by converting VPD
    
    rh <- VPD2RelHum(VPD=var_data$VPD, airtemp=var_data$Tair, 
                     vpd_units=att_data$VPD$units, tair_units=att_data$Tair$units)
      
    #Replace missing values
    na_ind <- which(is.na(var_data$RH))
    
    var_data$RH[na_ind]    <- rh[na_ind]
    var_data$RH_qc[na_ind] <- qc_val[na_ind]
    
    
    
    #-- Synthesise LWdown
    
    #Check that Tair units are in Kelvin and that RH is available
    if (att_data$Tair$units != "K") stop("Wrong LWdown units")
    if (length(which(vars == "RH")) == 0) stop("RH not available")
    
    lwdown <- SynthesizeLWdown(TairK = var_data$Tair, 
                               RH = var_data$RH)
    
    #Add corrected flux to variable data and QC flag
    na_ind <- which(is.na(var_data$LWdown))
    var_data$LWdown[na_ind]    <- lwdown[na_ind]  
    var_data$LWdown_qc[na_ind] <- qc_val[na_ind]
    
    #Add new gapfilled % to attribute data 
    att_data$LWdown$`Gap-filled_%` <- round(length(which(var_data$LWdown_qc > 0)) /
                                              length(var_data$LWdown_qc) * 100, digits=1)
    
    
    #-- Gapfill CO2

    #Linear prediction
    co2_pred <- linear_pred_co2(co2=var_data$CO2air, start_ind=1, 
                                end_ind=length(var_data$CO2air))
    

    #Replace all values with predicted and QC with post-processing value
    #(looks dodgy if only gapfill missing tsteps)
    var_data$CO2air    <- co2_pred[na_ind]
    var_data$CO2air_qc <- qc_val
    
    #Add information to metadata
    att_data$CO2air$CO2_correction <- "Linear fit, correcting all time steps"
    
    #Add new gapfilled % to attribute data 
    att_data$CO2air$`Gap-filled_%` <- round(length(which(var_data$CO2air_qc > 0)) /
                                              length(var_data$CO2air_qc) * 100, digits=1)
    
    
    
  } else if (site_code == "SD-Dem") {
    #Gapfill first year of CO2
    
    ind <- length(var_data$CO2air)/2 #corresponds to end of first yr
    
    #Linear prediction
    co2_pred <- linear_pred_co2(co2=var_data$CO2air, start_ind=ind+1, 
                                end_ind=length(var_data$CO2air))
    
    
    #Replace all values with predicted and QC with post-processing value
    #(looks dodgy if only gapfill missing tsteps)
    var_data$CO2air[1:ind]    <- co2_pred[1:ind]
    var_data$CO2air_qc[1:ind] <- qc_val
    
    #Add information to metadata
    att_data$CO2air$CO2_correction <- "Linear fit"
    
    #Add new gapfilled % to attribute data 
    att_data$CO2air$`Gap-filled_%` <- round(length(which(var_data$CO2air_qc > 0)) /
                                              length(var_data$CO2air_qc) * 100, digits=1)
    
    
  } else if (site_code == "US-ARM") {
    #Need to re-gapfill parts of LWdown
    
    inds <- c(20500:71000, 90500:93500)

    #Check that Tair units are in Kelvin and that RH is available
    if (att_data$Tair$units != "K") stop("Wrong LWdown units")
    if (length(which(vars == "RH")) == 0) stop("RH not available")
    
    lwdown <- SynthesizeLWdown(TairK = var_data$Tair[inds], 
                               RH = var_data$RH[inds])
    
    #Add corrected flux to variable data and QC flag
    var_data$LWdown[inds]    <- lwdown  
    var_data$LWdown_qc[inds] <- qc_val 
    
    #Add new gapfilled % to attribute data 
    att_data$LWdown$`Gap-filled_%` <- round(length(which(var_data$LWdown_qc > 0)) /
                                              length(var_data$LWdown_qc) * 100, digits=1)
    
    
    
  } else if (site_code == "US-NR1") {
    
    #Re-gapfill parts ofLWdown, and missing CO2 
    
    #-- Re-gapfill LWdown
    ind <- c(107500:137000)
    
    #Check that Tair units are in Kelvin and that RH is available
    if (att_data$Tair$units != "K") stop("Wrong LWdown units")
    if (length(which(vars == "RH")) == 0) stop("RH not available")
    
    lwdown <- SynthesizeLWdown(TairK = var_data$Tair, 
                               RH = var_data$RH)
    
    #Add corrected flux to variable data and QC flag
    var_data$LWdown[ind]    <- lwdown  
    var_data$LWdown_qc[ind] <- qc_val 
    
    #Add new gapfilled % to attribute data 
    att_data$LWdown$`Gap-filled_%` <- round(length(which(var_data$LWdown_qc > 0)) /
                                              length(var_data$LWdown_qc) * 100, digits=1)
    
    
    
    #-- CO2: fit linear trend to existing data and use to gapfill
    
    co2 <- var_data$CO2air
    
    #Linear prediction
    co2_pred <- linear_pred_co2(co2=co2, start_ind=1, end_ind=length(co2))
    
    #Replace missing with predicted and QC with post-processing value
    na_ind <- which(is.na(co2))
    
    var_data$CO2air[na_ind]    <- co2_pred[na_ind]
    var_data$CO2air_qc[na_ind] <- qc_val
    
    #Add information to metadata
    att_data$CO2air$CO2_correction <- "Linear fit, gapfilling missing time steps"
    
    #Add new gapfilled % to attribute data 
    att_data$CO2air$`Gap-filled_%` <- round(length(which(var_data$CO2air_qc > 0)) /
                                              length(var_data$CO2air_qc) * 100, digits=1)
    

    
  } 

  #If not fixes, returns original data
  
  return(list(var_data=var_data, att_data=att_data))
  
  
} #function






