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
    
    
    
    
  } else if (site_code == "ES-ES2") {
    
    #One day missing in met data, gapfill with previous day
    vars_to_fix <- c("Tair", "VPD", "Qair", "Precip", "SWdown", "Wind", "Psurf")
    
    for (v in vars_to_fix) {
      
      #Missing time steps
      na_ind <- which(is.na(var_data[[v]]))
      
      var_data[[v]][na_ind]    <- var_data[[v]][na_ind - (length(na_ind))]
      var_data[[paste0(v, "_qc")]][na_ind] <- qc_val
      
      #Add new gapfilled % to attribute data 
      att_data[[v]]$`Gap-filled_%` <- round(length(which(var_data[[paste0(v, "_qc")]] > 0)) /
                                              length(var_data[[paste0(v, "_qc")]]) * 100, digits=1)
      
    }
  
    
    
    
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
    var_data$RH_qc[na_ind] <- qc_val
    
    #Add new gapfilled % to attribute data 
    att_data$RH$`Gap-filled_%` <- round(length(which(var_data$RH_qc > 0)) /
                                              length(var_data$RH_qc) * 100, digits=1)
    
    
    
    #-- Synthesise LWdown
    
    #Check that Tair units are in Kelvin and that RH is available
    if (att_data$Tair$units != "K") stop("Wrong LWdown units")
    if (length(which(vars == "RH")) == 0) stop("RH not available")
    
    lwdown <- SynthesizeLWdown(TairK = var_data$Tair, 
                               RH = var_data$RH)
    
    #Add corrected flux to variable data and QC flag
    na_ind <- which(is.na(var_data$LWdown))
    var_data$LWdown[na_ind]    <- lwdown[na_ind]  
    var_data$LWdown_qc[na_ind] <- qc_val
    
    #Add new gapfilled % to attribute data 
    att_data$LWdown$`Gap-filled_%` <- round(length(which(var_data$LWdown_qc > 0)) /
                                              length(var_data$LWdown_qc) * 100, digits=1)
    
  
    
  } else if (site_code == "US-ARM") {
    #Need to re-gapfill parts of LWdown
    
    inds <- c(20500:71000, 90500:93500)

    #Check that Tair units are in Kelvin and that RH is available
    if (att_data$Tair$units != "K") stop("Wrong LWdown units")
    if (length(which(vars == "RH")) == 0) stop("RH not available")
    
    lwdown <- SynthesizeLWdown(TairK = var_data$Tair[inds], 
                               RH = var_data$RH[inds])
    
    #Add corrected flux to variable data and QC flag
    var_data$LWdown[inds]    <- lwdown[inds]  
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
    var_data$LWdown[ind]    <- lwdown[ind]  
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

    
  #Any sites with missing CO2 values
    
  #Missing CO2 values
  if (any(is.na(var_data$CO2air))) {
    
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
  
    
  #Missing RH values 
    
  #(some sites have missing RH but all VPD so convert between the two)
  if (any(is.na(var_data$RH))) {
    
    rh <- VPD2RelHum(VPD=var_data$VPD, airtemp=var_data$Tair, 
                     vpd_units=att_data$VPD$units, tair_units=att_data$Tair$units)
    
    #Replace missing values
    na_ind <- which(is.na(var_data$RH))
    
    var_data$RH[na_ind]    <- rh[na_ind]
    var_data$RH_qc[na_ind] <- qc_val
    
    #Add new gapfilled % to attribute data 
    att_data$RH$`Gap-filled_%` <- round(length(which(var_data$RH_qc > 0)) /
                                          length(var_data$RH_qc) * 100, digits=1)
    
  }
  
    
  #Missing LWdown values
    
  if (any(is.na(var_data$LWdown))) {
    
    #Check that Tair units are in Kelvin and that RH is available
    if (att_data$Tair$units != "K") stop("Wrong LWdown units")
    if (length(which(vars == "RH")) == 0) stop("RH not available")
    
    lwdown <- SynthesizeLWdown(TairK = var_data$Tair, 
                               RH = var_data$RH)
    
    #Replace missing with predicted and QC with post-processing value
    na_ind <- which(is.na(var_data$LWdown))
    
    var_data$LWdown[na_ind]    <- lwdown[na_ind]
    var_data$LWdown_qc[na_ind] <- qc_val
    
    #Add information to metadata
    att_data$LWdown$LWdown_correction <- "Synthesised from RH and Tair (Abramowitz method)"
    
    #Add new gapfilled % to attribute data 
    att_data$LWdown$`Gap-filled_%` <- round(length(which(var_data$LWdown_qc > 0)) /
                                              length(var_data$LWdown_qc) * 100, digits=1)
    
  }

    
    
  #Negative LAI
  lai_vars <- vars[which(grepl("LAI_", vars))]
  
  #Loop through products
  for (lai in lai_vars) {
    
    #Replace negative values with zero
    var_data[[lai]] <- replace(var_data[[lai]], var_data[[lai]] < 0, 0)
  }  
    
  

  #If not fixes, returns original data
  
  return(list(var_data=var_data, att_data=att_data))
  
  
} #function






