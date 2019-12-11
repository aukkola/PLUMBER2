
flux_corrections <- function(flux_nc, qc_info, outfile_flux, new_qc, 
                             qle_name="Qle", qh_name="Qh", rnet_name="Rnet", qg_name="Qg",
                             qle_cor_name="Qle_cor", qh_cor_name="Qh_cor")

{


  ################################################
  #########------ Flux corrections ------#########
  ################################################
  
  
  #######################
  ### Get timing info ###
  #######################
  
  
  #Get time vector
  time <- ncvar_get(flux_nc, "time")
  
  #Get time units
  time_units <- strsplit(ncatt_get(flux_nc, "time")$units, "seconds since ")[[1]][2]
  
  #Convert to Y-M-D h-m-s
  time_date <- as.POSIXct(time, origin=time_units, tz="GMT")
  
  #Get time interval (in fraction of day)
  tsteps_per_day <-  60*60*24 / (time[2] - time[1]) 
  
  
  #Find adjusted start and end time
  years <- as.numeric(format(time_date, "%Y"))
  
  
  
  
  #####################
  ### Get variables ###
  #####################
  
  
  ### Get all data variables with a time dimension ###
  
  #Get variable names
  vars <- names(flux_nc$var)
  
  #Load variable data
  var_data <- lapply(vars, function(x) ncvar_get(flux_nc, x))
  
  #Set names
  names(var_data) <- vars 
  
  
  #Get variable attributes
  att_data <- lapply(vars, function(x) ncatt_get(flux_nc, x))
  
  #Set names
  names(att_data) <- vars 
  
  
  
  
  
  #################################
  ### Energy balance correction ###
  #################################
  
  #If energy balance corrected values already exist (Fluxnet2015), skip this step
  if (!(qle_cor_name %in% vars)) {
    
      
    #Check that have all variables available
    if (all(c(qle_name, qh_name, rnet_name, qg_name) %in% vars)) {
      
      #Get corrected fluxes
      ebcf_corrected <- energy_balance_correction(qle=var_data[qle_name], 
                                                  qle_qc=var_data[paste0(qle_name, "_qc")], 
                                                  qh=var_data[qh_name], 
                                                  qh_qc=var_data[paste0(qh_name, "_qc")], 
                                                  rnet=var_data[rnet_name], 
                                                  qg=var_data[qg_name], 
                                                  qg_qc=var_data[paste0(qg_name, "_qc")], 
                                                  time=time_date, tstepsize=time[2] - time[1])
    }
  }
 
  
  
  
  #Add new variables to netcdf handle
  
  
  qle_cor_name
  "W/m2"
  
  
  
  
  
  
  ##########################
  ### Adjust time period ###
  ##########################
  
  #If need to adjust
  if (start_yr > 1 | end_yr < 0) {
    
    
    ### Adjust length of time-varying variables ##
    
    #Get dimensions for each variable
    dims <- lapply(vars, function(x) sapply(flux_nc[[s]][["var"]][[x]][["dim"]], function(dim) dim[["name"]]))
    
    # #Find which variables are time-varying
    var_inds <- which(sapply(dims, function(x) any(x == "time")))
    
    
    
    #Change dimensions and values for time-varying data
    for (v in vars[var_inds]) {
      
      #Change time dimension
      flux_nc$var[[v]]$varsize[3] <- length(time_var)
      
      #Change time values
      flux_nc$var[[v]]$dim[[3]]$vals <- time_var
      
      #Change time size
      flux_nc$var[[v]]$dim[[3]]$len <- length(time_var)
      
      #Change length
      flux_nc$var[[v]]$size[3] <- length(time_var)
      
      #Change values in var_data
      var_data[[v]] <- var_data[[v]][start_ind:end_ind]
      
      # #Change chunk size (no idea what this is but produces an error otherwise
      # #during nc_create)
      # met_nc[[s]]$var[[v]]$chunksizes <- NA
    }
    
    
    #Also adjust time dimension and units
    
    #Change time dimensions
    #Change values, length and unit
    flux_nc$dim$time$vals  <- time_var
    flux_nc$dim$time$len   <- length(time_var)
    
    flux_nc$dim$time$units <- paste0("seconds since ", new_start_year, "-01-01 00:00:00")
    
  }
  
  
  
  ##################### 
  ### Re-write file ###
  #####################
  
  
  ###--- Set dimensions ---###
  
  #Get dimensions from input file
  new_dims <- flux_nc$dim
  
  
  ###--- Define variables ---###
  
  #Get variables from input file
  new_vars <- flux_nc$var
  
  
  ###--- Set up new file ---###
  
  #New file handle
  out_nc <- nc_create(outfile_flux, vars=new_vars)
  
  
  ###--- Global attributes ---###
  
  #Get global attributes
  global_atts <- ncatt_get(flux_nc, varid=0)
  
  #Add new QC flag value to metadata    
  global_atts$QC_flag_descriptions <- paste0(global_atts$QC_flag_descriptions, 
                                             ", Post-processed: ", new_qc)
  
  #Add to file
  #For some reason this crashes if using lapply, loop works ok-ish    
  for(a in 1:length(global_atts)){
    ncatt_put(out_nc, varid=0, attname=names(global_atts)[a], 
              attval=unlist(global_atts[a]))
  }
  
  
  ###--- Variable data ---###
  
  #Write variables to output file
  for (v in names(new_vars)) {
    
    ncvar_put(nc=out_nc, varid=new_vars[[v]],
              vals=var_data[[v]])
  }
  
  
  ###--- Variable attributes ---###
  
  #Write attributes to output file
  for (v in names(att_data)) {
    
    for (a in names(att_data[[v]]))
      
      ncatt_put(nc=out_nc, varid=new_vars[[v]],
                attname=a, attval=att_data[[v]][[a]])
    
  }
  
  
  
  #Close output file
  nc_close(out_nc)
  
  #Close original file handle
  nc_close(flux_nc)
  

} #function