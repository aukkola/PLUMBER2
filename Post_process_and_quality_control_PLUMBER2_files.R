#Run this code after processing data through FluxnetLSM and adding 
#LAI time series from MODIS and Copernicus


library(ncdf4)
library(gsheet)


path <- "/srv/ccrc/data04/z3509830/Fluxnet_data/All_flux_sites_processed_PLUMBER2/"


#Source function
source(paste0(path, "/scripts/functions/site_exceptions.R"))



###########################
### Read QC information ###
###########################

#Read QC information from google sheet,
#e.g. start and end year and CO2 processing


qc_info <- gsheet2tbl("https://docs.google.com/spreadsheets/d/1bi9bbUpwzRycDJ16VTYGx8ApkDLV_IeSwEsdDJsX1SI/edit#gid=0")

qc_sites <- qc_info$Site_code



#####################
### Read CO2 data ###
#####################

#Use Mauna Loa annual CO2 record to gapfill CO2 records that are
#poor quality or missing


global_co2 <- read.table(paste0(path, "/Global_CO2_data/co2_annmean_mlo.txt"))


#Set new QC values for post-processing
#Use 101 as OzFlux using a lot of smaller values
new_qc <- 101

                       

######################
### Load site data ###
######################


#Met data

#Get sites
met_files <- list.files(paste0(path, "/all_sites_no_duplicates/Nc_files/Met_with_LAI/"), full.names=TRUE)

#Open file handles
met_nc <- lapply(met_files, nc_open, write=TRUE)


#Flux data

#Get sites
flux_files <- list.files(paste0(path, "/all_sites_no_duplicates/Nc_files/Flux/"), full.names=TRUE)

#Open file handles
flux_nc <- lapply(flux_files, nc_open, write=TRUE)



#Get site codes
site_codes <- sapply(met_nc, function(x) ncatt_get(x, varid=0, "site_code")$value)





### Loop through sites ###

for (s in 1:length(site_codes)) {
  
  
  #############################
  ### Check if keeping file ###
  #############################
  
  #site index for google sheet
  site_ind <- which(qc_sites == site_codes[s])
  
  #Get decision from google sheet
  status <- qc_info$decision[site_ind]
  
  
  ### If site not in google sheet ###
  
  if (length(site_ind) < 1) {
    
    stop(paste0("Site ", site_codes[s], " not found in google sheet"))

  
  ### If not keeping site ###
    
  } else if (status == "kill") {
    
    next
    
    
  ### Else process site ###
    
  } else {
    

    
    #########------ Met corrections ------#########
    
    
    #######################
    ### Get timing info ###
    #######################
    
    
    #Get time vector
    time <- ncvar_get(met_nc[[s]], "time")
    
    #Get time units
    time_units <- strsplit(ncatt_get(met_nc[[s]], "time")$units, "seconds since ")[[1]][2]
    
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
    vars <- names(met_nc[[s]]$var)
    
    #Load variable data
    var_data <- lapply(vars, function(x) ncvar_get(met_nc[[s]], x))
    
    #Set names
    names(var_data) <- vars 
    
    
    #Get variable attributes
    att_data <- lapply(vars, function(x) ncatt_get(met_nc[[s]], x))
    
    #Set names
    names(att_data) <- vars 
    
    
    
    
    ########################################
    ### Correct CO2 using global records ###
    ########################################
    
    
    #If replacing CO2 with global CO2
    
    if (qc_info$Global_CO2[site_ind]) {
     
      #loop through years
      co2_ts <- sapply(years, function(x) global_co2[which(global_co2 == x), 2])
      
      #Replace values in CO2 variable
      var_data$CO2air <- co2_ts 
      
      #Set QC values to missing
      var_data$CO2air_qc <- rep(NA, length(time))  

      #Replace gapfill percentage (now 100%)
      att_data$CO2air["Gap-filled_%"] <- 100
      
      
    }
    
     
    
    
    #############################################
    ### Check for additional site corrections ###
    #############################################
    
    #We have determined additional corrections to individual sites
    #mostly to CO2 and LWdown records
    
    #Check if any apply to this site
    
    site_fixes <- site_exceptions(site_codes[s], var_data, att_data, qc_val=new_qc)
    
    
    #Replace with new fixed data
    var_data <- site_fixes$var_data
    att_data <- site_fixes$att_data
    
    
    
    
    ##########################
    ### Adjust time period ###
    ##########################
    
    
    #Get years to process
    start_yr <- qc_info$Start_year[site_ind]
    
    end_yr   <- qc_info$End_year[site_ind]
    
    
    
    #If need to adjust
    if (start_yr > 1 | end_yr < 0) {
      
      #New start and end year
      new_start_year <- years[1] + start_yr -1 
      new_end_year   <- years[length(years)] - end_yr
      
      
      #Start and end indices
      start_ind <- which(years == new_start_year)[1]
      end_ind   <- tail(which(years == new_end_year), 1)
      
      
      #New time vector
      time_var <- seq(0, by=60*60*24 / tsteps_per_day, length.out=length(c(start_ind:end_ind)))
      
      
      
      ### Adjust length of time-varying variables ##
      
      #Get dimensions for each variable
      dims <- lapply(vars, function(x) sapply(met_nc[[s]][["var"]][[x]][["dim"]], function(dim) dim[["name"]]))
      
      # #Find which variables are time-varying
      var_inds <- which(sapply(dims, function(x) any(x == "time")))
      

      
      #Change dimensions and values for time-varying data
      for (v in vars[var_inds]) {
        
        #Change time dimension
        met_nc[[s]]$var[[v]]$varsize[3] <- length(time_var)
        
        #Change time values
        met_nc[[s]]$var[[v]]$dim[[3]]$vals <- time_var
        
        #Change time size
        met_nc[[s]]$var[[v]]$dim[[3]]$len <- length(time_var)
        
        #Change length
        met_nc[[s]]$var[[v]]$size[3] <- length(time_var)
        
        #Change values in var_data
        var_data[[v]] <- var_data[[v]][start_ind:end_ind]
        
        # #Change chunk size (no idea what this is but produces an error otherwise
        # #during nc_create)
        # met_nc[[s]]$var[[v]]$chunksizes <- NA
      }
      
      
      #Also adjust time dimension and units
      
      #Change time dimensions
      #Change values, length and unit
      met_nc[[s]]$dim$time$vals  <- time_var
      met_nc[[s]]$dim$time$len   <- length(time_var)
      
      met_nc[[s]]$dim$time$units <- paste0("seconds since ", new_start_year, "-01-01 00:00:00")
      
    }
    
    
    
    ##########################################
    ### Check for missing vals in met data ###
    ##########################################
    
    
    #Do a final check to make sure there are no missing values in
    #any met variables
    
    #Loop through variables
    for (v in vars[var_inds]) {
    
      #Check if any missing values
      if (any (is.na(var_data[[v]]))) {
        
        stop(paste0("Missing values in ", v, ", site: ", site_codes[s]))
      } 
    }
    
    
    
    #And check that LAI data is available
    lai_vars <- which(grepl("LAI", vars))
    
    #Should have two available, check that they are there
    
    if (length(lai_vars) != 2) {
      stop(paste0("LAI variables not available, check site: ", site_codes[s]))
    }
    
    
    
    
    ##################### 
    ### Re-write file ###
    #####################
    
    
    
    
    
    
    
    #Add new QC flag value to metadata    
    
    
    
    
    
    
    
    
    
    
    #########------ Flux corrections ------#########
    
    
    
    #################################
    ### Energy balance correction ###
    #################################
    
    
    
    
    
    
    
    
    
    
    ##########################
    ### Adjust time period ###
    ##########################
    
    
    
  

  } #if processing 

} #sites




                         
                         