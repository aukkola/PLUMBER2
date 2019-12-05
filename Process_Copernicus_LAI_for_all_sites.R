library(ncdf4)
library(raster)

#clear R environment
rm(list=ls(all=TRUE))

#Set path (where site LAI data lives)
path <- "/srv/ccrc/data04/z3509830/Fluxnet_data/All_flux_sites_processed/"


#Which fluxnet data to use?
flx_path <- "/srv/ccrc/data04/z3509830/Fluxnet_data/All_flux_sites_processed_PLUMBER2/"


#Create output folder for met files with processed LAI time series
outdir <- paste0(flx_path, "/all_sites_no_duplicates/Nc_files/Met_with_LAI/")

#Get sites
site_files <- list.files(outdir, full.names=TRUE)

#Open file handles
site_nc <- lapply(site_files, nc_open, write=TRUE)

#Get site codes
site_codes <- sapply(site_nc, function(x) ncatt_get(x, varid=0, "site_code")$value)


### Get LAI data ###

#Set LAI path
lai_path <- "/srv/ccrc/data51/z3466821/LAI/M0040186/processed/Monthly"

#Get files
lai_files <- list.files(lai_path, pattern="LAI_monthlymax_", full.names=TRUE)

#Should span 1999-2017
if(length(lai_files) != 19) stop("incorrect LAI files found")



#Set copernicus start and end year
copernicus_startyr <- 1999
copernicus_endyr   <- 2017


lai_outdir <- paste0(path, "/Copernicus_LAI_time_series/")

#Get site values for each time slice at a time
#(very slow to read in all the LAI data in one go)

#Initialise
site_tseries <- list()


#Get site coordinates
lat <- lapply(site_nc, ncvar_get, varid="latitude")
lon <- lapply(site_nc, ncvar_get, varid="longitude")

coords <- mapply(function(lon, lat) matrix(c(lon, lat), ncol=2),
                 lon=lon, lat=lat, SIMPLIFY=FALSE)



#No data is found for Cape Tribulation
#and the longitude in the Ozflux data file vs. their website doesn't match
#Using longitude on the website here so data is retrieved
if (length(which(site_codes == "AU-Ctr") > 0)){
  au_ctr <- which(site_codes == "AU-Ctr")
  coords[[au_ctr]][1] <- 145.3778
}
#Also no data at DK-Lva (checked that coordinates do match La Thuile metadata file)
#Adjust so finds data
if (length(which(site_codes == "DK-Lva") > 0)){
  dk_lva <- which(site_codes == "DK-Lva")
  coords[[dk_lva]][1] <- 12.1213
}


#These sites also need fixing...
#"NO-Blv" "IT-SRo" "IT-Pia" "IT-Noe" "IT-Cp2" "IT-Cpz" "DK-ZaH" "ES-ES1"
#Need to look for a function that finds nearest non-NA cell



#Loop through LAI time slices
for (y in 1:length(lai_files)) {
  
  #LAI output file (save as slow to process)
  outfile_lai <- paste0(lai_outdir, "/Copernicus_monthly_LAI_focal_mean_", 
                        y+copernicus_startyr-1, ".nc")
  
  
  #Check if file already processed
  if (file.exists(outfile_lai)) {
    
    lai_averaged <- brick(outfile_lai, varname="LAI")
    
  #Else process and save to file  
  } else {
    
    #Read LAI time slice
    lai <- brick(lai_files[y])
    
    #Average LAI so that each pixel represents the mean of itself 
    #and the neigbouring pixels
    lai_averaged <- brick(lapply(1:nlayers(lai), function(x) 
                          focal(lai[[x]], w=matrix(1,3,3), fun=mean)))
    
    writeRaster(lai_averaged, outfile_lai, varname="LAI", overwrite=TRUE)
    
  }
  
  
  
  #Extract site value
  
  for (s in 1:length(site_nc)) {
    
    #Initialise
    if(y == 1) site_tseries[[s]] <- vector()
    
    #Get time series !NB NEED TO FIX HERE TO REMOVE SCALING  1/9 (caused by error in focal function) !!!!!!!!!!!!!!!!!!!!!!!
    site_tseries[[s]] <- append(site_tseries[[s]], as.vector(extract(lai_averaged, coords[[s]]) /(1/9) ))
  }  
  
}




##########################
### Loop through sites ###
##########################


for (s in 1:length(site_nc)) {
  
  
  #Set time stamps (do inside loop as gets adjusted below)
  lai_time <- seq.Date(from=as.Date("1999-01-01"), by="month",
                       length.out=length(site_tseries[[1]]))
  
  
  #Smooth LAI time series with spline (and cap negative values)
  smooth_lai_ts = tryCatch(smooth.spline(lai_time, site_tseries[[s]])$y,
                           error= function(e) NA)
  if (all(is.na(smooth_lai_ts))) { 
    warning(paste0("could not process site", site_codes[s], ", missing values in LAI"))
    next 
  }
  
  smooth_lai_ts[smooth_lai_ts < 0] <- 0
  
  
  ###################################
  ### Match with site time series ###
  ###################################

  
  #Each year has 12 time steps
  
  #Get timing info for site
  site_start_time <- ncatt_get(site_nc[[s]], "time")$units
  site_time       <- ncvar_get(site_nc[[s]], "time")
  site_tstep_size <- 86400 / (site_time[2] - site_time[1]) 
  
  #Extract year
  startyr    <- as.numeric(substr(site_start_time, start=15, stop=18))
  obs_length <- length(site_time)
  nyr        <- round(obs_length/(site_tstep_size*365))
  endyr      <- startyr + nyr - 1 
  
  
  #Need to create climatological average if site time series starts before 
  #or ends after Copernicus data
  if (startyr < copernicus_startyr | endyr > copernicus_endyr) {
    
    
    ### Construct Copernicus climatology ###
    
    #Copernicus has 12 time steps per year
    no_tsteps <- length(which(grepl(copernicus_startyr, lai_time)))
    
    #Initialise
    copernicus_clim <- vector(length=no_tsteps)
    
    for (c in 1:no_tsteps) {
      
      #Indices for whole years
      inds <- seq(c, by=no_tsteps, length.out=length(lai_time)/no_tsteps)
      
      #Calculate average for time step
      copernicus_clim[c] <- mean(smooth_lai_ts[inds])
      
    }
    
    #Add climatological values to smoothed lai time series
    
    #Overwrite original time series with new extended data
    
    #If start year earlier than Copernicus
    if (startyr < copernicus_startyr) {
      
      #LAI data
      smooth_lai_ts <- append(rep(copernicus_clim, copernicus_startyr - startyr), 
                              smooth_lai_ts)
      
      #Time vector
      lai_time <- append(seq.Date(as.Date(paste0(startyr, "-01-01")), by="month", 
                         length.out=no_tsteps * (copernicus_startyr - startyr)),
                         lai_time)
                                  
    }
    
    #If end year later than Copernicus
    if (endyr > copernicus_endyr) {
      
      #LAI data
      smooth_lai_ts <- append(smooth_lai_ts,
                              rep(copernicus_clim, endyr - copernicus_endyr))
      
      
      #Time vector
      lai_time <- append(lai_time,
                         seq.Date(as.Date(paste0(copernicus_endyr+1, "-01-01")), by="month", 
                         length.out=no_tsteps * (endyr - copernicus_endyr)))
      
    }
  }
  
  
  #Find modis time step corresponding to site start time
  start_ind <- which(lai_time == paste0(startyr, "-01-01"))
  end_ind   <- tail(which(grepl(endyr, lai_time)), 1) #Last index of end year
  
  #Extract MODIS time steps matching site
  copernicus_ts_for_site   <- smooth_lai_ts[start_ind:end_ind]
  copernicus_time_for_site <- lai_time[start_ind:end_ind]
  
  
  #Repeat modis time series to create a time series matching site time step
  copernicus_tseries <- vector()
  
  
  #Loop through time steps
  for (t in 1:length(copernicus_time_for_site)) {
    
    #Last time step
    if (t == length(copernicus_time_for_site)) {
      
      #Use the number of time steps that ensures final time series matches the length of site data
      copernicus_tseries <- append(copernicus_tseries, rep(copernicus_ts_for_site[t], 
                                   length(site_time) - length(copernicus_tseries)))
      
      #All other time steps  
    } else {
      
      time_diff <- copernicus_time_for_site[t+1] - copernicus_time_for_site[t]
      
      #Repeat each days estimate by the number of days and time steps per day
      copernicus_tseries <- append(copernicus_tseries, rep(copernicus_ts_for_site[t], 
                                   as.numeric(time_diff * site_tstep_size)))
      
    }
    
  }
  
  
  #Check that the number of time steps match
  if (length(copernicus_tseries) != length(site_time)) stop("MODIS and site time steps don't match")
  
  
  
  
  ######################################
  ### Add LAI time series to NC file ###
  ######################################
  
  #Save to file
  
  # Define variable:
  laivar = ncvar_def('LAI_Copernicus', '-', list(site_nc[[s]]$dim[[1]], site_nc[[s]]$dim[[2]], site_nc[[s]]$dim[[3]]),
                     missval=-9999,longname='Copernicus Global Land Service leaf area index')
  # Add variable and then variable data:
  site_nc[[s]] = ncvar_add(site_nc[[s]], laivar)
  ncvar_put(site_nc[[s]], 'LAI_Copernicus', copernicus_tseries)
  
  #Close file handle
  nc_close(site_nc[[s]])
  

  
}

  