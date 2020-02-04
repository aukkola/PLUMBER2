library(ncdf4)
library(raster)
library(zoo)

#clear R environment
rm(list=ls(all=TRUE))

#Set path (where site LAI data lives)
path <- "/srv/ccrc/data04/z3509830/Fluxnet_data/All_flux_sites_processed/"


#Which fluxnet data to use?
flx_path <- "/srv/ccrc/data04/z3509830/Fluxnet_data/All_flux_sites_processed_PLUMBER2/"


#----- function
#find nearest non-NA grid cell
sample_raster_NA <- function(r, xy)
{
        apply(X = xy, MARGIN = 1, 
        FUN = function(xy) r[which.min(replace(distanceFromPoints(r, xy), is.na(r), NA))])
}
#-----  
  
  
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


#"NO-Blv" is so far north, no Copernicus data available


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
    year_lai_vals <- as.vector(extract(lai_averaged, coords[[s]]) /(1/9) )
    
    
    #If didn't find any, find nearest non-NA pixel
    if (any (is.na(year_lai_vals))) {
      
      #Crop to smaller size to speed up distance calculation
      temp_lai <- crop(lai_averaged, extent(c(coords[[s]][1]-5, coords[[s]][1]+5,
                                              coords[[s]][2]-5, coords[[s]][2]+5)))

      
      #lapply to layers
      year_lai_vals <- sapply(1:nlayers(temp_lai), function(l) sample_raster_NA(temp_lai[[l]], coords[[s]]))
      
      
    }
    
    site_tseries[[s]] <- append(site_tseries[[s]], year_lai_vals)  
      
    
  }  
  
}



#Check that no missing LAI values
if (any(sapply(site_tseries, function(x) any(is.na(x))))) {
  stop("Missing LAI values present")
}







##########################
### Loop through sites ###
##########################


for (s in 1:length(site_nc)) {
  
  
  #Set time stamps (do inside loop as gets adjusted below)
  lai_time <- seq.Date(from=as.Date("1999-01-01"), by="month",
                       length.out=length(site_tseries[[s]]))
  
  
  #Smooth LAI time series with spline (and cap negative values)
  smooth_lai_ts = tryCatch(smooth.spline(lai_time, site_tseries[[s]])$y,
                           error= function(e) NA)
  if (all(is.na(smooth_lai_ts))) { 
    warning(paste0("could not process site", site_codes[s], ", missing values in LAI"))
    next 
  }
  
  smooth_lai_ts[smooth_lai_ts < 0] <- 0
  
  
  
  ##########################
  ### Create climatology ###
  ##########################
  

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
  
  
  
  
  ###################################
  ### Calculate running anomalies ###
  ###################################
  
  #Initialise
  lai_clim_anomalies <- rep(NA, length(lai_time))
  
  #Repeat climatology for whole time series
  copernicus_clim_all <- rep_len(copernicus_clim, length(lai_time))
  
  
  #Calculate running mean anomaly (+/- 6 months either side of each time step)
  anomaly <- rollmean(smooth_lai_ts - copernicus_clim_all, k=12, fill=NA)
  
  #Add rolling mean anomaly to climatology
  lai_clim_anomalies <- copernicus_clim_all + anomaly  
  
  
  #Check if remaining NA values from missing time steps, gapfill if found
  if (any(is.na(lai_clim_anomalies))) {
    
    #Find missing values
    missing <- which(is.na(lai_clim_anomalies))
    
    #Repeat climatology for all years and gapfill time series
    clim_all_yrs <- rep(copernicus_clim, floor(length(lai_time)/no_tsteps))
    lai_clim_anomalies[missing] <- clim_all_yrs[missing]
    
  }
  
  
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
  
  
    
  #Add climatological values to smoothed lai time series
  
  #Overwrite original time series with new extended data
  
  #If start year earlier than Copernicus
  if (startyr < copernicus_startyr) {
    
    #LAI data
    lai_clim_anomalies <- append(rep(copernicus_clim, copernicus_startyr - startyr), 
                                 lai_clim_anomalies)
    
    #Time vector
    lai_time <- append(seq.Date(as.Date(paste0(startyr, "-01-01")), by="month", 
                       length.out=no_tsteps * (copernicus_startyr - startyr)),
                       lai_time)
                                
  }
  
  #If end year later than Copernicus
  if (endyr > copernicus_endyr) {
    
    #LAI data
    lai_clim_anomalies <- append(smooth_lai_ts,
                            rep(lai_clim_anomalies, endyr - copernicus_endyr))
    
    
    #Time vector
    lai_time <- append(lai_time,
                       seq.Date(as.Date(paste0(copernicus_endyr+1, "-01-01")), by="month", 
                       length.out=no_tsteps * (endyr - copernicus_endyr)))
    
  }

  
  
  #Find modis time step corresponding to site start time
  start_ind <- which(lai_time == paste0(startyr, "-01-01"))
  end_ind   <- tail(which(grepl(endyr, lai_time)), 1) #Last index of end year
  
  #Extract MODIS time steps matching site
  copernicus_ts_for_site   <- lai_clim_anomalies[start_ind:end_ind]
  copernicus_time_for_site <- lai_time[start_ind:end_ind]
  
  
  #Repeat Copernicus time series to create a time series matching site time step
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

  