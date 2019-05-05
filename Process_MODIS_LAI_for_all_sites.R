library(ncdf4)


#clear R environment
rm(list=ls(all=TRUE))

#Set path
path <- "/srv/ccrc/data04/z3509830/Fluxnet_data/All_flux_sites_processed/"


#Create output folder for met files with processed LAI time series
outdir <- paste0(path, "/all_sites_no_duplicates/Nc_files/Met_with_LAI/")

dir.create(outdir, recursive=TRUE)


#Copy Met files to the new folder
file.copy(from = list.files(paste0(path, "/all_sites_no_duplicates/Nc_files/Met/"), full.names=TRUE), 
          to = outdir, overwrite=TRUE)


#Get sites
site_files <- list.files(outdir, full.names=TRUE)

#Open file handles
site_nc <- lapply(site_files, nc_open, write=TRUE)

#Get site codes
site_codes <- sapply(site_nc, function(x) ncatt_get(x, varid=0, "site_code")$value)




#MODISTools R package only lets you to extract 1km around site coordinates.
#This is greater than the average site footprint.
#Code returns data for 5 x 5 pixels, manually verified than the pixel numbers
#provided (1-25) go by rows, i.e.

# [ 1 ] [ 2 ]  [ 3 ] [ 4 ]  [ 5 ]
# [ 6 ] [ 7*]  [ 8*] [ 9*]  [10 ]
# [11 ] [12*]  [ce*] [14*]  [15 ]
# [16 ] [17*]  [18*] [19*]  [20 ]
# [21 ] [22 ]  [23 ] [24 ]  [25 ]

# Extracting pixels in the centre and immediately around it (*)
# These correspond to a radius of 500m around site coordinates

pixel_no <- c(7, 8, 9, 12, 13, 14, 17, 18, 19)

#Use only good quality data
#Random bit integer format, ask Martin if need to work these out again...

qc_flags <- c(0, 2, 24 ,26, 32, 34, 56, 58)



#Loop through sites
for (s in 1:length(site_codes)) {
  
  modis_path <- paste0(path, "/MODIS_LAI_time_series/Raw_data")
  
  lai_file <- list.files(modis_path, pattern=paste0(site_codes[s], "_MCD15A2H_Lai_500m_"), 
                         full.names=TRUE)
  
  qc_file <- list.files(modis_path, pattern=paste0(site_codes[s], "_MCD15A2H_FparLai_QC"), 
                         full.names=TRUE)
  
  sd_file <- list.files(modis_path, pattern=paste0(site_codes[s], "_MCD15A2H_LaiStdDev_500m_"), 
                         full.names=TRUE)
  
  #Check if found site file
  if (length(lai_file) == 0) {
    print(paste0("Could not find site: ", site_codes[s]))
    next
  }
  
  #Check that only found one file per variable
  if (any (c(length(lai_file), length(qc_file), length(sd_file)) > 1) ) {
    stop("Found too many lai files, check")
    
  }
  
  #Read data
  lai <- read.csv(lai_file, header=TRUE)
  qc  <- read.csv(qc_file, header=TRUE)
  sd  <- read.csv(sd_file, header=TRUE)
  
  #First mask out bad LAI values
  lai_masked <- lai 
  lai_masked$value[!(qc$value %in% qc_flags)] <- NA
  
  #Number of time steps
  no_tsteps <- dim(lai)[1] / max(lai$pixel)
  
  
  #Extract 3 x3 pixels
  lai_pixel <- matrix(nrow=no_tsteps, ncol=length(pixel_no))
  sd_pixel  <- matrix(nrow=no_tsteps, ncol=length(pixel_no))
  qc_pixel  <- matrix(nrow=no_tsteps, ncol=length(pixel_no))

  #Save time stamps
  lai_time <- as.Date(lai$calendar_date[which(lai$pixel == pixel_no[1])])
  
  #Loop through pixels
  for (p in 1:length(pixel_no)) {
    
    #Get time series for pixel and scale using scale factor
    lai_pixel[,p] <- lai$value[which(lai$pixel == pixel_no[p])] * lai$scale[1]
    sd_pixel[,p]  <- sd$value[which(sd$pixel == pixel_no[p])] * sd$scale[1]
    qc_pixel[,p]  <- qc$value[which(sd$pixel == pixel_no[p])] 

  }
  
  
  ##################################
  ### Mask out poor quality data ###
  ##################################
  
  #Mask out where QC flag 
  lai_pixel <- replace(lai_pixel, !(qc_pixel %in% qc_flags), NA)
  sd_pixel  <- replace(sd_pixel, !(qc_pixel %in% qc_flags), NA)
  
  #Also mask out where sd really low (likely cloud effects)
  sd_pixel  <- replace(sd_pixel, sd_pixel < 0.1, NA)
  lai_pixel <- replace(lai_pixel, is.na(sd_pixel), NA)
  
  ### Average each time step ###
  
  #Initialise lai time series
  lai_ts <- vector(length=no_tsteps)
  #lai_ts_mean <- vector(length=no_tsteps)
  
  #Loop through time steps
  for (t in 1:no_tsteps) {
    
    #If no values available
    if (all(is.na(lai_pixel[t,]))) {
      
      lai_ts[t] <- NA
    
    #If values available
    } else {
      
      #Weight grid cell estimates by their standard deviation
      #Following Martin's method (https://github.com/mdekauwe/get_MODIS_LAI_australia/blob/master/build_modis_climatology.py),
      #but normalising by sum of standard deviations
      
      sd_vals <- sd_pixel[t,]
        
      weights = (1/sd_vals**2) / sum(1/sd_vals**2, na.rm=TRUE)
      
      #Check that weights sum up to 1 (because of a precision issue presumably,
      #rounding to 5 decimals, otherwise might not equal 1 even when correct)
      if (round(sum(weights, na.rm=TRUE), 5) != 1) stop("Weighting not correct")
      
      #Calculate weighted average
      lai_ts[t] <- weighted.mean(lai_pixel[t,], w=weights, na.rm=TRUE)
      
      #lai_ts_mean[t] <- mean(lai_pixel[t,], na.rm=T)
    }
  }
  

  ######################################
  ### Gapfill and smooth with spline ###
  ######################################
  
  #Set x values
  x <- 1:length(lai_ts)
  
  #Define spline function
  func = splinefun(x=x, y=lai_ts, method="fmm",  ties = mean)
  
  #Gapfill with spline (and cap negative values)
  lai_spline <- func(seq(min(x), max(x), by=1))
  lai_spline[lai_spline < 0] <- 0
  
  #Smooth with spline (and cap negative values)
  smooth_lai_ts = smooth.spline(x, lai_spline)$y
  smooth_lai_ts[smooth_lai_ts < 0] <- 0
  
  
  # #Test
  # plot(lai_spline, type='l')
  # lines(lai_ts, col='red')
  # lines(smooth_lai_ts, col='blue')

  
  

  
  #######################################
  ### Add missing time steps in MODIS ###
  #######################################
  
  #Some MODIS time series are missing time steps, not sure why.
  #Add these missing time steps by using the mean of previous and next value
  
  
  #Get modis timing information
  modis_startyr <- as.numeric(format(lai_time[1], "%Y"))
  modis_endyr   <- as.numeric(format(lai_time[length(lai_time)], "%Y"))
  
  
  #Check if LAI has missing time steps
  if (any (diff.Date(lai_time) > 8)) {

    #Loop through years
    all_tsteps <- vector()
    for (y in modis_startyr:modis_endyr){
      
      all_tsteps <- append(all_tsteps, seq.Date(as.Date(paste0(y, "-01-01")), 
                           by=lai_time[2]-lai_time[1], length.out=46))
      
    }
    
    #Remove extra tstesp for final year
    all_tsteps <- all_tsteps[-which(all_tsteps > lai_time[length(lai_time)])]
    
    #Find missing tsteps
    missing <- which(!(all_tsteps %in% lai_time))

    
    #Create new time series, where missing tsteps set to NA
    new_ts <- smooth_lai_ts
    
    for (n in missing) { 
      
      #Add new value
      new_ts <- append(new_ts, NA, after=n-1) 
      
      #Then gapfill missing values by taking mean of previous and next values
      
      #!!!! N.B. this still returns a warning, need to fix !!!!!!!!!!!!!!!!!!!!!!!!!!
      
      previous <- new_ts[max(which(!is.na(new_ts[1:max(1, (n-1))])))] #make sure can't go negative
      nextval <- new_ts[n+min(which(!is.na(new_ts[min(length(new_ts), n+1):length(new_ts)])))] #make sure can't exceed vec length
      
      new_ts[n] <- mean(previous, nextval, na.rm=TRUE)
      }
    
    
    
    
    #Replace time vector
    lai_time <- all_tsteps
    smooth_lai_ts <- new_ts

  }
  
  
  
  
  ###################################
  ### Match with site time series ###
  ###################################
  
  
  #Each year has 46 time steps, but first and last year are incomplete
  
  #Get timing info for site
  site_start_time <- ncatt_get(site_nc[[s]], "time")$units
  site_time       <- ncvar_get(site_nc[[s]], "time")
  site_tstep_size <- 86400 / (site_time[2] - site_time[1]) 
  
  #Extract year
  startyr    <- as.numeric(substr(site_start_time, start=15, stop=18))
  obs_length <- length(site_time)
  nyr        <- round(obs_length/(site_tstep_size*365))
  endyr      <- startyr + nyr - 1 
  
  
  #Need to create climatological average if site time series starts before modis
  #(use 2003 as modis start year to test this as modis starts halfway through 2002)
  if (startyr < modis_startyr+1) {
    
    ### Construct modis climatology ###
    
    #First find time steps for first (incomplete) year
    first_year <- which(grepl(modis_startyr, lai_time))
    
    #Then find time steps for last (incomplete) year
    last_year <- which(grepl(modis_endyr, lai_time))
    
    #MODIS has 46 time steps per year
    no_tsteps <- length(which(grepl(modis_startyr +1 , lai_time)))

    #Initialise
    modis_clim <- vector(length=no_tsteps)
    
    for (c in 1:no_tsteps) {
      
      #Indices for whole years
      inds <- seq(c + length(first_year), by=no_tsteps, length.out=floor(length(lai_time)/no_tsteps))
        
      #Add first and last year if applicable
      if (c > (no_tsteps - length(first_year))) { inds <- append(inds, first_year[c-length(first_year)])}
      
      if( c <= length(last_year)) { inds <- append(inds, last_year[c]) }
      
      #Calculate average for time step
      modis_clim[c] <- mean(smooth_lai_ts[inds])
      
    }
    
    #Add climatological values to smoothed lai time series
    
    #Overwrite original time series with new extended data
    smooth_lai_ts <- append(c(rep(modis_clim, modis_startyr - startyr), 
                              modis_clim[1:(no_tsteps - length(first_year))]),
                              smooth_lai_ts)
    
    
    #Add new time stamps
    extended_lai_time <- vector()
    for (y in startyr:modis_startyr) {
      
      extended_lai_time <- append(extended_lai_time, seq.Date(as.Date(paste0(y, "-01-01")), 
               by=lai_time[2]-lai_time[1], length.out=no_tsteps))
      
    }
    
    #Overwrite lai_time with new time series
    lai_time <- append(extended_lai_time[1:(length(extended_lai_time) - length(first_year))],
                                lai_time)
    
  }
  
  
  #Find modis time step corresponding to site start time
  start_ind <- which(lai_time == paste0(startyr, "-01-01"))
  end_ind   <- tail(which(grepl(endyr, lai_time)), 1) #Last index of end year
  
  #Extract MODIS time steps matching site
  modis_ts_for_site   <- smooth_lai_ts[start_ind:end_ind]
  modis_time_for_site <- lai_time[start_ind:end_ind]
  
  
  #Repeat modis time series to create a time series matching site time step
  modis_tseries <- vector()
  

  #Loop through time steps
  for (t in 1:length(modis_time_for_site)) {
    
    #Last time step
    if (t == length(modis_time_for_site)) {
    
      #Use the number of time steps that ensures final time series matches the length of site data
      modis_tseries <- append(modis_tseries, rep(modis_ts_for_site[t], length(site_time) - length(modis_tseries)))
      
    #All other time steps  
    } else {
      
      time_diff <- modis_time_for_site[t+1] - modis_time_for_site[t]
      
      #Repeat each days estimate by the number of days and time steps per day
      modis_tseries <- append(modis_tseries, rep(modis_ts_for_site[t], time_diff * site_tstep_size))
      
    }
    
  }
  
  
  #Check that the number of time steps match
  if (length(modis_tseries) != length(site_time)) stop("MODIS and site time steps don't match")
  
  
  
  ######################################
  ### Add LAI time series to NC file ###
  ######################################
  
  #Save to file
  
  # Define variable:
  laivar = ncvar_def('LAI_MODIS', '-', list(site_nc[[s]]$dim[[1]], site_nc[[s]]$dim[[2]], site_nc[[s]]$dim[[3]]),
                          missval=-9999,longname='MODIS 8-daily LAI')
  # Add variable and then variable data:
  site_nc[[s]] = ncvar_add(site_nc[[s]], laivar)
  ncvar_put(site_nc[[s]], 'LAI_MODIS', modis_tseries)
  
  #Close file handle
  nc_close(site_nc[[s]])
  
  
}





# 
# # 
# # x <- 1:length(lai_ts)
# # func = splinefun(x=x, y=lai_ts, method="fmm",  ties = mean)
# # lai_spline <- func(seq(min(x), max(x), 1))
# # 
# # test = smooth.spline(x, lai_spline)
# # 
# # 
# # 
# # 
# # lai_block <- matrix(lai_spline, ncol=2, byrow=T)
# # lai_block <- apply(lai_block, MARGIN=1, max)
# # x1 <- x[seq(1, by=2, length.out=length(lai_block))]
# # test1 = smooth.spline(x1, lai_block, spar=0.1)
# # test1 = smooth.spline(x1, lai_block)
# # 
# # 
# # ind <- c(1:400)
# 
# plot(lai_spline[ind], type='l', col="red", xlab="Time", ylab='LAI')
# #lines(lai_ts_mean, col='blue')
# lines(lai_ts[ind], col='blue')
# #lines(lai_ts[1:200], col='black')
# #lines(test$y, col="green")
# 
# lines(x1[ind], test1$y[ind], col="green", lwd=2)
# 
# legend("top", legend=c("original", "gapfilled with splinefun", "smoothed with spline"), fill=c("blue", "red", "green"))
# 
# 


