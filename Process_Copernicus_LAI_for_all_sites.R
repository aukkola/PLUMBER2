library(ncdf4)
library(raster)

#clear R environment
rm(list=ls(all=TRUE))

#Set path
path <- "/srv/ccrc/data04/z3509830/Fluxnet_data/All_flux_sites_processed/"

#Create output folder for met files with processed LAI time series
outdir <- paste0(path, "/all_sites_no_duplicates/Nc_files/Met_with_LAI/")

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

lai <- brick(lapply(lai_files, brick))


#Loop through sites
for (s in 1:length(site_codes)) {
  
  
  #Set time stamps
  lai_time <- seq.Date(from=as.Date("1999-01-01"), by="month",
                       length.out=nlayers(lai))
  
  
  
  #Get site coordinates
  lat <- ncvar_get(site_nc[[s]], "latitude")
  lon <- ncvar_get(site_nc[[s]], "longitude")
  
  coords <- matrix(c(lon, lat), ncol=2)
  
  
  #Get time series data
  
  site_lai <- extract(lai, y=coords, )
  
  
  
  
  
  
  
  
  
  
  
  ######################################
  ### Add LAI time series to NC file ###
  ######################################
  
  #Save to file
  
  # Define variable:
  laivar = ncvar_def('LAI_Copernicus', '-', list(site_nc[[s]]$dim[[1]], site_nc[[s]]$dim[[2]], site_nc[[s]]$dim[[3]]),
                     missval=-9999,longname='Copernicus Global Land Service leaf area index')
  # Add variable and then variable data:
  site_nc[[s]] = ncvar_add(site_nc[[s]], laivar)
  ncvar_put(site_nc[[s]], 'LAI_MODIS', modis_tseries)
  
  #Close file handle
  nc_close(site_nc[[s]])
  
}






