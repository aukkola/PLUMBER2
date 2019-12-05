#Adding two LAI estimates to site met files

#1) MODIS LAI (product MCD15A2H)
#2) Copernicus Global Land Service LAI (http://land.copernicus.eu/global/products/lai)

#Using monthly climatology for time periods where
#no time-varying LAI available

library(raster)
library(ncdf4)

#clear R environment
rm(list=ls(all=TRUE))


#Set path
path <- "/srv/ccrc/data04/z3509830/Fluxnet_data/All_flux_sites_processed_PLUMBER2"


#List site files
site_files <- list.files(paste0(path, "/all_sites_no_duplicates/Nc_files/Met/"),
                         full.names=TRUE)


#Open file handles
site_nc <- lapply(site_files, nc_open, write=TRUE)

#Get site codes
site_codes <- sapply(site_nc, function(x) ncatt_get(x, varid=0, "site_code")$value)
  
#Site latitude and longitude
site_lat <- sapply(site_nc, ncvar_get, "latitude")
site_lon <- sapply(site_nc, ncvar_get, "longitude")


#Loop through sites
for (s in 1:length(site_nc)) {
  
  
  #################
  ### MODIS LAI ###
  #################
  
  modis_file <- list.files(paste0(path, "/MODIS_LAI_time_series"))
  
  
  
  
  
  
  ######################
  ### Copernicus LAI ###
  ######################
  
  #Copernicus data was downloaded and processed by Arden Burrell
  
  #Set path
  cop_path <- "/srv/ccrc/data51/z3466821/LAI/M0040186/processed/Monthly"
  
  
  
  
}






