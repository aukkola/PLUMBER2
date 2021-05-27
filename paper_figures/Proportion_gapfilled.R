library(ncdf4)
library(maptools)
library(viridis)
library(raster)


#clear R environment
rm(list=ls(all=TRUE))


path <- "/srv/ccrc/data04/z3509830/Fluxnet_data/All_flux_sites_processed_PLUMBER2/"



#--Define function ------

#Get length of time period
get_tperiod <- function(time) {
  
  len <- length(time)
  
  floor(len / (86400 / (  time[2] - time[1])) / 365)
  
}
#-----------------


######################
### Screened sites ###
######################


#List site files
site_files <- list.files(paste0(path, "/Post-processed_PLUMBER2_outputs/Nc_files/Met/"), 
                         full.names=TRUE)

#Open file handles
nc <- lapply(site_files, nc_open)


#Site code
site_code <- sapply(nc, function(x) ncatt_get(x, varid=0, "site_code")$value)

#Get length of time period (in years)
time_period <- sapply(nc, function(x) get_tperiod(ncvar_get(x, "time")))

start_yr <- sapply(nc, function(x) as.numeric(substr(ncatt_get(x, "time", "units")$value, 15, 18)))

end_yr <- start_yr + time_period - 1


#Get QC flags 
tair   <- lapply(nc, ncvar_get, "Tair_qc")
precip <- lapply(nc, function(x) tryCatch(ncvar_get(x, "Precip_qc"), error=function(e) NA))
swdown <- lapply(nc, ncvar_get, "SWdown_qc") 
qair   <- lapply(nc, function(x) tryCatch(ncvar_get(x, "Qair_qc"), error=function(e) NA))
wind   <- lapply(nc, ncvar_get, "Wind_qc") 
co2    <- lapply(nc, ncvar_get, "CO2air_qc") 
lwdown <- lapply(nc, ncvar_get, "LWdown_qc") 


#Close file connections
lapply(nc, nc_close)



##########################
### Non-screened sites ###
##########################



#List site files
all_site_files <- list.files(paste0(path, "/all_sites_no_duplicates/Nc_files/Met/"), 
                             full.names=TRUE)

#Open file handles
nc_all <- lapply(all_site_files, nc_open)

site_code_all <- sapply(nc_all, function(x) ncatt_get(x, varid=0, "site_code")$value)

included_ind <- which(site_code_all %in% site_code)


#Get QC flags 
tair_all   <- lapply(nc_all[included_ind], ncvar_get, "Tair_qc")
precip_all <- lapply(nc_all[included_ind], function(x) tryCatch(ncvar_get(x, "Precip_qc"), error=function(e) NA))
swdown_all <- lapply(nc_all[included_ind], ncvar_get, "SWdown_qc") 
qair_all   <- lapply(nc_all[included_ind], function(x) tryCatch(ncvar_get(x, "Qair_qc"), error=function(e) NA))
wind_all   <- lapply(nc_all[included_ind], ncvar_get, "Wind_qc") 
co2_all    <- lapply(nc_all[included_ind], ncvar_get, "CO2air_qc") 
lwdown_all <- lapply(nc_all[included_ind], ncvar_get, "LWdown_qc") 


#Close file connections
lapply(nc_all, nc_close)



##############################
### Calculate % gap-filled ###
##############################


#Level of gap-filling in screened data

tair_screened   <- sapply(tair, function(x) length(which(x > 0)) / length(x) * 100)
precip_screened <- sapply(precip, function(x) length(which(x > 0)) / length(x) * 100)
swdown_screened <- sapply(swdown, function(x) length(which(x > 0)) / length(x) * 100)
qair_screened   <- sapply(qair, function(x) length(which(x > 0)) / length(x) * 100)
wind_screened   <- sapply(wind, function(x) length(which(x > 0)) / length(x) * 100)
lwdown_screened <- sapply(lwdown, function(x) length(which(x > 0)) / length(x) * 100)
co2_screened    <- sapply(co2, function(x) length(which(x > 0)) / length(x) * 100)



#Level of gap-filling in non_screened data

tair_non_screened   <- sapply(tair_all, function(x) length(which(x > 0)) / length(x) * 100)
precip_non_screened <- sapply(precip_all, function(x) length(which(x > 0)) / length(x) * 100)
swdown_non_screened <- sapply(swdown_all, function(x) length(which(x > 0)) / length(x) * 100)
qair_non_screened   <- sapply(qair_all, function(x) length(which(x > 0)) / length(x) * 100)
wind_non_screened   <- sapply(wind_all, function(x) length(which(x > 0)) / length(x) * 100)
lwdown_non_screened <- sapply(lwdown_all, function(x) length(which(x > 0)) / length(x) * 100)
co2_non_screened    <- sapply(co2_all, function(x) length(which(x > 0)) / length(x) * 100)

co2_missing_nonscrn <- sapply(co2_all, function(x) length(which(is.na(x))) / length(x) * 100)


### Key variables ###
# 
# #Tair
round(length(which(unlist(tair) > 0)) / length(unlist(tair)) * 100, digits=1)
round(length(which(unlist(tair_all) > 0)) / length(unlist(tair_all)) * 100, digits=1)

# #Precip
round(length(which(unlist(precip) > 0)) / length(unlist(precip)) * 100, digits=1)
round(length(which(unlist(precip_all) > 0)) / length(unlist(precip_all)) * 100, digits=1)

#SWdown
round(length(which(unlist(swdown) > 0)) / length(unlist(swdown)) * 100, digits=1)
round(length(which(unlist(swdown_all) > 0)) / length(unlist(swdown_all)) * 100, digits=1)

#Qair
round(length(which(unlist(qair) > 0)) / length(unlist(qair)) * 100, digits=1)
round(length(which(unlist(qair_all) > 0)) / length(unlist(qair_all)) * 100, digits=1)

#Wind
round(length(which(unlist(wind) > 0)) / length(unlist(wind)) * 100, digits=1)
round(length(which(unlist(wind_all) > 0)) / length(unlist(wind_all)) * 100, digits=1)

#LWdown
round(length(which(unlist(lwdown) > 0)) / length(unlist(lwdown)) * 100, digits=1)
round(length(which(unlist(lwdown_all) > 0)) / length(unlist(lwdown_all)) * 100, digits=1)


#CO2
round(length(which(unlist(co2) > 0)) / length(unlist(co2)) * 100, digits=1)
round(length(which(unlist(co2_all) > 0)) / length(unlist(co2_all)) * 100, digits=1)


#CO2 missing in original data
round(mean(co2_missing_nonscrn), digits=1)




### All met variables ###

all_screened <- list(tair, precip, swdown, qair, wind, co2, lwdown)

all_non_screened <- list(tair_all, precip_all, swdown_all, 
                         qair_all, wind_all, co2_all, lwdown_all)



#Gapfilling
screened_gapfilled     <- length(which(unlist(all_screened) > 0)) / length(unlist(all_screened)) * 100

non_screened_gapfilled <- length(which(unlist(all_non_screened) > 0)) / length(unlist(all_non_screened)) * 100





