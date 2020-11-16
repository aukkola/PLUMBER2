library(ncdf4)
library(maptools)
library(viridis)


#clear R environment
rm(list=ls(all=TRUE))


path <- "/srv/ccrc/data04/z3509830/Fluxnet_data/All_flux_sites_processed_PLUMBER2/"

#Source functions
source(paste0(path, "/scripts/functions/cut_results.R"))


#Load world shapefile

world <- readShapePoly(paste0(path, "/../../World_shapefile/World"))



#--Define function ------

#Get length of time period
get_tperiod <- function(time) {
  
  len <- length(time)
  
  floor(len / (86400 / (  time[2] - time[1])) / 365)
  
}
#-----------------



############################
### Get site information ###
############################


#List site files
site_files <- list.files(paste0(path, "/Post-processed_PLUMBER2_outputs/Nc_files/Met/"), 
                         full.names=TRUE)

#Open file handles
nc <- lapply(site_files, nc_open)


#Get coordinates
lat <- sapply(nc, ncvar_get, varid='latitude')
lon <- sapply(nc, ncvar_get, varid='longitude')


#Get veg type (a few sites don't have short veg code, set manually later)
veg_type <- sapply(nc, function(x) tryCatch(ncvar_get(x, varid='IGBP_veg_short'), error=function(e) NA))


#Get length of time period (in years)
time_period <- sapply(nc, function(x) get_tperiod(ncvar_get(x, "time")))

start_yr <- sapply(nc, function(x) as.numeric(substr(ncatt_get(x, "time", "units")$value, 15, 18)))

end_yr <- start_yr + time_period - 1


#Close file connections
lapply(nc, nc_close)




################
### Plotting ###
################


breaks <- seq(0.5, 26, by=5)


cols <- rev(viridis(length(breaks)-1)) #colorRampPalette(c("#c7e9b4", "#7fcdbb", "#41b6c4",
                   #"#1d91c0", "#225ea8", "#0c2c84"))


classes <- cut_results(time_period, breaks)

plot_col <- cols[classes]


plot(world, col="grey95", border="grey50")


#Selected sites coloured by the number of site years
points(lon, lat, col=plot_col, cex=0.5, pch=20)

#Excluded sites in grey
points()

legend("bottomleft", legend=c("1-5", "6-10", "11-15", "16-20", "21", "Excluded"),
       col=c(cols, "grey50"), pch=c(rep(20, length(breaks)-1), 18), bty="n")














