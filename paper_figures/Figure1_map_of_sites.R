library(ncdf4)
library(maptools)
library(viridis)
library(raster)


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

###########################
### CRU temp and precip ###
###########################

#Use 1990-2017

ind <- 109:444
precip <- mean(brick(paste0("/srv/ccrc/data04/z3509830/Obs_precip_products/CRU_TS4.02/Data_1981_2017/",
                            "CRU_TS4.02_monthly_total_precipitation_1981_2017.nc"))[[ind]]) * 12

ind <- 1069:1404
temp <- mean(brick(paste0("/srv/ccrc/data04/z3509830/Obs_Tair_SWdown_products/CRU_TS4.02/Data_1901_2017/",
                          "Tair/CRU_TS4.02_monthly_mean_air_temperature_1901_2017.nc"))[[ind]])



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

#Site code
site_code <- sapply(nc, function(x) ncatt_get(x, varid=0, "site_code")$value)

#Get length of time period (in years)
time_period <- sapply(nc, function(x) get_tperiod(ncvar_get(x, "time")))

start_yr <- sapply(nc, function(x) as.numeric(substr(ncatt_get(x, "time", "units")$value, 15, 18)))

end_yr <- start_yr + time_period - 1


#Get Tair and precip
tair_site   <- lapply(nc, ncvar_get, "Tair")
precip_site <- lapply(nc, ncvar_get, "Precip")


#Close file connections
lapply(nc, nc_close)



##########################
### Get excluded sites ###
##########################


#List site files
all_site_files <- list.files(paste0(path, "/all_sites_no_duplicates/Nc_files/Met/"), 
                         full.names=TRUE)

#Open file handles
nc_all <- lapply(all_site_files, nc_open)

site_code_all <- sapply(nc_all, function(x) ncatt_get(x, varid=0, "site_code")$value)

excluded_ind <- which(!(site_code_all %in% site_code))


#Get coordinates
lat_excluded <- sapply(nc_all[excluded_ind], ncvar_get, varid='latitude')
lon_excluded <- sapply(nc_all[excluded_ind], ncvar_get, varid='longitude')



#Close file connections
lapply(nc_all, nc_close)



################
### Plotting ###
################


layout(mat=matrix(c(1,1,1,2,3,4,5,6,7), ncol=3, byrow=TRUE),
       heights=c(1,0.5,0.5))
#layout.show(7)

par(mai=c(0.2,0.2,0.2,0.2))
par(omi=c(0.2,0.2,0.2,0.2))


breaks <- seq(0.5, 26, by=5)


cols <- rev(viridis(length(breaks)-1)) #colorRampPalette(c("#c7e9b4", "#7fcdbb", "#41b6c4",
                   #"#1d91c0", "#225ea8", "#0c2c84"))


classes <- cut_results(time_period, breaks)

plot_col <- cols[classes]


plot(world, col="grey95", border="grey50")


#Selected sites coloured by the number of site years
points(lon, lat, col=plot_col, cex=0.5, pch=20)

#Excluded sites in grey
points(lon_excluded, lat_excluded, pch=18, cex=0.5, col="grey30")

legend("bottomleft", legend=c("1-5", "6-10", "11-15", "16-20", "21", "Excluded"),
       col=c(cols, "grey30"), pch=c(rep(20, length(breaks)-1), 18), bty="n")

#America, Australia and Europe insets

#Outlines on map

#America
xmin_am <- -130
xmax_am <- -70
ymin_am <- 20
ymax_am <- 60
polygon(c(xmin_am, xmax_am, xmax_am, xmin_am),
        c(ymin_am, ymin_am, ymax_am, ymax_am), col=NA, border="black")

#Europe
xmin_eu <- -13
xmax_eu <- 45
ymin_eu <- 30
ymax_eu <- 70

polygon(c(xmin_eu, xmax_eu, xmax_eu, xmin_eu),
        c(ymin_eu, ymin_eu, ymax_eu, ymax_eu), col=NA, border="black")

#Australia
xmin_au <- 110
xmax_au <- 155
ymin_au <- -46
ymax_au <- -6

polygon(c(xmin_au, xmax_au, xmax_au, xmin_au),
        c(ymin_au, ymin_au, ymax_au, ymax_au), col=NA, border="black")


#Subplots

#America
plot(crop(world, extent(c(xmin_am, xmax_am, ymin_am, ymax_am))), col="grey95", border="grey50")
points(lon, lat, col=plot_col, cex=1, pch=20)
polygon(c(xmin_am, xmax_am, xmax_am, xmin_am),
        c(ymin_am, ymin_am, ymax_am, ymax_am), col=NA, border="black")

#Europe
plot(crop(world, extent(c(xmin_eu, xmax_eu, ymin_eu, ymax_eu))), col="grey95", border="grey50")
points(lon, lat, col=plot_col, cex=1, pch=20)
polygon(c(xmin_eu, xmax_eu, xmax_eu, xmin_eu),
        c(ymin_eu, ymin_eu, ymax_eu, ymax_eu), col=NA, border="black")

#Australia
plot(crop(world, extent(c(xmin_au, xmax_au, ymin_au, ymax_au))), col="grey95", border="grey50")
points(lon, lat, col=plot_col, cex=1, pch=20)
polygon(c(xmin_au, xmax_au, xmax_au, xmin_au),
        c(ymin_au, ymin_au, ymax_au, ymax_au), col=NA, border="black")




### Record length ###

#TODO: change into barplot !!!!!!!!!!!!

#Data length (histogram)
hist(time_period, xlab="Length of time period", ylab="Number of sites",
     main="")


### Vegetation type ###

#Biome type (histogram)
veg_counts <- sort(sapply(unique(veg_type), function(x) length(which(veg_type == x))),
                   decreasing=TRUE)
veg_classes <- gsub(" ", "", names(veg_counts))
barplot(veg_counts, names.arg=veg_classes, ylab="Number of sites")


### Climate envelope ###

#TODO: ADD excluded sites !!!!!!!!!!!!!!!!!!

#Climate envelope (MAT, MAP scatter)

#CRU data
temp_vals <- values(temp)
pr_vals   <- values(precip)


#Site data
mean_tair_site  <- sapply(tair_site, function(x) mean(x) - 273.15)
mean_precip_site <- sapply(precip_site, function(x) mean(x) * 60*60*24*365)



plot(temp_vals, pr_vals, cex=0.1, col="grey70", pch=20,
     xlab="Mean annual  temperature", ylab="Mean annual precipitation (mm)")


points(mean_tair_site, mean_precip_site, col="black", cex=0.6, pch=20)


