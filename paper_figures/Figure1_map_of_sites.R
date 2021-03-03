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


#Set output directory
outdir <- paste0(path, "/PLUMBER2_paper_figs/")
dir.create(outdir)

#Output file
outfile <- paste0(outdir, "/Map_of_sites.png")

png(outfile, height=8, width=8, res=500, unit="in")



layout(mat=matrix(c(1,1,1,2,3,4,5,6,7), ncol=3, byrow=TRUE),
       heights=c(0.55,0.4,0.5))
#layout.show(7)

par(mai=c(0,0.2,0,0.2))
par(omi=c(0.4,0.2,0,0.1))


breaks <- seq(0.5, 26, by=5)


cols <- rev(viridis(length(breaks)-1)) #colorRampPalette(c("#c7e9b4", "#7fcdbb", "#41b6c4",
                   #"#1d91c0", "#225ea8", "#0c2c84"))


classes <- cut_results(time_period, breaks)

plot_col <- cols[classes]


plot(crop(world, extent(c(-180, 180, -55, 85))), col="grey95", border="grey50",
     ylim=c(-55, 85))


#Excluded sites in grey
points(lon_excluded, lat_excluded, pch=18, cex=0.9, col="black")

#Selected sites coloured by the number of site years
points(lon, lat, col=plot_col, cex=0.5, pch=20)


legend(x=-180, y=0, legend=c("1-5", "6-10", "11-15", "16-20", "21", "Excluded"),
       col=c(cols, "black"), pch=c(rep(20, length(breaks)-1), 18), bty="n", cex=1.2)

#panel number
mtext(side=3, "a)", cex=1.1, font=2, line=-2, adj=0, xpd=NA)



### America, Australia and Europe insets ###

#Outlines on map

#America
xmin_am <- -127
xmax_am <- -67
ymin_am <- 20
ymax_am <- 60
polygon(c(xmin_am, xmax_am, xmax_am, xmin_am),
        c(ymin_am, ymin_am, ymax_am, ymax_am), col=NA, border="black")

#Europe
xmin_eu <- -15
xmax_eu <- 45
ymin_eu <- 30
ymax_eu <- 70

polygon(c(xmin_eu, xmax_eu, xmax_eu, xmin_eu),
        c(ymin_eu, ymin_eu, ymax_eu, ymax_eu), col=NA, border="black")

#Australia
xmin_au <- 100
xmax_au <- 160
ymin_au <- -46
ymax_au <- -6

polygon(c(xmin_au, xmax_au, xmax_au, xmin_au),
        c(ymin_au, ymin_au, ymax_au, ymax_au), col=NA, border="black")


#Subplots

lims <- list(America=c(xmin_am, xmax_am, ymin_am, ymax_am),
             Europe=c(xmin_eu, xmax_eu, ymin_eu, ymax_eu),
             Australia=c(xmin_au, xmax_au, ymin_au, ymax_au))

#plot subplots
for (l in 1:length(lims)) {
  
  plot(crop(world, extent(lims[[l]])), col="grey95", 
       border="grey50", ylim=lims[[l]][3:4], xlim=lims[[l]][1:2],
       xaxs="i", yaxs="i")
  
  #Included sites
  points(lon, lat, col=plot_col, cex=1, pch=20, xpd=FALSE)
  
  #Excluded sites
  points(lon_excluded, lat_excluded, pch=18, cex=0.9, col="black")
  
  #Box around plot
  polygon(c(lims[[l]][1], lims[[l]][2], lims[[l]][2], lims[[l]][1]),
          c(lims[[l]][3], lims[[l]][3], lims[[l]][4], lims[[l]][4]), col=NA, border="black")
  
  #label
  mtext(side=3, adj=0.03, line=-4, font=2, cex=0.9, paste0(l, ".")) 
  
  
}


#Reset mai
par(mai=c(0.2,0.4,0.3,0.2))


### Record length ###

#TODO: change into barplot !!!!!!!!!!!!

#Data length (get frequencies fromhistogram)
hist_data <- hist(time_period, plot=FALSE)

bars <- barplot(height=hist_data$counts,  ylab="",
        xlab="", ylim=c(0, 50), col="#3690c0")

#Add x-axis
axis(side=1, at=c(bars[,1]-0.6, bars[nrow(bars),1]+0.7), labels=hist_data$breaks)

#y- and x-label
mtext(side=2, "Number of sites", line=2.5)
mtext(side=1, "Length of time period (years)", line=2.5)

#panel number
mtext(side=3, "b)", cex=1.1, font=2, line=1, adj=-0.15, xpd=NA)



### Vegetation type ###

#Biome type (histogram)
veg_counts <- sort(sapply(unique(veg_type), function(x) length(which(veg_type == x))),
                   decreasing=TRUE)
veg_classes <- gsub(" ", "", names(veg_counts))
barplot(veg_counts, names.arg=veg_classes, ylab="", las=2, col="#3690c0")

#y-label
mtext(side=2, "Number of sites", line=2.5)

#panel number
mtext(side=3, "c)", cex=1.1, font=2, line=1, adj=-0.15, xpd=NA)



### Climate envelope ###

#TODO: ADD excluded sites !!!!!!!!!!!!!!!!!!

#Climate envelope (MAT, MAP scatter)

#CRU data
temp_vals <- values(temp)
pr_vals   <- values(precip)


#Density colours
dens_cols  <- densCols(temp_vals, pr_vals, 
                       colramp=colorRampPalette(c("#d0d1e6", "#045a8d")),
                       nbin=1000)



#Site data
mean_tair_site  <- sapply(tair_site, function(x) mean(x) - 273.15)
mean_precip_site <- sapply(precip_site, function(x) mean(x) * 60*60*24*365)

#CRU points
plot(temp_vals, pr_vals, cex=0.1, col=dens_cols, pch=20,
     xlab="", ylab="")

#y- and x-label
mtext(side=2, "Mean annual precipitation (mm)", line=2.5)
mtext(side=1, "Mean annual  temperature", line=2.5)


#Excluded sites


#Site points
points(mean_tair_site, mean_precip_site, col="black", cex=0.6, pch=20)



#panel number
mtext(side=3, "d)", cex=1.1, font=2, line=1, adj=-0.15, xpd=NA)



dev.off()
