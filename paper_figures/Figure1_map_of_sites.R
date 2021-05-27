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
precip_site <- lapply(nc, ncvar_get, "Precip") ### CHANGE TO avPrecip variable !!!!!!


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

#Get length of time period (in years)
time_period_all <- sapply(nc_all, function(x) get_tperiod(ncvar_get(x, "time")))


#Get coordinates
lat_excluded <- sapply(nc_all[excluded_ind], ncvar_get, varid='latitude')
lon_excluded <- sapply(nc_all[excluded_ind], ncvar_get, varid='longitude')


#Get Tair and precip
tair_excluded_site   <- lapply(nc_all[excluded_ind], ncvar_get, "Tair")
precip_excluded_site <- lapply(nc_all[excluded_ind], ncvar_get, "Precip") ### CHANGE TO avPrecip variable !!!!!!


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


#Exclude colour
excl_col <- "#c51b7d"

land_col <- "grey97"


### World map ###


breaks <- c(seq(0.5, 20, by=5), 22)




cols <- colorRampPalette(c("#c7e9b4", "#41b6c4", "#225ea8",
                           "#081d58"))(length(breaks)-1)
  
  #colorRampPalette(c("#7fcdbb", "#41b6c4",
  #                        "#1d91c0", "#225ea8", "#0c2c84"))(length(breaks)-1) #rev(viridis(length(breaks)-1)) 

classes <- cut_results(time_period, breaks[2:(length(breaks)-1)])

plot_col <- cols[classes]


plot(crop(world, extent(c(-180, 180, -55, 85))), col=land_col, border="grey50",
     ylim=c(-55, 85))


#Excluded sites
points(lon_excluded, lat_excluded, pch=18, cex=0.9, col=excl_col)

#Selected sites coloured by the number of site years
points(lon, lat, col=plot_col, cex=0.75, pch=20)


legend(x=-180, y=0, legend=c("1-5", "6-10", "11-15", "16-21", "Excluded"),
       col=c(cols, excl_col), pch=c(rep(20, length(breaks)-1), 18), bty="n", cex=1.2,
       title="Site years")

#panel number
mtext(side=3, "a)", cex=1.1, font=2, line=-2, adj=0, xpd=NA)



### America, Australia and Europe insets ###

#Outlines on map

#America
xmin_am <- -127
xmax_am <- -67
ymin_am <- 20
ymax_am <- 60

#Europe
xmin_eu <- -15
xmax_eu <- 45
ymin_eu <- 30
ymax_eu <- 70

#Australia
xmin_au <- 100
xmax_au <- 160
ymin_au <- -46
ymax_au <- -6

#Collate coordinates
lims <- list(America=c(xmin_am, xmax_am, ymin_am, ymax_am),
             Europe=c(xmin_eu, xmax_eu, ymin_eu, ymax_eu),
             Australia=c(xmin_au, xmax_au, ymin_au, ymax_au))


#Polygons on map for subplots
for (l in 1:length(lims)) {
  polygon(c(lims[[l]][1], lims[[l]][2], lims[[l]][2], lims[[l]][1]),
          c(lims[[l]][3], lims[[l]][3], lims[[l]][4], lims[[l]][4]), col=NA, border="black")

  #label
  text(x=lims[[l]][1], y=lims[[l]][4], font=2, adj=c(-0.25,1.25), cex=1.1, paste0(l, ".")) 
  
}


#plot subplots
for (l in 1:length(lims)) {
  
  plot(crop(world, extent(lims[[l]])), col=land_col, 
       border="grey50", ylim=lims[[l]][3:4], xlim=lims[[l]][1:2],
       xaxs="i", yaxs="i")
  
  
  #Excluded sites
  ind <- which(lon_excluded >= lims[[l]][1] & lon_excluded <= lims[[l]][2] &
               lat_excluded >= lims[[l]][3] & lat_excluded <= lims[[l]][4]) #need this or R plots points outside domain
  
  points(lon_excluded[ind], lat_excluded[ind], pch=18, cex=0.9, col=excl_col)

  #Included sites
  ind <- which(lon >= lims[[l]][1] & lon <= lims[[l]][2] &
                 lat >= lims[[l]][3] & lat <= lims[[l]][4]) #need this or R plots points outside domain
  
  points(lon[ind], lat[ind], col=plot_col, cex=1.25, pch=20, xpd=FALSE)

    
  #Box around plot
  polygon(c(lims[[l]][1], lims[[l]][2], lims[[l]][2], lims[[l]][1]),
          c(lims[[l]][3], lims[[l]][3], lims[[l]][4], lims[[l]][4]), col=NA, border="black")
  
  #label
  mtext(side=3, adj=0.03, line=-4, font=2, cex=0.9, paste0(l, ".")) 
  
  
}


#Reset mai
par(mai=c(0.2,0.35,0.4,0.1))


### Record length ###

#Data length (get frequencies fromhistogram)
hist_data <- hist(time_period, plot=FALSE)

bars <- barplot(height=hist_data$counts,  ylab="",
        xlab="", ylim=c(0, 50), col="#3690c0", las=2)

#Add x-axis
axis(side=1, at=c(bars[,1]-0.6, bars[nrow(bars),1]+0.7), labels=hist_data$breaks, las=2)

#y- and x-label
mtext(side=2, "Number of sites", line=2.5, cex=0.9)
mtext(side=1, "Length of time period (years)", line=2.5, cex=0.9)

#panel number
mtext(side=3, "b)", cex=1.1, font=2, line=1, adj=-0.15, xpd=NA)



### Vegetation type ###

#Biome type (histogram)
veg_counts <- sort(sapply(unique(veg_type), function(x) length(which(veg_type == x))),
                   decreasing=TRUE)
veg_classes <- gsub(" ", "", names(veg_counts))
barplot(veg_counts, names.arg=veg_classes, ylab="", las=2, col="#3690c0")

#y-label
mtext(side=2, "Number of sites", line=2.5, cex=0.9)

#panel number
mtext(side=3, "c)", cex=1.1, font=2, line=1, adj=-0.15, xpd=NA)



### Climate envelope ###

#Climate envelope (MAT, MAP scatter)

#CRU data
temp_vals <- values(temp)
pr_vals   <- values(precip)


#Density colours
dens_cols  <- densCols(temp_vals, pr_vals, 
                       colramp=colorRampPalette(c("#d0d1e6", "#045a8d")),
                       nbin=1000)



#Site data

#Included
mean_tair_site   <- sapply(tair_site, function(x) mean(x) - 273.15)
mean_precip_site <- sapply(precip_site, function(x) mean(x) * 60*60*24*365)

#Excluded

mean_tair_excluded_site   <- sapply(tair_excluded_site, function(x) mean(x) - 273.15)
mean_precip_excluded_site <- sapply(precip_excluded_site, function(x) mean(x) * 60*60*24*365)



#CRU points
plot(temp_vals, pr_vals, cex=0.1, col=dens_cols, pch=20,
     xlab="", ylab="", ylim=c(0,8000))

#y- and x-label
mtext(side=2, "MAP (mm)", line=2.5, cex=0.9)
mtext(side=1, expression("MAT"~"("*degree*"C)"), line=2.5, cex=0.9)


#Excluded sites
points(mean_tair_excluded_site, mean_precip_excluded_site, col=excl_col, cex=0.6, pch=18)

#Site points
points(mean_tair_site, mean_precip_site, col="black", cex=0.6, pch=20)

#Legend
legend("topleft", c("Included", "Excluded"), col=c("black", excl_col), pch=c(20, 18), bty="n")

#panel number
mtext(side=3, "d)", cex=1.1, font=2, line=1, adj=-0.15, xpd=NA)



dev.off()
