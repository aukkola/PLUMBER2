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



#Close file connections
lapply(nc, nc_close)



##############################
### Get original site data ###
##############################


#List site files
all_site_files <- list.files(paste0(path, "/all_sites_no_duplicates/Nc_files/Met/"), 
                             full.names=TRUE)

#Open file handles
nc_all <- lapply(all_site_files, nc_open)

site_code_all <- sapply(nc_all, function(x) ncatt_get(x, varid=0, "site_code")$value)

included_ind <- which(site_code_all %in% site_code)

#Get length of time period (in years)
time_period_all <- sapply(nc_all, function(x) get_tperiod(ncvar_get(x, "time")))


#Get coordinates
lat_excluded <- sapply(nc_all[included_ind], ncvar_get, varid='latitude')
lon_excluded <- sapply(nc_all[included_ind], ncvar_get, varid='longitude')


#Close file connections
lapply(nc_all, nc_close)


#Check that site order matching before plotting
if(any(!(site_code == site_code_all[included_ind]))) stop("Site order different")

#Calculate number of excluded years
excluded_yrs <- time_period_all[included_ind] - time_period



################
### Plotting ###
################


#Set output directory
outdir <- paste0(path, "/PLUMBER2_paper_figs/")
dir.create(outdir)

#Output file
outfile <- paste0(outdir, "/Map_of_excluded_site_yrs.png")

png(outfile, height=3, width=8, res=500, unit="in")



layout(mat=matrix(c(1,2), ncol=2, byrow=TRUE),
       widths=c(0.7,0.3))
#layout.show(2)

par(mai=c(0,0,0,0.5))
par(omi=c(0.2,0,0,0.1))


land_col <- "grey97"


### World map ###


breaks <- c(-0.05, 0.5, 3.5, 6.5, 9.5, 13)


cols <- colorRampPalette(c("#c7e9b4", "#41b6c4", "#225ea8",
                           "#081d58"))(length(breaks)-1)
# cols <- rev(viridis(length(breaks)-1)) #colorRampPalette(c("#c7e9b4", "#7fcdbb", "#41b6c4",
# #"#1d91c0", "#225ea8", "#0c2c84"))



classes <- cut_results(excluded_yrs, breaks[2:(length(breaks)-1)])

plot_col <- cols[classes]


plot(crop(world, extent(c(-180, 180, -55, 85))), col=land_col, border="grey50",
     ylim=c(-55, 85))


#Selected sites coloured by the number of site years
points(lon, lat, col=plot_col, cex=0.55, pch=20)


legend(x=-180, y=15, legend=c("0-1", "2-3", "4-5", "6-7", "8-12"),
       col=cols, pch=c(rep(20, length(breaks)-1), 18), bty="n", cex=0.7,
       title="Excluded years", xpd=NA)

#panel number
mtext(side=3, "a)", cex=1.1, font=2, line=-2, adj=0, xpd=NA)



#Reset mai
par(mai=c(0.6,0.2,0.6,0))


### Histogram ###


#Data length (get frequencies from histogram)
hist_data <- hist(excluded_yrs, breaks=seq(-0.5, by=1, length.out=14), plot=FALSE)

bars <- barplot(height=hist_data$counts,  ylab="",
                xlab="", ylim=c(0, max(hist_data$counts)), col="#3690c0", las=2)

#Add x-axis
axis(side=1, at=bars[,1], 
     labels=hist_data$breaks[1:(length(hist_data$breaks)-1)]+0.5, las=2)

#y- and x-label
mtext(side=2, "Number of sites", line=2.5, cex=0.9)
mtext(side=1, "Excluded years", line=2, cex=0.9)

#panel number
mtext(side=3, "b)", cex=1.1, font=2, line=1, adj=-0.15, xpd=NA)



dev.off()


# To work out stats for regions and datasets:
#
# ### America, Australia and Europe insets ###
# 
# #Outlines on map
# 
# #America
# xmin_am <- -127
# xmax_am <- -67
# ymin_am <- 20
# ymax_am <- 60
# 
# #Europe
# xmin_eu <- -15
# xmax_eu <- 45
# ymin_eu <- 30
# ymax_eu <- 70
# 
# #Australia
# xmin_au <- 100
# xmax_au <- 160
# ymin_au <- -46
# ymax_au <- -6
# 
# #Collate coordinates
# lims <- list(America=c(xmin_am, xmax_am, ymin_am, ymax_am),
#              Europe=c(xmin_eu, xmax_eu, ymin_eu, ymax_eu),
#              Australia=c(xmin_au, xmax_au, ymin_au, ymax_au))
# 
# for (l in 1:length(lims)){
#   ind <- which(lon_excluded >= lims[[l]][1] & lon_excluded <= lims[[l]][2] &
#                lat_excluded >= lims[[l]][3] & lat_excluded <= lims[[l]][4]) #need this or R plots points outside domain
# 
#   print(paste0("l = ", l))
#   print(median(excluded_yrs[ind]))
#   print(round(mean(excluded_yrs[ind]),digits=1))
#   
#   
# }
# 
# 
# datasets <- c("FLUXNET2015", "OzFlux", "LaThuile")
# for (d in 1:length(datasets)) {
#   
#   ind_sites <- which(grepl(datasets[d], site_files))
#   
#   print(datasets[d])
#   print(median(excluded_yrs[ind_sites]))
#   print(round(mean(excluded_yrs[ind_sites]),digits=1))
#   
# }
# 
# 
# 
# 
# 
