#Run this code after processing data through FluxnetLSM and adding 
#LAI time series from MODIS and Copernicus


library(ncdf4)
library(gsheet)
library(maptools)

#clear R environment
rm(list=ls(all=TRUE))

#Set path
path <- "/srv/ccrc/data04/z3509830/Fluxnet_data/All_flux_sites_processed_PLUMBER2/"


#Shapefile
world_file <- paste(path, "/../../World_shapefile/World", sep="")
world <- readShapePoly(world_file)


#Source functions
source(paste0(path, "/scripts/functions/site_exceptions.R"))
source(paste0(path, "/scripts/functions/diagnostic_plot.R"))
source(paste0(path, "/scripts/functions/plot_timeseries.R"))
source(paste0(path, "/scripts/functions/energy_balance_correction.R"))
source(paste0(path, "/scripts/functions/met_corrections.R"))
source(paste0(path, "/scripts/functions/flux_corrections.R"))
source(paste0(path, "/scripts/functions/gapfill_postprocess.R"))




#Output path
outpath <- paste0(path, "/Post-processed_PLUMBER2_outputs/")

dir.create(outpath)



###########################
### Read QC information ###
###########################

#Read QC information from google sheet,
#e.g. start and end year and CO2 processing


qc_info <- gsheet2tbl("https://docs.google.com/spreadsheets/d/1bi9bbUpwzRycDJ16VTYGx8ApkDLV_IeSwEsdDJsX1SI/edit#gid=0")

qc_sites <- qc_info$Site_code



#####################
### Read CO2 data ###
#####################

#Use Mauna Loa annual CO2 record to gapfill CO2 records that are
#poor quality or missing


global_co2 <- read.table(paste0(path, "/Global_CO2_data/co2_annmean_mlo.txt"))


#Set new QC values for post-processing
#Use 101 as OzFlux using a lot of smaller values
new_qc <- 101



######################
### Load site data ###
######################


#Met data

#Get sites
met_files <- list.files(paste0(path, "/all_sites_no_duplicates/Nc_files/Met_with_LAI/"), full.names=TRUE)

#Open file handles
met_nc <- lapply(met_files, nc_open, write=TRUE)


#Get site codes
site_codes <- sapply(met_nc, function(x) ncatt_get(x, varid=0, "site_code")$value)

lat <- sapply(met_nc, function(x) ncvar_get(x, "latitude"))
lon <- sapply(met_nc, function(x) ncvar_get(x, "longitude"))


#Close met file handles
lapply(met_nc, nc_close)


#status
status <- unlist(sapply(site_codes, function(x) qc_info$decision[which(qc_info$Site_code == x)]))


#Then find sites to process (i.e. decision is not kill)
good_sites <- which(status != "kill")


qc_info_list <- lapply(site_codes, function(x) qc_info[which(qc_info$Site_code ==x),])


#Map
plot(world, col="grey90")


#Loop through sites
for (s in good_sites) {
  
  points(lon[s], lat[s], pch=20, cex=0.8, col="red")
}


