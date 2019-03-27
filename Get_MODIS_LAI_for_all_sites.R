#Need to run on Hurricane !! R version too old on other servers

library(MODISTools)
library(ncdf4)


#clear R environment
rm(list=ls(all=TRUE))

#Set path
path <- "/srv/ccrc/data04/z3509830/Fluxnet_data/All_flux_sites_processed/"


#Get sites
site_files <- list.files(paste0(path, "/all_sites_no_duplicates/Nc_files/Met/"),
                         full.names=TRUE)

#Open file handles
site_nc <- lapply(site_files, nc_open)

#Get site codes
site_codes <- sapply(site_nc, function(x) ncatt_get(x, varid=0, "site_code")$value)

#Get latitude
site_lat <- sapply(site_nc, ncvar_get, "latitude")

#Get longitude
site_lon <- sapply(site_nc, ncvar_get, "longitude")


#Remove duplicates
rm_ind <- duplicated(site_codes)

site_codes <- site_codes[-rm_ind]
site_lat   <- site_lat[-rm_ind]
site_lon   <- site_lon[-rm_ind]

#Collate into data frame
sites_to_fetch     <- data.frame("site_name" = site_codes)
sites_to_fetch$lat <- site_lat
sites_to_fetch$lon <- site_lon
  

#Create output folder for LAI time series
outdir <- paste0(path, "/MODIS_LAI_time_series/Raw_data/")

dir.create(outdir, recursive=TRUE)


#Check that site hasn't already been processed

#MODIS/Terra+Aqua Leaf Area Index/FPAR 8-Day L4
product <- "MCD15A2H"

#LAI band
band_lai <- "Lai_500m"

#LAI SD band
band_sd <- "LaiStdDev_500m"

#LAI QC band
band_qc <- "FparLai_QC"

#Cells to obtain around site
#n.b. tried using 0.5 km but this makes the code crash (???)
km <- 1


#Get LAI
mt_batch_subset(df=sites_to_fetch, product=product, band=band_lai, start = "2000-01-01",
                end = format(Sys.time(), "%Y-%m-%d"), km_lr = km, km_ab = km,
                out_dir = outdir, ncores = 10, internal=FALSE)


#Get LAI SD
mt_batch_subset(df=sites_to_fetch, product=product, band=band_sd, start = "2000-01-01",
                end = format(Sys.time(), "%Y-%m-%d"), km_lr = km, km_ab = km,
                out_dir = outdir, ncores = 10, internal=FALSE)


#Get LAI QC
mt_batch_subset(df=sites_to_fetch, product=product, band=band_qc, start = "2000-01-01",
                end = format(Sys.time(), "%Y-%m-%d"), km_lr = km, km_ab = km,
                out_dir = outdir, ncores = 10, internal=FALSE)



# 
# 
# subsets <- mt_batch_subset(df = df,
#                            product = "MOD11A2",
#                            band = "LST_Day_1km",
#                            internal = TRUE,
#                            start = "2004-01-01",
#                            end = "2004-03-31")
# 
# df=sites_to_fetch[1:5,]
# test=mt_batch_subset(df=df, product=product        , band=band_lai, start = "2000-01-01",
#                 end = format(Sys.time(), "%Y-%m-%d"), km_lr = 0.5, km_ab = 0.5,
#                 out_dir = outdir, ncores = 10)
# 
# 
# #
# 
# subset <- mt_subset(product = "MOD13Q1",
#                     lat = 42.534171,
#                     lon = -72.179003,
#                     band = "250m_16_days_NDVI",
#                     start = "2004-01-01",
#                     end = "2005-12-30",
#                     km_lr = 1,
#                     km_ab = 1,
#                     site_name = "testsite")
