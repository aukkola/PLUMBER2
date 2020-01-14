#Need to run on Hurricane !! R version too old on other servers

library(MODISTools)
library(ncdf4)


#clear R environment
rm(list=ls(all=TRUE))

#Set path
path <- "/srv/ccrc/data04/z3509830/Fluxnet_data/All_flux_sites_processed_PLUMBER2/"


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

if (any(rm_ind)) {
  
  site_codes <- site_codes[-rm_ind]
  site_lat   <- site_lat[-rm_ind]
  site_lon   <- site_lon[-rm_ind]
  
}

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

print("Batch processing LAI")

#Get LAI
mt_batch_subset(df=sites_to_fetch, product=product, band=band_lai, start = "2000-01-01",
                end = format(Sys.time(), "%Y-%m-%d"), km_lr = km, km_ab = km,
                out_dir = outdir, ncores = 10, internal=FALSE)


print("Batch processing SD")

#Get LAI SD
mt_batch_subset(df=sites_to_fetch, product=product, band=band_sd, start = "2000-01-01",
                end = format(Sys.time(), "%Y-%m-%d"), km_lr = km, km_ab = km,
                out_dir = outdir, ncores = 10, internal=FALSE)


print("Batch processing QC")

#Get LAI QC
mt_batch_subset(df=sites_to_fetch, product=product, band=band_qc, start = "2000-01-01",
                end = format(Sys.time(), "%Y-%m-%d"), km_lr = km, km_ab = km,
                out_dir = outdir, ncores = 10, internal=FALSE)



print("Batch processing finished")

#For some reason, batch processing is not retrieving all sites
#Get these separately 


sites_fetched_lai <- sapply(list.files(outdir, pattern="_Lai_500m_"),
                        function(x) strsplit(x, "_")[[1]][1])

sites_fetched_sd <- sapply(list.files(outdir, pattern="_LaiStdDev_500m_"),
                            function(x) strsplit(x, "_")[[1]][1])


sites_fetched_qc <- sapply(list.files(outdir, pattern="_FparLai_QC_"),
                            function(x) strsplit(x, "_")[[1]][1])


#Find missing sites
missing_ind_lai <- which(!(sites_to_fetch$site_name %in% sites_fetched_lai))
missing_ind_sd  <- which(!(sites_to_fetch$site_name %in% sites_fetched_sd))
missing_ind_qc  <- which(!(sites_to_fetch$site_name %in% sites_fetched_qc))


#Get missing LAI
if (length(missing_ind_lai) > 0) {
  for (s in missing_ind_lai) {

    print(paste0("LAI: processing site ", sites_to_fetch$site_name[s]))
    
    #Get coordinates
    lat=sites_to_fetch$lat[s]
    lon=sites_to_fetch$lon[s]
    
    #Get LAI
    mt_subset(product = product, lat = lat, lon = lon,
              band = band_lai, start = "2000-01-01",
              end = format(Sys.time(), "%Y-%m-%d"),
              km_lr = km, km_ab = km,
              site_name = as.character(sites_to_fetch$site_name[s]),
              out_dir = outdir, internal=FALSE)
    
  }
}

    
#Get missing SD
if (length(missing_ind_sd) > 0) {
  for (s in missing_ind_sd) {
    
    print(paste0("SD: processing site ", sites_to_fetch$site_name[s]))
    
    #Get coordinates
    lat=sites_to_fetch$lat[s]
    lon=sites_to_fetch$lon[s]
    
    #Get LAI SD
    mt_subset(product = product, lat = lat, lon = lon,
              band = band_sd, start = "2000-01-01",
              end = format(Sys.time(), "%Y-%m-%d"),
              km_lr = km, km_ab = km,
              site_name = as.character(sites_to_fetch$site_name[s]),
              out_dir = outdir, internal=FALSE)
    
  }
}


#Get missing QC
if (length(missing_ind_qc) > 0) {
  for (s in missing_ind_qc) {
    
    print(paste0("QC: processing site ", sites_to_fetch$site_name[s]))
    
    #Get coordinates
    lat=sites_to_fetch$lat[s]
    lon=sites_to_fetch$lon[s]
    
    #Get LAI QC
    mt_subset(product = product, lat = lat, lon = lon,
              band = band_qc, start = "2000-01-01",
              end = format(Sys.time(), "%Y-%m-%d"),
              km_lr = km, km_ab = km,
              site_name = as.character(sites_to_fetch$site_name[s]),
              out_dir = outdir, internal=FALSE)
    
  }
}


    
    
    
    
    



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

# lat=74.4733
# lon=-20.5503
# lat=-14.063300  
# lon=131.31810
# subset <- mt_subset(product = "MCD15A2H",
#                     lat = lat,
#                     lon = lon,
#                     band = "Lai_500m",
#                     start = "2004-01-01",
#                     end = "2006-01-01",
#                     km_lr = 1,
#                     km_ab = 1,
#                     site_name = "testsite")


# 
# mt_subset(product = product, lat = lat, lon = lon,
#           band = band_lai, start = "2000-01-01",
#           end = "2003-01-01",
#           km_lr = km, km_ab = km,
#           site_name = as.character(sites_to_fetch$site_name[s]),
#           out_dir = outdir, internal=FALSE)
# 



