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



################
### Get data ###
################


#List site files
site_files <- list.files(paste0(path, "/Post-processed_PLUMBER2_outputs/Nc_files/Flux/"), 
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


#Get original and corrected fluxes

qle_cor <- lapply(nc, function(x) tryCatch(ncvar_get(x, "Qle_cor"), error=function(e) NA))
qh_cor  <- lapply(nc, function(x) tryCatch(ncvar_get(x, "Qh_cor"), error=function(e) NA))
  
ind <- which(sapply(qle_cor, length) > 1) #sites with corrected fluxes

qle_cor <- qle_cor[ind]
qh_cor  <- qh_cor[ind]

#Only get original data for sites for which corrected data is available   

qle <- lapply(nc[ind], ncvar_get, "Qle")
qh  <- lapply(nc[ind], ncvar_get, "Qh")

rnet <- lapply(nc[ind], ncvar_get, "Rnet")
qg  <- lapply(nc[ind], function(x) tryCatch(ncvar_get(x, "Qg"),  error=function(e) NA))


ebc <- lapply(1:length(qle), function(x) if (length(qg[[x]]) > 1) (qle[[x]] + qh[[x]]) / (rnet[[x]] + qg[[x]]) else NA)


#Number of sites corrected
length(ind)
       
#Average difference in corrected and original data across all sites
mean(unlist(qle_cor) / unlist(qle), na.rm=TRUE)

#Qh has some 0 values leading to Inf when dividing by it, fix 
diff <- unlist(qh_cor) / unlist(qh)

mean(diff[which(!is.infinite(diff))], na.rm=TRUE)




#Average difference in corrected and original data across individual sites

#Qle
ratio_qle <- lapply(1:length(qle_cor), function(x) qle_cor[[x]]/ qle[[x]])

#Remove Inf (a couple of sites return NA because all Qle_cor values are missing)
sites_qle <- sapply(ratio_qle,  function(x) mean(x[which(!is.infinite(x))], na.rm=TRUE))

#Qh
ratio_qh <- lapply(1:length(qh_cor), function(x) qh_cor[[x]]/ qh[[x]])

#Remove Inf (a couple of sites return NA because all Qle_cor values are missing)
sites_qh <- sapply(ratio_qh,  function(x) mean(x[which(!is.infinite(x))], na.rm=TRUE))




#Proportion of missing data
#Qle
length(which(is.na(unlist(qle_cor)))) / length(unlist(qle_cor)) * 100 #all sites
  
length(which(is.na(unlist(qle)))) / length(unlist(qle)) * 100

sapply(qle_cor, function(x) length(which(is.na(x))) / length(x) * 100) #individual sites

sapply(qle, function(x) length(which(is.na(x))) / length(x) * 100)


#Qh
length(which(is.na(unlist(qh_cor)))) / length(unlist(qh_cor)) * 100

length(which(is.na(unlist(qh)))) / length(unlist(qh)) * 100






#Variance

var_qle_cor <- sapply(qle_cor, function(x) sd(x, na.rm=TRUE) / mean(x, na.rm=TRUE))

var_qle <- sapply(qle, function(x) sd(x, na.rm=TRUE) / mean(x, na.rm=TRUE))



#Average ECB error
ebc_all <- unlist(ebc)

mean(ebc_all[which(!is.infinite(ebc_all))], na.rm=TRUE)
median(ebc_all[which(!is.infinite(ebc_all))], na.rm=TRUE)


mean(sapply(ebc, function(x) mean(x[which(!is.infinite(x))], na.rm=TRUE)), na.rm=TRUE)




breaks <- c(0, 0.5, 0.75, 1, 1.25, 1.5, 3)


cols <- colorRampPalette(c("#c7e9b4", "#41b6c4", "#225ea8",
                           "#081d58"))(length(breaks)-1)
# cols <- rev(viridis(length(breaks)-1)) #colorRampPalette(c("#c7e9b4", "#7fcdbb", "#41b6c4",
# #"#1d91c0", "#225ea8", "#0c2c84"))

land_col <- "grey97"

classes <- cut_results(sites_qle, breaks[2:(length(breaks)-1)])

plot_col <- cols[classes]


plot(crop(world, extent(c(-180, 180, -55, 85))), col=land_col, border="grey50",
     ylim=c(-55, 85))


#Selected sites coloured by the number of site years
points(lon[ind], lat[ind], col=plot_col, cex=0.55, pch=20)





#Level of gap-filling in all data
#Qle
qc <- lapply(nc, ncvar_get, "Qle_qc")

range(sapply(qc, function(x) length(which(x > 0)) / length(x) * 100))

length(which(unlist(qc) > 0)) / length(unlist(qc)) * 100


#Qh
qc <- lapply(nc, ncvar_get, "Qh_qc")

range(sapply(qc, function(x) length(which(x > 0)) / length(x) * 100))

length(which(unlist(qc) > 0)) / length(unlist(qc)) * 100


#NEE
qc <- lapply(nc, ncvar_get, "NEE_qc")

range(sapply(qc, function(x) length(which(x > 0)) / length(x) * 100))

length(which(unlist(qc) > 0)) / length(unlist(qc)) * 100


