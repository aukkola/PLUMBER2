#Run this code after processing data through FluxnetLSM and adding 
#LAI time series from MODIS and Copernicus


library(ncdf4)
library(gsheet)
library(parallel)

path <- "/srv/ccrc/data04/z3509830/Fluxnet_data/All_flux_sites_processed_PLUMBER2/"


#Source function
source(paste0(path, "/scripts/functions/site_exceptions.R"))
source(paste0(path, "/scripts/functions/diagnostic_plot.R"))
source(paste0(path, "/scripts/functions/plot_timeseries.R"))
source(paste0(path, "/scripts/functions/energy_balance_correction.R"))
source(paste0(path, "/scripts/functions/met_corrections.R"))
source(paste0(path, "/scripts/functions/flux_corrections.R"))



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


#Flux data

#Get sites
flux_files <- list.files(paste0(path, "/all_sites_no_duplicates/Nc_files/Flux/"), full.names=TRUE)

#Open file handles
flux_nc <- lapply(flux_files, nc_open, write=TRUE)



#Get site codes
site_codes <- sapply(met_nc, function(x) ncatt_get(x, varid=0, "site_code")$value)



#Create output directories

outdir_met <- paste0(outpath, "/Nc_files/Met/")
dir.create(outdir_met, recursive=TRUE)

outdir_flux <- paste0(outpath, "/Nc_files/Flux/")
dir.create(outdir_flux, recursive=TRUE)

outdir_plot <- paste0(outpath, "/Diagnostic_plots/")
dir.create(outdir_plot, recursive=TRUE)



#Create output file names
#Met
outfiles_met <- sapply(met_files, function(x)paste0(outdir_met, "/", basename(x)))

#Flux
outfiles_flux <- sapply(flux_files, function(x)paste0(outdir_flux, "/", basename(x)))



### Check which sites to process ###

status <- unlist(sapply(site_codes, function(x) qc_info$decision[which(qc_sites == x)]))
  

#Check that all sites available
if (any (is.na(status))) {
  
  missing_sites <- site_codes[which(is.na(status))]
  
  stop("Sites missing in google doc")
}


#Then find site to process (i.e. decision is not kill)

good_sites <- which(status != "kill")

# 
# #For site exceptions, REMOVE AFTER
# site <- which(site_codes == "DK-Sor")
# met_nc <- met_nc[[site]]
# site_code <- site_codes[site]
# qc_info <- qc_info[site,]
# 
# 



#Initialise parallel cores

cl <- makeCluster(getOption('cl.cores', 12))

clusterExport(cl, 'good_sites')
clusterExport(cl, 'qc_info')
clusterExport(cl, 'new_qc')
clusterExport(cl, 'outfiles_met')
clusterExport(cl, 'outfiles_flux')

clusterExport(cl, 'met_corrections')
clusterExport(cl, 'flux_corrections')
clusterExport(cl, 'energy_balance_correction')

clusterEvalQ(cl, library(ncdf4))

  

#######################
### Met corrections ###
#######################



#Met corrections
# parLapply(cl, good_sites, function(x) met_corrections(met_nc=met_nc[[x]], outfile_met=outfiles_met[x], 
#                                                       site_code=site_codes[x], 
#                                                       qc_info=qc_info[which(qc_sites == site_codes[x]),], 
#                                                       new_qc=new_qc))


clusterMap(cl, function(met, out, qc) met_corrections(met_nc=met, outfile_met=out, 
                                                       qc_info=qc, new_qc=new_qc),
           met=met_nc[good_sites], out=outfiles_met[good_sites], 
           qc=qc_info[which(qc_sites %in% site_codes[good_sites]),])




mapply(function(met, out, qc) met_corrections(met_nc=met, outfile_met=out, 
                                                      qc_info=qc, new_qc=new_qc),
           met=met_nc[good_sites], out=outfiles_met[good_sites], 
           qc=qc_info[which(qc_sites %in% site_codes[good_sites]),])




# lapply(good_sites, function(x) met_corrections(met_nc=met_nc[[x]], outfile_met=outfiles_met[x],
#                                                       site_code=site_codes[x],
#                                                       qc_info=qc_info[which(qc_sites == site_codes[x]),],
#                                                       new_qc=new_qc))


########################
### Flux corrections ###
########################






# parLapply(cl, good_sites, function(x) flux_corrections(flux_nc=flux_nc[[x]], outfile_flux=outfiles_flux[x], 
#                                                       qc_info=qc_info[which(qc_sites == site_codes[x]),], 
#                                                       new_qc=new_qc))

clusterMap(cl, function(flx, out, qc) flux_corrections(flux_nc=flx, outfile_flux=out, 
                                            qc_info=qc, new_qc=new_qc),
           flx=flux_nc[good_sites], out=outfiles_flux[good_sites], 
           qc=qc_info[which(qc_sites %in% site_codes[good_sites]),])




#flux_nc <- lapply(flux_files, nc_open, write=TRUE)

# source(paste0(path, "/scripts/functions/flux_corrections.R"))
# 
lapply(good_sites, function(x) flux_corrections(flux_nc=flux_nc[[x]], outfile_flux=outfiles_flux[x],
                                                       qc_info=qc_info[which(qc_sites == site_codes[x]),],
                                                       new_qc=new_qc))

# 
# 
# lapply(flux_nc, nc_close)
# 

#############################################
#########------ Plot new data ------#########
#############################################
    

    diagnostic_plot(site_code=site_codes[s], outdir=outdir_plot,  
                    met_file=outfile_met, flux_file=outfile_flux)







stopCluster(cl)   
  



  




                         
                         