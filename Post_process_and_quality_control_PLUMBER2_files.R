#Run this code after processing data through FluxnetLSM and adding 
#LAI time series from MODIS and Copernicus


library(ncdf4)
library(gsheet)
library(parallel)
library(chron)

#clear R environment
rm(list=ls(all=TRUE))

#Set path
path <- "/srv/ccrc/data04/z3509830/Fluxnet_data/All_flux_sites_processed_PLUMBER2/"


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


#Flux data

#Get sites
flux_files <- list.files(paste0(path, "/all_sites_no_duplicates/Nc_files/Flux/"), full.names=TRUE)

#Open file handles
#flux_nc <- lapply(flux_files, nc_open, write=TRUE)



#Get site codes
site_codes <- sapply(met_nc, function(x) ncatt_get(x, varid=0, "site_code")$value)

#Close met file handles
lapply(met_nc, nc_close)



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

if (any(!(site_codes %in% qc_info$Site_code))) {
  stop("Site missing from google sheet")
}


status <- unlist(sapply(site_codes, function(x) qc_info$decision[which(qc_info$Site_code == x)]))
  




#Then find sites to process (i.e. decision is not kill)
good_sites <- which(status != "kill")

# #Good site names
# qc_sites_good <- qc_sites[good_sites]
# 
# good_inds <- sapply(qc_sites_good, function(x) which(site_codes == x))
# 

#Turn into a list so can be passed to functions better
#qc_info_list <- lapply(good_sites, function(x) qc_info[x,])


qc_info_list <- lapply(site_codes, function(x) qc_info[which(qc_info$Site_code ==x),])



#Initialise parallel cores
cl <- makeCluster(getOption('cl.cores', 12))

clusterExport(cl, 'qc_info_list')
clusterExport(cl, 'new_qc')
clusterExport(cl, 'global_co2')
clusterExport(cl, 'outdir_plot')

clusterExport(cl, 'met_corrections')
clusterExport(cl, 'flux_corrections')
clusterExport(cl, 'energy_balance_correction')
clusterExport(cl, 'hrs')
clusterExport(cl, 'add_months')
clusterExport(cl, 'SynthesizeLWdown')
clusterExport(cl, 'site_exceptions')
clusterExport(cl, 'VPD2RelHum')
clusterExport(cl, 'calc_esat')
clusterExport(cl, 'linear_pred_co2')
clusterExport(cl, 'diagnostic_plot')
clusterExport(cl, 'Timeseries')

clusterEvalQ(cl, library(ncdf4))
clusterEvalQ(cl, library(chron))

  

#######################
### Met corrections ###
#######################


#Met corrections

clusterMap(cl, function(met, out, qc) met_corrections(infile_met=met, outfile_met=out,
                                                      qc_info=qc, new_qc=new_qc, global_co2=global_co2),
           met=met_files[good_sites], out=outfiles_met[good_sites],
           qc=qc_info_list[good_sites])

# 
# 
# # For testing individual site:
# s=69
# met_corrections(infile_met=met_files[s], outfile_met=outfiles_met[s],
#                  qc=qc_info_list[[which(qc_sites %in% site_codes[s])]],
#                 new_qc=new_qc, global_co2=global_co2)
# 
# 
# mapply(function(met, out, qc) met_corrections(infile_met=met, outfile_met=out,
#                                                       qc_info=qc, new_qc=new_qc, global_co2=global_co2),
#            met=met_files[good_sites], out=outfiles_met[good_sites],
#            qc=qc_info_list[good_sites])


########################
### Flux corrections ###
########################


#Flux corrections
clusterMap(cl, function(flx, out, qc) flux_corrections(infile_flux=flx, outfile_flux=out,
                                                      qc_info=qc, new_qc=new_qc),
           flx=flux_files[good_sites], out=outfiles_flux[good_sites],
           qc=qc_info_list[which(qc_sites %in% site_codes[good_sites])])






#############################################
#########------ Plot new data ------#########
#############################################
    

clusterMap(cl, function(site, met_file, flux_file) diagnostic_plot(site_code=site, outdir=outdir_plot,
                                                                   met_file=met_file, flux_file=flux_file),
           site=site_codes[good_sites], met_file=met_files[good_sites], flux_file=flux_files[good_sites])





stopCluster(cl)   
  



  




                         
                         