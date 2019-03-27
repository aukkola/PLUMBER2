#clear R environment
rm(list=ls(all=TRUE))

path <- "/srv/ccrc/data45/z3509830/CABLE_runs/CABLE_site_runs/"


#List experiments
exp <- list.files(paste0(path, "/Outputs/"))


#List files 
exp_files <- lapply(exp, function(x) list.files(paste0(path, "/Outputs/", x, "/outputs/")))


#Find common files
common_sites <- intersect(exp_files[[1]], exp_files[[2]])


#Copy these to a new directory
new_dir <- paste0(path, "/CABLE_outputs_subset/")

dir.create(new_dir)


for (k in 1:length(exp)) {
  
  #Create new directory
  dir.create(paste(new_dir, exp[k], sep="/"))
  
  #Copy files
  lapply(common_sites, function(x) file.copy(from=paste0(path, "/Outputs/", exp[k], "/outputs/", x),
                                             to=paste(new_dir, exp[k], x, sep="/")))
  
  
}


#Then copy obs
obs_dir <- "/srv/ccrc/data04/z3509830/Fluxnet_data/All_flux_sites_processed/all_sites_no_duplicates/Nc_files/Flux"

patterns <- substr(common_sites, 1, 16)


#Create new directory
dir.create(paste0(new_dir, "/Observations/"))

#Copy files
lapply(patterns, function(x) file.copy(from=list.files(obs_dir, pattern=x, full.names=TRUE),
                                           to=paste0(new_dir, "/Observations/")))







