#devtools::install_github("aukkola/FluxnetLSM", ref="master") #Package broken for some reason, must be installed locally

setwd("/srv/ccrc/data04/z3509830/Fluxnet_data//All_flux_sites_processed_no_filtering/FluxnetLSM")
#not sure why no-lock option is needed but won't install otherwise
install.packages(".", repos=NULL, type='source', INSTALL_opts = c('--no-lock')) 


library(FluxnetLSM)  # convert_fluxnet_to_netcdf
library(parallel)
library(ncdf4)

#clear R environment
rm(list=ls(all=TRUE))


#Thresholds for missing and gap-filled time steps
missing_met  <- 0   #max. percent missing (must be set)
missing_flux <- 100 

gapfill_met_tier1 <- 100  #max. gapfilled percentage
gapfill_met_tier2 <- 100
gapfill_flux      <- 100

min_yrs      <- 1   #min. number of consecutive years


#Set main path
path <- "/srv/ccrc/data04/z3509830/Fluxnet_data/"

#Set output path for all data
out_path <- paste0(path, "/All_flux_sites_processed_no_filtering")

#Number of cluster
ncl <- 12

### Process all datasets separately and then find non-duplicate sites ###


####################
### FLUXNET 2015 ###
####################

#Outputs will be saved to this directory
out_path_flx <- paste0(out_path, "/FLUXNET2015_sites/") 

#Remove path
unlink(out_path_flx, recursive = TRUE)

### Hourly and Halfhourly ###
tstep <- c("Hourly", "Halfhourly")

#Loop through time steps
for(k in 1:length(tstep)){
  
  in_path  <- paste(path, "/FLUXNET2016/Original_data/",
                    tstep[k], "_qc_fixed", sep="")
  era_path <- paste(path, "/FLUXNET2016/Original_data/EraInterim/",
                    tstep[k], sep="")
  
  
  # Input Fluxnet data files (using FULLSET in this example, se R/Helpers.R for details)
  infiles <- get_fluxnet_files(in_path)
  
  #Retrieve site codes
  site_codes <- sapply(infiles, get_path_site_code)
  
  #Retrieve dataset versions
  datasetversions <- sapply(infiles, get_fluxnet_version_no)
  
  # Find ERA-files corresponding to site codes
  ERA_files     <- sapply(site_codes, function(x) get_fluxnet_erai_files(era_path, site_code=x))
  
  
  #Stop if didn't find ERA files
  if(any(sapply(ERA_files, length)==0)){
    stop("No ERA files found, amend input path")
  }
  
  
  ### Process files ###
  
  #Initialise clusters (using 2 cores here)
  cl <- makeCluster(getOption('cl.cores', ncl))
  
  #Import variables to cluster
  clusterExport(cl, "out_path_flx")
  clusterExport(cl, "convert_fluxnet_to_netcdf")
  if(exists("missing_met"))   {clusterExport(cl, "missing_met")}
  if(exists("missing_flux"))   {clusterExport(cl, "missing_flux")}
  if(exists("gapfill_met_tier1"))   {clusterExport(cl, "gapfill_met_tier1")}
  if(exists("gapfill_met_tier2"))   {clusterExport(cl, "gapfill_met_tier2")}
  if(exists("gapfill_flux"))   {clusterExport(cl, "gapfill_flux")}
  if(exists("min_yrs"))   {clusterExport(cl, "min_yrs")}

  
  #Loops through sites
  clusterMap(cl=cl, function(w,x,y,z) tryCatch(convert_fluxnet_to_netcdf(infile=w, site_code=x, out_path=out_path_flx,
                                                                         datasetversion=z, met_gapfill="ERAinterim", 
                                                                         flux_gapfill="statistical", era_file=y,
                                                                         missing_met=missing_met, missing_flux=missing_flux,
                                                                         gapfill_met_tier1=gapfill_met_tier1,
                                                                         gapfill_met_tier2=gapfill_met_tier2,
                                                                         gapfill_flux=gapfill_flux, min_yrs=min_yrs, 
                                                                         model="CABLE", check_range_action="warn",
                                                                         include_all_eval=TRUE),
                                               error = function(e) NULL),
             w=infiles, x=site_codes, y=ERA_files, z=datasetversions)
  
  stopCluster(cl)
  
  
}




#################
### La Thuile ###
#################


#Input path
in_path  <- paste(path, "/LaThuile/Original_data/raw_data", sep="/")

#Output path
out_path_lt <- paste0(out_path, "/LaThuile_sites/") 

#Remove path
unlink(out_path_lt, recursive = TRUE)

#Input Fluxnet data files (using FULLSET in this example, se R/Helpers.R for details)
infiles <- get_fluxnet_files(in_path, datasetname="LaThuile")


#Retrieve site codes
site_codes <- unique(sapply(infiles, get_path_site_code))


#Reorganise input files by site
infiles <- lapply(site_codes, function(site) get_fluxnet_files(in_path, site_code=site,
                                                               datasetname="LaThuile"))


### Process files ###

#Initialise clusters (using 2 cores here)
cl <- makeCluster(getOption('cl.cores', ncl))

#Import variables to cluster
clusterExport(cl, "out_path_lt")
clusterExport(cl, "convert_fluxnet_to_netcdf")
if(exists("missing_met"))   {clusterExport(cl, "missing_met")}
if(exists("missing_flux"))   {clusterExport(cl, "missing_flux")}
if(exists("gapfill_met_tier1"))   {clusterExport(cl, "gapfill_met_tier1")}
if(exists("gapfill_met_tier2"))   {clusterExport(cl, "gapfill_met_tier2")}
if(exists("gapfill_flux"))   {clusterExport(cl, "gapfill_flux")}
if(exists("min_yrs"))   {clusterExport(cl, "min_yrs")}

#Loops through sites
clusterMap(cl=cl, function(w,x) tryCatch(convert_fluxnet_to_netcdf(infile=w, site_code=x, out_path=out_path_lt,
                                                                   datasetname="LaThuile",
                                                                   met_gapfill="statistical", 
                                                                   flux_gapfill="statistical",
                                                                   missing_met=missing_met, missing_flux=missing_flux,
                                                                   gapfill_met_tier1=gapfill_met_tier1,
                                                                   gapfill_met_tier2=gapfill_met_tier2,
                                                                   gapfill_flux=gapfill_flux, min_yrs=min_yrs, 
                                                                   model="CABLE", include_all_eval=TRUE,
                                                                   check_range_action="warn",
                                                                   copyfill=365, regfill=365),
                                                                error = function(e) NULL),
                                                                w=infiles, x=site_codes)
stopCluster(cl)




# 
# mapply( function(w,x) convert_fluxnet_to_netcdf(infile=w, site_code=x, out_path=out_path_lt,
#                                                                    datasetname="LaThuile",
#                                                                    met_gapfill="statistical", 
#                                                                    flux_gapfill="statistical",
#                                                                    missing=missing, gapfill_all=gapfill_all,
#                                                                    min_yrs=min_yrs, model="CABLE",
#                                                                    include_all_eval=TRUE,
#                                                                    check_range_action="truncate",
#                                                                    copyfill=365, regfill=365),
#                                      
#            w=infiles, x=site_codes)


##############
### OzFlux ###
##############

#Input path
in_path <- paste0(path, '/OzFlux/Original_data/')

#Outputs will be saved to this directory
out_path_pre <- paste(path, "/OzFlux/Pre-processed_OzFlux_data/", sep="/")
out_path_oz  <- paste(out_path, "/OzFlux_sites/", sep="/")

#Remove path
unlink(out_path_pre, recursive = TRUE)
unlink(out_path_oz, recursive = TRUE)


#Find sitefiles
site_files <- list.files(in_path, full.names=TRUE)


### Pre-process files ###
lapply(site_files, preprocess_OzFlux, outpath=out_path_pre)


#Fluxnet site codes (having to set these manually for now, should add to pre-processing)
site_codes <- list(AdelaideRiver         = "AU-Ade",
                   AliceSpringsMulga     = "AU-ASM",           #NEW
                   Calperum              = "AU-Cpr", 
                   CapeTribulation       = "AU-Ctr",
                   CowBay                = "AU-Cow",         
                   CumberlandPlain       = "AU-Cum",
                   DalyPasture           = "AU-DaP", #Site code only DaP on ozflux website
                   DalyUncleared         = "AU-DaS",
                   DryRiver              = "AU-Dry",
                   Emerald               = "AU-Emr", #Not provided on ozflux website, maybe Arcturus? Taking site code from site metadata file  
                   FoggDam               = "AU-Fog",           #NEW
                   Gingin                = "AU-Gin",  
                   GreatWesternWoodlands = "AU-GWW",
                   HowardSprings         = "AU-How",   
                   Litchfield            = "AU-Lit",            #NEW
                   #Loxton                = "AU-Lox",   #Less than one year of data
                   Otway                 = "AU-Otw",           
                   RedDirtMelonFarm      = "AU-RDF", #Not provided on ozflux website, taking site code from site metadata file 
                   Ridgefield            = "AU-Rgf",            #NEW
                   RiggsCreek            = "AU-Rig",
                   RobsonCreek           = "AU-Rob",            #NEW
                   Samford               = "AU-Sam",
                   SturtPlains           = "AU-Stp",
                   TiTreeEast            = "AU-TTE",             #NEW
                   Tumbarumba            = "AU-Tum",
                   WallabyCreek          = "AU-Wac",            #NEW
                   Warra                 = "AU-Wrr",             #NEW
                   Whroo                 = "AU-Whr",
                   WombatStateForest     = "AU-Wom",
                   Yanco                 = "AU-Ync"
)


#Get sites (listed manually above for now, otherwise can't get site code)
sites <- names(site_codes)

#Find input files
in_files_oz <- unlist(sapply(sites, function(x) list.files(path=out_path_pre, pattern=x, full.names=TRUE)))


#Check that have as many original and pre-processed sites
if (length(in_files_oz) != length(site_files)) stop("Check why sites don't match")



### Process data ###

#Initialise clusters (using 2 cores here)
cl <- makeCluster(getOption('cl.cores', ncl))

#Import variables to cluster
clusterExport(cl, "out_path_oz")
clusterExport(cl, "convert_fluxnet_to_netcdf")
if(exists("missing_met"))   {clusterExport(cl, "missing_met")}
if(exists("missing_flux"))   {clusterExport(cl, "missing_flux")}
if(exists("gapfill_met_tier1"))   {clusterExport(cl, "gapfill_met_tier1")}
if(exists("gapfill_met_tier2"))   {clusterExport(cl, "gapfill_met_tier2")}
if(exists("gapfill_flux"))   {clusterExport(cl, "gapfill_flux")}
if(exists("min_yrs"))   {clusterExport(cl, "min_yrs")}

#Loops through sites
clusterMap(cl=cl, function(w,x) tryCatch(convert_fluxnet_to_netcdf(infile=w, site_code=x, out_path=out_path_oz,
                                                                   datasetname="OzFlux", met_gapfill="statistical", 
                                                                   flux_gapfill="statistical",
                                                                   missing_met=missing_met, missing_flux=missing_flux,
                                                                   gapfill_met_tier1=gapfill_met_tier1,
                                                                   gapfill_met_tier2=gapfill_met_tier2,
                                                                   gapfill_flux=gapfill_flux, min_yrs=min_yrs,
                                                                   include_all_eval=TRUE, check_range_action="warn", 
                                                                   model="CABLE"),
                                                                   error = function(e) NULL),
                                                                w=in_files_oz, x=site_codes)
stopCluster(cl)






convert_fluxnet_to_netcdf(infile=in_files_oz[1], site_code=site_codes[1], out_path=out_path_oz,
                          datasetname="OzFlux", met_gapfill="statistical", 
                          flux_gapfill="statistical",
                          missing_met=missing_met, missing_flux=missing_flux,
                          gapfill_met_tier1=gapfill_met_tier1,
                          gapfill_met_tier2=gapfill_met_tier2,
                          gapfill_flux=gapfill_flux, min_yrs=min_yrs,
                          include_all_eval=TRUE, check_range_action="warn", 
                          model="CABLE")
                          


###################################
### Combine non-duplicate sites ###
###################################

### Find sites and files ###

#List all datasets
datasets <- list.files(out_path, full.names=TRUE)

#Find all site files
all_files <- unlist(sapply(datasets, function(x) list.files(paste0(x, "/Nc_files/Met/"))))

#Find dataset for each file
all_files_dataset <- sapply(all_files, function(x) strsplit(x, "_")[[1]][3])

#Get site codes
all_sites <- sapply(all_files, substr, start=1, stop=6)

#Find unique sites
unique_sites <- unique(all_sites)


### Create output directory for non-duplicate sites ###

out_path_all <- paste0(out_path, "/all_sites_no_duplicates/")

#Remove old directory if it exists
unlink(out_path_all, recursive=TRUE)

#Create output paths for met and flux
dir.create(paste0(out_path_all, "/Nc_files/Met/"), recursive=TRUE)
dir.create(paste0(out_path_all, "/Nc_files/Flux/"), recursive=TRUE)
dir.create(paste0(out_path_all, "/Nc_files/Figures/"), recursive=TRUE)



### Loop through unique sites ###

for (s in 1:length(unique_sites)) {
  
  #Find number of site replicates
  ind <- which(all_sites == unique_sites[s])
  
  
  #If only one instance, copy that file
  if (length(ind) == 1) {
    
    #Met
    file.copy(from=list.files(paste0(out_path, "/", all_files_dataset[ind],"_sites/Nc_files/Met"),
              pattern=all_sites[ind], full.names=TRUE), to=paste0(out_path_all, "/Nc_files/Met"))
            
    #Flux
    file.copy(from=list.files(paste0(out_path, "/", all_files_dataset[ind],"_sites/Nc_files/Flux"),
              pattern=all_sites[ind], full.names=TRUE), to=paste0(out_path_all, "/Nc_files/Flux"))
    
    
  #If available in multiple datasets
  } else {
    
    #Prioritise OzFlux
    if(any(all_files_dataset[ind] == "OzFlux")) {
      dataset_to_use <- "OzFlux"
      
    #Then use FLUXNET2015
    } else if (any(all_files_dataset[ind] == "FLUXNET2015")) {
      dataset_to_use <- "FLUXNET2015"
    
    #And finally La Thuile 
    } else if (any(all_files_dataset[ind] == "LaThuile")) {
      dataset_to_use <- "LaThuile"
      
    #Dataset not recognised (stop)
    } else {
      stop ("Dataset not recognised")
    }
    
    
    
    #Copy files
    
    #Met
    file.copy(from=list.files(paste0(out_path, "/", dataset_to_use,"_sites/Nc_files/Met"),
                              pattern=all_sites[ind[1]], full.names=TRUE), to=paste0(out_path_all, "/Nc_files/Met"))
    
    #Flux
    file.copy(from=list.files(paste0(out_path, "/", dataset_to_use,"_sites/Nc_files/Flux"),
                              pattern=all_sites[ind[1]], full.names=TRUE), to=paste0(out_path_all, "/Nc_files/Flux"))
    
    #Figures
    file.copy(from=list.files(paste0(out_path, "/", dataset_to_use,"_sites/Nc_files/Figures"),
                              pattern=all_sites[ind[1]], full.names=TRUE), to=paste0(out_path_all, "/Nc_files/Figures"),
              recursive = TRUE)
    
    
  }
}







