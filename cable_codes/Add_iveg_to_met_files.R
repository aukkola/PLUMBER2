library(ncdf4)
library(readr)


#clear R environment
rm(list=ls(all=TRUE))

#Set path
path <- "/srv/ccrc/data04/z3509830/Fluxnet_data/All_flux_sites_processed_PLUMBER2/"


#Open Nc files
met_files_old <- list.files(paste0(path, "/Post-processed_PLUMBER2_outputs/Nc_files/Met"),
                        pattern=".nc", full.names = TRUE)


#Copy to a new folder
outdir <- paste0(paste0(path, "/Post-processed_PLUMBER2_outputs/Nc_files/Met_CABLE/"))
dir.create(outdir)

file.copy(met_files_old, outdir)

#Get new files
met_files <- list.files(outdir, pattern=".nc", full.names = TRUE)


met_nc <- lapply(met_files, function(x) nc_open(x, write=TRUE))


#Get site codes
site_codes <- sapply(met_nc, function(x) ncatt_get(x, varid=0, "site_code")$value)


#Read site metadata
metadata_url <- "https://raw.githubusercontent.com/aukkola/FluxnetLSM/master/data/Site_metadata.csv"
site_data    <- read_csv(url(metadata_url), col_names=c("SiteCode", "Exclude", "Exclude_reason",
                                                        "CABLE_PFT", "CABLE_patchfrac", "Source_CABLE_PFT" ), 
                         col_types=list("c", "l", "c", "c", "c", "c"))





#Loop through sites

#CABLE doesn't like files with 1 tile, not writing patchfrac to those files

for (s in 1:length(site_codes)) {
  
  
  #Convert iveg to numeric
  iveg_char <- site_data$CABLE_PFT[which(site_data$SiteCode == site_codes[s])]
  
  iveg <- as.numeric(strsplit(iveg_char, ",")[[1]])
  
  
  #Convert patchfrac to numeric
  patchfrac_char <- site_data$CABLE_patchfrac[which(site_data$SiteCode == site_codes[s])]
  
  patchfrac <- as.numeric(strsplit(patchfrac_char, ",")[[1]])
  
  
  
  
  #Define patch dimension
  dimptc <- ncdim_def("patch", "", patchfrac, unlim=FALSE)
  
  #Get x and y dimensions from file
  dimx <- met_nc[[s]]$dim$x
  dimy <- met_nc[[s]]$dim$y

  
  
  #Define iveg and patchfrac variables
  
  if (length(patchfrac) > 1) {
    iveg_var <- ncvar_def(name="iveg", units="-", dim=list(dimx, dimy,dimptc),
                          missval=-9999, prec="double")  
    
  } else {
    iveg_var <- ncvar_def(name="iveg", units="-", dim=list(dimx, dimy),
                          missval=-9999, prec="double")  

  }
  
  

  #Define patchfrac variable
  patchfrac_var <- ncvar_def(name="patchfrac", units="-", dim=list(dimx, dimy, dimptc),
                             missval=-9999, prec="double")
  
  
  #Add variables to file
  met_nc[[s]] <- ncvar_add(met_nc[[s]], iveg_var, verbose=FALSE)
  
  if (length(patchfrac) > 1) {
    met_nc[[s]] <- ncvar_add(met_nc[[s]], patchfrac_var, verbose=FALSE)
  }
  
  
  #Add values
  
  
  if (length(patchfrac) > 1) {
    ncvar_put(met_nc[[s]], varid=iveg_var, vals=iveg, start=c(1,1,1), count=c(1,1, length(patchfrac)))
    
    ncvar_put(met_nc[[s]], varid=patchfrac_var, vals=patchfrac, start=c(1,1,1), count=c(1,1, length(patchfrac)))
    
  } else {
    
    ncvar_put(met_nc[[s]], varid=iveg_var, vals=iveg, start=c(1,1), count=c(1,1))
    
  }
  
  #Close file handle
  nc_close(met_nc[[s]])
  
}












