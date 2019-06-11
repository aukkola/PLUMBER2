library(ncdf4)


#clear R environment
rm(list=ls(all=TRUE))


#Set path
path <- "/srv/ccrc/data04/z3509830/Fluxnet_data/All_flux_sites_processed_no_filtering/"

#Source function
source(paste0(path, "/scripts/functions/plot_timeseries.R"))


#Variables to plot
vars <- c("SWdown", "LWdown", "Precip", "Tair", "Qair", "Wind", "CO2air", 
          "Rnet", "Qle", "Qh", "Qg", "Ebal")
  
#Variable types
var_type <- c("Met", "Met", "Met", "Met", "Met", "Met", "Flux",
              "Flux", "Flux", "Flux", "Flux", "Flux") 



#Find sites
site_files <- list.files(paste0(path, "/all_sites_no_duplicates/Nc_files/Met"),
                          full.names=TRUE)

#Open file handles
site_nc_met <- lapply(site_files, nc_open)

#Get site codes
site_codes <- sapply(site_nc_met, function(x) ncatt_get(x, varid=0, "site_code")$value)


#Set output directory
outdir <- paste0(path, "/PLUMBER2_site_selection_figs/")
dir.create(outdir)



#Loop through sites
for (s in 1:length(site_codes)) {
  
  
  #Get matching flux file
  nc_flux <- nc_open(list.files(paste0(path, "/all_sites_no_duplicates/Nc_files/Flux"), 
                                pattern=site_codes[s], full.names=TRUE))
  
  
  #Get time (used to var lengths etc)
  time <- ncvar_get(nc_flux, "time")
  
  #Get timing information
  timing <- GetTimingNcfile(nc_flux)  
  
  
  #Set up file
  outfile <- paste0(outdir, "/", site_codes[s], "_timeseries_plot.png")
  png(outfile, height=1300, width=2800, res=50, pointsize=20)
  
  par(mfcol=c(3,4))
  
  
  #Loop through variables
  for (v in 1:length(vars)) {
    
    #Get data
    
    #Energy balance
    if (vars[v] == "Ebal") {
      
      rnet <- tryCatch(ncvar_get(nc_flux, "Rnet"), error=function(e) rep(0, length(time)))
      qle  <- tryCatch(ncvar_get(nc_flux, "Qle"), error=function(e) rep(0, length(time)))
      qh   <- tryCatch(ncvar_get(nc_flux, "Qh"), error=function(e) rep(0, length(time)))
      
      #Look for ground heat flux, set to zero if not found
      qg <- tryCatch(ncvar_get(nc_flux, "Qg"), error=function(e) rep(0, length(time)))
      
      #Calculate ratio (Rnet-G) / (Qle + Qh)
      var_data <- (rnet - qg) / (qle + qh)
      
      qc_data <- NA
      
      data_unit <- "-"
      
      
      #Cap values to 0-5
      var_data[var_data > 1.5] <- 1.5
      var_data[var_data < 0]  <- 0
      var_data[which(is.infinite(var_data))] <- 1.5
      
      
      
    #All other variables
    } else {
      
      #Set file for met variable
      if (var_type[v] == "Met") {
        nc <- site_nc_met[[s]]
        
      #Set file for  flux variable  
      } else if (var_type[v] == "Flux") { 
        nc <- nc_flux
      }
      
        #Get variable data
        var_data <- tryCatch(ncvar_get(nc, vars[v]), error=function(e) rep(NA, length(time)))
        
        #Get QC data
        qc_data <- tryCatch(ncvar_get(nc, paste0(vars[v], "_qc")), error=function(e) NA)
        
        #Get data units
        data_unit <- tryCatch(ncatt_get(nc, varid=vars[v], "units")$value, error=function(e) NA)
    
    }
      
    #Convert precip and Tair units
    if (vars[v] == "Precip") {
      
      var_data <- var_data * (time[2] - time[1])
      data_unit <- paste("mm/", (time[2] - time[1])/60, "min", sep="")
    }
    
    if (vars[v] == "Tair") {
      
      var_data <- var_data - 273.15
      data_unit <- "deg C"
    }
    
    

    #Extract QC data and replace all gap-filled values with 0
    # and measured with 1 (opposite to Fluxnet but what PALS expects)
    if (any(!is.na(qc_data))) {
      
      qc_data[qc_data > 0]  <- 2 #replace gap-filled values with a temporary value
      qc_data[qc_data == 0] <- 1 #set measured to 1
      qc_data[qc_data == 2] <- 0 #set gap-filled to 0
      
      #If first value missing, set to measured (to avoid an error when PALS
      #checks if first value -1)
      if(is.na(qc_data[1])){ qc_data[1] <- 0}
      
      
      #Else set to PALS option corresponding to no QC data
    } else {
      qc_data <- matrix(-1, nrow = 1, ncol = 1)
    }
    
    
  
    #Plot file
    Timeseries(obslabel=site_codes[s], tsdata=as.matrix(var_data),
               varname=vars[v],
               ytext=paste(vars[v], " (", data_unit, ")", sep=""), 
               legendtext=vars[v],
               plotcex=2 , timing=timing, 
               smoothed = FALSE, winsize = 1, 
               plotcolours="black",
               vqcdata = as.matrix(qc_data),
               na.rm=TRUE)
  
  
  
  } #vars
  
  dev.off()
  
} #sites



  
  




  



