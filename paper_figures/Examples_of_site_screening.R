library(ncdf4)


#clear R environment
rm(list=ls(all=TRUE))


#Set path
path <- "/srv/ccrc/data04/z3509830/Fluxnet_data/All_flux_sites_processed_PLUMBER2/"

#Source function
source(paste0(path, "/scripts/functions/plot_timeseries_no_metrics.R"))


#Variables to plot
vars <- c("SWdown", "Precip", "Tair", "Qair", "Wind")

#Sites to plot
sel_sites <- c("AU-Lit", "BE-Bra", "US-Tw2")

#Find sites
site_files <- list.files(paste0(path, "/all_sites_no_duplicates/Nc_files/Met_with_LAI/"),
                         full.names=TRUE)

#Open file handles
site_nc_met <- lapply(site_files, nc_open)

#Get site codes
site_codes <- sapply(site_nc_met, function(x) ncatt_get(x, varid=0, "site_code")$value)


#Set output directory
outdir <- paste0(path, "/PLUMBER2_paper_figs/")
dir.create(outdir)


#Set up file
outfile <- paste0(outdir, "/Site_screening_timeseries_plots.png")
#png(outfile, height=11, width=8, res=1200, unit="in", pointsize=5)

png(outfile, height=3200, width=2520, res=50, pointsize=40)


par(mfcol=c(5,3))
par(mai=c(0.6, 1.8, 0.8, 0.2))
par(omi=c(1.1, 1.5, 1.8, 0.1))


#indices for selected sites
sel_ind <- sapply(sel_sites, function(x) which(site_codes == x))

#Loop through selected sites
for (s in sel_ind) {
  
  
  #Get time (used to var lengths etc)
  time <- ncvar_get(site_nc_met[[s]], "time")
  
  #Get timing information
  timing <- GetTimingNcfile(site_nc_met[[s]])  
  
  #QC information
  QC_measured <- 0
  if(substr(site_codes[s], 1, 3) == "AU-") QC_measured <- append(QC_measured, 10)
  
  
  
  #Loop through variables
  for (v in 1:length(vars)) {
    
    
    #Get data
    nc <- site_nc_met[[s]]
        
       
    #Get variable data
    var_data <- tryCatch(ncvar_get(nc, vars[v]), error=function(e) rep(NA, length(time)))
    
    #Get QC data
    qc_data <- tryCatch(ncvar_get(nc, paste0(vars[v], "_qc")), error=function(e) NA)
    
    #Get data units
    data_unit <- tryCatch(ncatt_get(nc, varid=vars[v], "units")$value, error=function(e) NA)
    
    
    #Variable y-label
    ylab <- paste(vars[v], " (", data_unit, ")", sep="")
    
    #Convert precip from kg/m2/s to mm/timestep
    if (vars[v] == "Precip") {
      
      var_data <- var_data * (time[2] - time[1])
      ylab <- paste(vars[v], " (mm/", (time[2] - time[1])/60, "min)", sep="")
    }
    
    #Tair (Kelvin to C)
    if (vars[v] == "Tair") {
      
      var_data <- var_data - 273.15
      ylab <- expression(bold("T"["air"]~"("*degree*C*")"))
    }
    
    #Swdown needs subscript
    if (vars[v] == "SWdown") {
      ylab <- expression(bold("SW"["down"]~"("*"W/m"^2*")"))
    }
    
    #Qair needs subscript
    if (vars[v] == "Qair") {
      ylab <- expression(bold("Q"["air"]~"(kg/kg)"))
    }
    
    
    
    #Extract QC data and replace all gap-filled values with 0
    # and measured with 1 (opposite to Fluxnet but what PALS expects)
    if (any(!is.na(qc_data))) {
      
      
      qc_data[!(qc_data %in% QC_measured)]  <- 2 #replace gap-filled values with a temporary value
      qc_data[qc_data %in% QC_measured]     <- 1 #set measured to 1
      qc_data[qc_data == 2] <- 0 #set gap-filled to 0
      
      #If first value missing, set to measured (to avoid an error when PALS
      #checks if first value -1)
      if(is.na(qc_data[1])){ qc_data[1] <- 0}
      
      
      #Else set to PALS option corresponding to no QC data
    } else {
      qc_data <- matrix(-1, nrow = 1, ncol = 1)
    }
    
    
    x_axis <- FALSE
    if(v == length(vars)) x_axis <- TRUE
    
    #Plot file
    Timeseries(obslabel="", tsdata=as.matrix(var_data),
               varname="",
               ytext="", 
               legendtext=vars[v],
               plotcex=2 , timing=timing, 
               smoothed = FALSE, winsize = 1, 
               plotcolours="black",
               vqcdata = as.matrix(qc_data),
               na.rm=TRUE, x_axis=x_axis)
    
    
    #Add labels
    
    #Site code
    if (v == 1) mtext(side=3, site_codes[s], line=2, cex=1.5, font=2)
    
    #Y-label
    if (s==sel_ind[1]) mtext(side=2, ylab,
                    line=4, cex=1.5, font=2)
    
  } #vars
} #sites


dev.off()


#Close file handles
lapply(site_nc_met, nc_close)








