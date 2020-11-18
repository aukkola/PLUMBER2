#Need to run on Hurricane !! R version too old on other servers

library(MODISTools)
library(ncdf4)


#clear R environment
rm(list=ls(all=TRUE))

#Set path
path <- "/srv/ccrc/data04/z3509830/Fluxnet_data/All_flux_sites_processed_PLUMBER2/"





sites_to_fetch <- data.frame(site_name=c("AU-Cpr", "AU-DaS", "AU-Dry", "AU-Gin",
                                         "AU-GWW", "AU-How", "AU-Lit", "AU-Otw",
                                         "AU-TTE", "BW-Ma1", "CA-NS6", "CA-NS7",
                                         "CA-SF3", "SD-Dem", "US-Whs"))

sites_to_fetch$lat <- c(-34.00206, -14.159283, -15.2588, -31.375,
                        -30.1913, -12.4952, -13.17904, -38.53234, 
                        -22.287, -19.9165, 55.91669846, 56.63579941,
                        54.09159851, 13.2829, 31.74379921)

sites_to_fetch$lon <- c(140.58912, 131.388, 132.3706, 115.65,
                        120.6541, 131.15, 130.7945, 142.81681,
                        133.64, 23.56033, -98.96440125, -99.94830322,
                        -106.0053, 30.4783, -110.0522)

                        	


#MODIS/Terra+Aqua Leaf Area Index/FPAR 8-Day L4
product <- "MOD44B"


#Cells to obtain around site
#n.b. tried using 0.5 km but this makes the code crash (???)
km <- 1


#Get LAI
treefrac <- mt_batch_subset(df=sites_to_fetch, product=product, band="Percent_Tree_Cover", start = "2000-01-01",
                            end = format(Sys.time(), "%Y-%m-%d"), km_lr = km, km_ab = km,
                            ncores = 10, internal=TRUE)

sd <- mt_batch_subset(df=sites_to_fetch, product=product, band="Percent_Tree_Cover_SD", start = "2000-01-01",
                               end = format(Sys.time(), "%Y-%m-%d"), km_lr = km, km_ab = km,
                               ncores = 10, internal=TRUE)

qc <- mt_batch_subset(df=sites_to_fetch, product=product, band="Quality", start = "2000-01-01",
                               end = format(Sys.time(), "%Y-%m-%d"), km_lr = km, km_ab = km,
                               ncores = 10, internal=TRUE)

nonveg <- mt_batch_subset(df=sites_to_fetch, product=product, band="Percent_NonVegetated", start = "2000-01-01",
                          end = format(Sys.time(), "%Y-%m-%d"), km_lr = km, km_ab = km,
                          ncores = 10, internal=TRUE)

nonveg_sd <- mt_batch_subset(df=sites_to_fetch, product=product, band="Percent_NonVegetated_SD", start = "2000-01-01",
                             end = format(Sys.time(), "%Y-%m-%d"), km_lr = km, km_ab = km,
                             ncores = 10, internal=TRUE)



# Extracting pixels in the centre and immediately around it (*)
# These correspond to a radius of 500m around site coordinates

pixel_no <- c(7, 8, 9, 12, 13, 14, 17, 18, 19)

qc_flags <- 0


##########################
### Loop through sites ###
##########################


mean_treefrac <- vector(length=nrow(sites_to_fetch))
mean_nonveg   <- vector(length=nrow(sites_to_fetch))



for (s in 1:nrow(sites_to_fetch)) {
  
  tf_site <- treefrac[which(treefrac$site == sites_to_fetch$site_name[s]),]
  sd_site <- sd[which(treefrac$site == sites_to_fetch$site_name[s]),]
  qc_site <- qc[which(treefrac$site == sites_to_fetch$site_name[s]),]
  
  nonveg_site    <- nonveg[which(treefrac$site == sites_to_fetch$site_name[s]),]
  nonveg_sd_site <- nonveg_sd[which(treefrac$site == sites_to_fetch$site_name[s]),]
  
  
  
  
  #Number of time steps
  
  #If data files have been acquired on different days from MODIS
  #server, they might be different lengths so take minimum. adjust further down
  
  no_tsteps <- min(nrow(tf_site), nrow(sd_site), nrow(qc_site)) / max(tf_site$pixel)
  
  
  
  #Extract 3 x3 pixels
  treefrac_pixel  <- matrix(nrow=no_tsteps, ncol=length(pixel_no))
  sd_pixel        <- matrix(nrow=no_tsteps, ncol=length(pixel_no))
  qc_pixel        <- matrix(nrow=no_tsteps, ncol=length(pixel_no))
  nonveg_pixel    <- matrix(nrow=no_tsteps, ncol=length(pixel_no))
  nonveg_sd_pixel <- matrix(nrow=no_tsteps, ncol=length(pixel_no))
  
  #Save time stamps
  lai_time <- as.Date(tf_site$calendar_date[which(tf_site$pixel == pixel_no[1])])
  
  #Loop through pixels
  for (p in 1:length(pixel_no)) {
    
    #Get time series for pixel and scale using scale factor (and adjust for different lengths with min_dim)
    treefrac_pixel[,p] <- tf_site$value[which(tf_site$pixel == pixel_no[p])][1:no_tsteps] * 
                          as.numeric(tf_site$scale[1])
    sd_pixel[,p]       <- sd_site$value[which(sd_site$pixel == pixel_no[p])][1:no_tsteps] * 
                          as.numeric(sd_site$scale[1])
    qc_pixel[,p]       <- qc_site$value[which(sd_site$pixel == pixel_no[p])][1:no_tsteps] 
    
    nonveg_pixel[,p]   <- nonveg_site$value[which(nonveg_site$pixel == pixel_no[p])][1:no_tsteps] * 
                          as.numeric(nonveg_site$scale[1])
    nonveg_sd_pixel[,p] <- nonveg_sd_site$value[which(nonveg_sd_site$pixel == pixel_no[p])][1:no_tsteps] * 
                           as.numeric(nonveg_sd_site$scale[1])
    
    
  }
  
  
  
  
  ##################################
  ### Mask out poor quality data ###
  ##################################
  
  #Mask out where QC flag 
  treefrac_pixel  <- replace(treefrac_pixel, !(qc_pixel %in% qc_flags), NA)
  sd_pixel        <- replace(sd_pixel, !(qc_pixel %in% qc_flags), NA)
  
  nonveg_pixel    <- replace(nonveg_pixel, !(qc_pixel %in% qc_flags), NA)
  nonveg_sd_pixel <- replace(nonveg_sd_pixel, !(qc_pixel %in% qc_flags), NA)
  
  
  #Also mask out where sd really low (likely cloud effects)
  sd_pixel        <- replace(sd_pixel, sd_pixel < 0.1, NA)
  treefrac_pixel  <- replace(treefrac_pixel, is.na(sd_pixel), NA)
  
  nonveg_sd_pixel <- replace(nonveg_sd_pixel, nonveg_sd_pixel < 0.1, NA)
  nonveg_pixel  <- replace(nonveg_pixel, is.na(nonveg_sd_pixel), NA)
  
  
  #Set fill values to missing
  treefrac_pixel <- replace(treefrac_pixel, treefrac_pixel > 100, NA)
  nonveg_pixel   <- replace(nonveg_pixel, nonveg_pixel > 100, NA)
  
  
  ### Average each time step ###
  
  #Initialise lai time series
  treefrac_ts <- vector(length=no_tsteps)
  nonveg_ts   <- vector(length=no_tsteps)
  
  
  #Loop through time steps
  for (t in 1:no_tsteps) {
    
    ### Tree fraction ###
    
    #If no values available
    if (all(is.na(treefrac_pixel[t,]))) {
      
      treefrac_ts[t] <- NA

      #If values available
    } else {
      
      #Weight grid cell estimates by their standard deviation
      #Following Martin's method (https://github.com/mdekauwe/get_MODIS_LAI_australia/blob/master/build_modis_climatology.py),
      #but normalising by sum of standard deviations
      
      sd_vals <- sd_pixel[t,]
      
      weights = (1/sd_vals**2) / sum(1/sd_vals**2, na.rm=TRUE)
      
      #Check that weights sum up to 1 (because of a precision issue presumably,
      #rounding to 5 decimals, otherwise might not equal 1 even when correct)
      if (round(sum(weights, na.rm=TRUE), 5) != 1) stop("Weighting not correct")
      
      #Calculate weighted average
      treefrac_ts[t] <- weighted.mean(treefrac_pixel[t,], w=weights, na.rm=TRUE)
      
    }
  }
  
  mean_treefrac[s] <- mean(treefrac_ts, na.rm=TRUE)
  
  
  #Loop through time steps
  for (t in 1:no_tsteps) {
    
    ### Non-vegetated ###
    
    #If no values available
    if (all(is.na(nonveg_pixel[t,]))) {
      
      nonveg_ts[t] <- NA
      
      #If values available
    } else {
      
      #Weight grid cell estimates by their standard deviation
      #Following Martin's method (https://github.com/mdekauwe/get_MODIS_LAI_australia/blob/master/build_modis_climatology.py),
      #but normalising by sum of standard deviations
      
      sd_vals <- nonveg_sd_pixel[t,]
      
      weights = (1/sd_vals**2) / sum(1/sd_vals**2, na.rm=TRUE)
      
      #Check that weights sum up to 1 (because of a precision issue presumably,
      #rounding to 5 decimals, otherwise might not equal 1 even when correct)
      if (round(sum(weights, na.rm=TRUE), 5) != 1) stop("Weighting not correct")
      
      #Calculate weighted average
      nonveg_ts[t] <- weighted.mean(nonveg_pixel[t,], w=weights, na.rm=TRUE)
      
    }
  }
  
  mean_nonveg[s] <- mean(nonveg_ts, na.rm=TRUE)
}

  







