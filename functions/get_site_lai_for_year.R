#Function to get Copernicus LAI for each year

get_site_lai_for_year <- function(coords, lai_averaged) {
  
  #Get time series !NB NEED TO FIX HERE TO REMOVE SCALING  1/9 (caused by error in focal function) !!!!!!!!!!!!!!!!!!!!!!!
  year_lai_vals <- as.vector(extract(lai_averaged, coords) /(1/9) )
  
  #If didn't find any, find nearest non-NA pixel
  if (any (is.na(year_lai_vals))) {
    
    #Crop to smaller size to speed up distance calculation
    temp_lai <- crop(lai_averaged, extent(c(coords[1]-3, coords[1]+3,
                                            coords[2]-3, coords[2]+3)))
    
    
    #lapply to layers
    year_lai_vals <- sapply(1:nlayers(temp_lai), function(l) sample_raster_NA(temp_lai[[l]], coords))
    
    
  }
  
  return(year_lai_vals)
}