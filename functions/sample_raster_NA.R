
#function for findind nearest non-NA grid cell
sample_raster_NA <- function(r, xy)
{
  apply(X = xy, MARGIN = 1, 
        FUN = function(xy) r[which.min(replace(distanceFromPoints(r, xy), is.na(r), NA))])
}

