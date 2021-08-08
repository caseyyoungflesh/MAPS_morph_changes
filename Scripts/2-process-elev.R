#######################
# Mosaic and process elevation data
#
# # GMTED2010 Global Elevation Product - 30 arc seconds (~ 1km)
#https://www.usgs.gov/centers/eros/science/usgs-eros-archive-digital-elevation-global-multi-resolution-terrain-elevation?qt-science_center_objects=0#qt-science_center_objects
# Region specified on Earth Explorer and downloaded using the EE bulk download tool
# Only mean elevation for each pixel retained to save space (each metric [min, max, std, ...] was downloaded as a different .tiff)
#######################


# input -------------------------------------------------------------------

dir <- 'XXXX'
MAPS_process_date <- '2021-07-28'


# load packages -----------------------------------------------------------

library(gdalUtils)
library(sp)


# mosaic GMTED tifs -------------------------------------------------------

setwd(paste0(dir, 'Data/GMTED2010_mean/'))

#mosaic GMTED2010 DEM tiles together with gdal
# system('gdalbuildvrt elev_mosaic.vrt *.tif')
# system('gdal_translate elev_mosaic.vrt elev_mosaic.tif')


# extract elev at stations ------------------------------------------------

#read in mosaic tif
elev_mosaic <- raster::raster('elev_mosaic.tif')

#extract elevation at each data point
setwd(paste0(dir, '/Data'))
maps_data <- readRDS(paste0('MAPS-processed-', MAPS_process_date, '.rds'))

#convert lat/long to spdf
maps_data_spdf <- sp::SpatialPointsDataFrame(cbind(maps_data$lng, maps_data$lat), 
                                             proj4string = elev_mosaic@crs, maps_data)

#extract elevation from elevation raster
maps_data$GMTED_elev <- raster::extract(elev_mosaic, maps_data_spdf, df = FALSE)


# save RDS ----------------------------------------------------------------

setwd(paste0(dir, '/Data'))
saveRDS(maps_data, paste0('MAPS-master-', MAPS_process_date, '.rds'))

