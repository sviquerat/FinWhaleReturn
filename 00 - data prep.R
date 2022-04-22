rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #only works in RStudio!
SCRIPTDIR<-file.path(getwd(),'SCRIPT')
for (f in list.files(SCRIPTDIR,pattern = '*.R',full.names = T)){ source(f)}

#### SPATIAL DATA #####
#### load, prepare and save data used during the analysis
#### includes manually created (clipped) survey boundaries 

library(raster)

shelf<-rgdal::readOGR(file.path(SPDIR,'shelf_edge_Clip.gpkg')) #Shelf edge from Herr et al. 2016
shelf_poly<-rgdal::readOGR(file.path(SPDIR,'shelf_edge_Clip_Polygon.gpkg')) #same as previous data, as polygon
predBoundary<-rgdal::readOGR(file.path(SPDIR,'Prediction_Boundary.gpkg')) #Boundary for prediction grid
coastline<-rgdal::readOGR(file.path(SPDIR,'Antarctica_Clip.gpkg')) #'high' resolution coastline of Antarctica
coastline_simple <- rgeos::gSimplify(coastline, tol = 1000, topologyPreserve = T) #simplify coastline

# clip prediction area with coastline data
prediction <- rgeos::gDifference(predBoundary, coastline_simple) #throws warning and error in latest rgeos (0.5-5) but works fine
prediction<-sp::SpatialPolygonsDataFrame(prediction,data.frame(id=1))
area_km2<-rgeos::gArea(prediction) #calculate area of clipped prediction area in km²
prediction$area_km2<-area_km2

# create subregions of islands of interest
extent<-as(raster::extent(-2700000,-2530000,1780000,1920000),'SpatialPolygons')
extent@polygons[[1]]@ID<-'Elephant Island'
extent2<-as(raster::extent(-2718000, -2600000, 1557000,1678000),'SpatialPolygons')
extent2@polygons[[1]]@ID<-'King George Island'
extent3<-as(raster::extent(-2699000, -2620000, 1400000,1480000),'SpatialPolygons')
extent3@polygons[[1]]@ID<-'Livingstone Island'
raster::crs(extent)<-raster::crs(extent2)<-raster::crs(extent3)<-ANT_POL_STEREO
p<-list("Elephant" = extent,"KingGeorge" = extent2,"Livingstone" = extent3)
islands = SpatialPolygons(lapply(p, function(x){x@polygons[[1]]}))
raster::crs(islands)<-ANT_POL_STEREO

rgdal::writeOGR(prediction,dsn=file.path(SPDIR,'clipped_survey_area.gpkg'),layer='Clipped survey area',driver = 'GPKG', overwrite_layer = T)
save(prediction,coastline_simple, shelf, file=file.path(DATRESDIR,'PS112_spatialData.RData'),compress='gzip')

#### SURVEY DATA ####
#### load distance sampling data set and augment with covariates

load(file.path(DATADIR,'PS112_HELI_DATA.RData'))
seg<-segmentate(data,'transect','track_length_km',LIMIT=5)$data
data$seg_label<-seg$seg_label
data$seg_length_km<-seg$seg_length

save(data,data.fin,file=file.path(DATRESDIR,'PS112_modified_survey_data.RData'),compress='gzip')

# visualise all components of data prep:
plot(prediction,col=adjustcolor('lightblue',.2)) #Polygon of preidction area
#lines(islands,col=adjustcolor('red',.5),lty=2) #Extent of island subregions
lines(coastline_simple) #Simplified antarctic coastline
lines(shelf,col='orange',lwd=2) #Polygon of Antarctic Shelf edge
points(data$y~data$x,col='red',pch=16,cex=.1) #line transect survey data
points(data.fin$y~data.fin$x,col='green',pch=16,cex=.5) #valid sightings
axis(1)
axis(2)
box()
legend('topright',legend=c(
  'Clipped survey area','Shelf edge', 'Line transect data','fin whale records on effort'),
  col=c(adjustcolor('lightblue',1),'orange','red','green'),
  pch=c(17,NA,16,16),lty=c(NA,2,NA,NA),cex=.7)
title(main='Summary of data created in this script')
