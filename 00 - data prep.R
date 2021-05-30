rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #only works in RStudio!
SCRIPTDIR<-file.path(getwd(),'SCRIPT')
for (f in list.files(SCRIPTDIR,full.names = T)){ source(f)}

#### SPATIAL DATA #####
library(raster)

shelf<-rgdal::readOGR(file.path(SPDIR,'shelf_edge_Clip.gpkg'))
shelf_poly<-rgdal::readOGR(file.path(SPDIR,'shelf_edge_Clip_Polygon.gpkg'))
coastline<-rgdal::readOGR(file.path(SPDIR,'Antarctica_Clip.gpkg'))
predBoundary<-rgdal::readOGR(file.path(SPDIR,'Prediction_Boundary.gpkg'))
outerBoundary<-rgdal::readOGR(file.path(SPDIR,'Survey_Area_Mask.gpkg'))

coastline_simple <- rgeos::gSimplify(coastline, tol = 1000, topologyPreserve = T)
prediction <- rgeos::gDifference(predBoundary, coastline_simple)
prediction<-sp::SpatialPolygonsDataFrame(prediction,data.frame(id=1))

#prediction_km<-sp::spTransform(prediction,ANT_POL_STEREO_km)
area_km2<-rgeos::gArea(prediction)
prediction$area_km2<-area_km2

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
save(prediction,islands,coastline_simple, shelf, file=file.path(DATRESDIR,'PS112_spatialData.RData'),compress='gzip')

#### SURVEY DATA ####
load(file.path(DATADIR,'PS112_HELI_DATA.RData'))
LL<-SpatialPoints(cbind(data$lon,data$lat),WGS84)
LL<-sp::spTransform(LL,crs(ANT_POL_STEREO))
data$x<-LL@coords[,1]
data$y<-LL@coords[,2]
LL<-SpatialPointsDataFrame(LL,data)
rgdal::writeOGR(LL,dsn=file.path(SPDIR,'PS112_Survey_data.gpkg'),layer='PS112 Survey Data',driver = 'GPKG', overwrite_layer = T)
seg<-segmentate(data,'transect','track_length_km',LIMIT=5)$data
data$seg_label<-seg$seg_label
data$seg_length_km<-seg$seg_length

save(data,file=file.path(DATRESDIR,'PS112_modified_survey_data.RData'),compress='gzip')

# visualise all components of data prep:
plot(prediction,col=adjustcolor('lightblue',.2)) #Polygon of preidction area
plot(islands,add=T,border=adjustcolor('red',.5)) #Extent of island subregions
plot(coastline_simple,add=T) #Simplified antarctic coastline
plot(shelf,add=T,col='green',lwd=2) #Polygon of Antarctic Shelf edge
plot(LL,add=T,col='red',pch=16,cex=.4) #
axis(1)
axis(2)
box()
