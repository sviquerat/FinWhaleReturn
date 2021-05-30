rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #only works in RStudio!
SCRIPTDIR<-file.path(getwd(),'SCRIPT')
for (f in list.files(SCRIPTDIR,full.names = T)){ source(f)}

#### SPATIAL DATA #####
library(raster)
library(rgeos)

# SPRESDIR<-file.path(SPRESDIR,'ENV')
# dir.create(SPRESDIR,recursive=T,showWarnings = F)
# GFXDIR<-file.path(ENVDIR,'GFX')
# dir.create(GFXDIR,recursive=T,showWarnings = F)

cellsize<-5000

boundary<-rgdal::readOGR(file.path(SPDIR,'Prediction_Boundary.gpkg'))
shelf<-rgdal::readOGR(file.path(SPDIR,'shelf_edge_Clip.gpkg'))
shelf_poly<-rgdal::readOGR(file.path(SPDIR,'shelf_edge_Clip_Polygon.gpkg'))
surveyarea<-rgdal::readOGR(file.path(SPDIR,'Survey_Area_Mask.gpkg'))
coastline<-rgdal::readOGR(file.path(SPDIR,'Antarctica_Clip.gpkg'))

depth_i<-raster::raster(file.path(SPDIR,'ibcso_v1_bed.tif')) #60?S
depth_e<-raster::raster(file.path(SPDIR,'ETOPO1_DEM.tif')) #60?S
b<-spTransform(boundary,crs(depth_i))
depth_i<-crop(depth_i,b)

depth <- SpatialPixelsDataFrame(depth_i,data.frame(depth=depth_i@data@values))
raster::crs(depth)<-raster::crs(depth_i)

pos<-depth@coords[is.na(depth$depth),]
LL<-sp::SpatialPoints(cbind(pos[,1],pos[,2]),crs(depth))
LL<-sp::spTransform(LL,crs(depth_e))
depth$depth[is.na(depth$depth)]<-extract(depth_e,LL)

slope<-terrain(raster(depth),opt='slope')
aspect<-terrain(raster(depth),opt='aspect')
TPI<-terrain(raster(depth),opt='TPI')
TRI<-terrain(raster(depth),opt='TRI')
roughness<-terrain(raster(depth),opt='roughness')

slope <- SpatialPixelsDataFrame(slope,data.frame(slope=slope@data@values))
aspect <- SpatialPixelsDataFrame(aspect,data.frame(aspect=aspect@data@values))
TPI <- SpatialPixelsDataFrame(TPI,data.frame(TPI=TPI@data@values))
TRI <- SpatialPixelsDataFrame(TRI,data.frame(TRI=TRI@data@values))
roughness <- SpatialPixelsDataFrame(roughness,data.frame(roughness=roughness@data@values))

raster::crs(aspect)<-raster::crs(TPI)<-raster::crs(slope)<-raster::crs(TRI)<-raster::crs(roughness)<-raster::crs(depth)

pts<-sp::spsample(depth,type='regular',cellsize=cellsize)
raster::crs(pts)<-raster::crs(depth)

##### calculation of distance to shelf edge & distance to coastline
ds<-c()
pb<-txtProgressBar(min=0,max=nrow(pts@coords),style=3)
for (r in 1:nrow(pts@coords)){
  setTxtProgressBar(pb,r)
  cc<-pts[r,]
  dd = rgeos::gDistance(shelf,spTransform(cc,crs(shelf)))/1000 #distance to shelf is in km!
  ds<-c(ds,dd)
}
close(pb)

pts$dist2shelf <- ds
t<-sp::over(sp::spTransform(pts,raster::crs(shelf_poly)),shelf_poly)
on_shelf<-which(!is.na(t$OBJECTID))

# add distinction between on the shelf and off the shelf
pts$dist2shelf[on_shelf]<-pts$dist2shelf[on_shelf]*-1 #points on shelf have negative distance!
pts$on_shelf<-'oceanic'
pts$on_shelf[pts$dist2shelf<=0]<-'shelf'

ds<-c()
pb<-txtProgressBar(min=0,max=nrow(pts@coords),style=3)
for (r in 1:nrow(pts@coords)){
  setTxtProgressBar(pb,r)
  cc<-pts[r,]
  dd = rgeos::gDistance(coastline,spTransform(cc,crs(coastline)))/1000 #distance to coastline is in km!
  ds<-c(ds,dd)
}
close(pb)

pts$dist2coast <- ds
t<-over(spTransform(pts,crs(coastline)),coastline)

# add distinction between points on land and off land
on_land<-which(!is.na(t$OBJECTID))
pts$dist2coast[on_land]<-pts$dist2coast[on_land]*-1 #points on land have negative distance!
pts<-spTransform(pts,ANT_POL_STEREO)

dist2coast <- SpatialPixelsDataFrame(pts,data.frame(dist2coast_km=pts$dist2coast))
dist2shelf <- SpatialPixelsDataFrame(pts,data.frame(dist2shelf_km=pts$dist2shelf))
shelf_habitat <- SpatialPixelsDataFrame(pts,data.frame(on_off_shelf=pts$on_shelf))

crs(dist2coast)<-crs(dist2shelf)<-crs(shelf_habitat)<-crs(pts)

cov<-list(depth=depth,slope=slope,aspect=aspect,dist2coast=dist2coast,dist2shelf = dist2shelf,shelf_habitat = shelf_habitat,
          TRI=TRI,TPI=TPI,roughness=roughness,samplepoints=pts)

for (i in 1:(length(cov))){
  if (names(cov)[i]=='samplepoints'){next}
  r<-raster(cov[[i]])
  name<-names(cov)[i]
  raster::writeRaster(r,file.path(ENVDIR,paste0('env_',name)),format='GTiff',overwrite=T)
}
save(cov, file = file.path(DATRESDIR,'PS112_environmental_covars.RData'),compress='gzip')

