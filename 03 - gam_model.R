rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #only works in RStudio!
SCRIPTDIR<-file.path(getwd(),'SCRIPT')
for (f in list.files(SCRIPTDIR,full.names = T)){ source(f)}

#### Density surface model #####
library(mgcv)
library(sp)
library(raster)

load(file.path(DATRESDIR,'PS112_dsm_data.RData'))
load(file.path(DATRESDIR,'PS112_ds_data.RData'))
load(file.path(DATRESDIR,'PS112_predGrid.RData'))

predGrid<-predGrid_data$predGrid
data<-dsm_data$data
data$depth<-data$depth
m1<-gam(G~s(x,y), data = data, offset=log(data$effort_km2),family = tw())
m2<-gam(G~s(x,y) + s(depth), data = data, offset=log(data$effort_km2),family = tw())
m3<-gam(G~s(x,y) + s(TPI), data = data, offset=log(data$effort_km2),family = tw())
m4<-gam(G~s(x,y) + s(TRI), data = data, offset=log(data$effort_km2),family = tw())
m5<-gam(G~s(x,y) + s(aspect), data = data, offset=log(data$effort_km2),family = tw())
m6<-gam(G~s(x,y) + s(roughness), data = data, offset=log(data$effort_km2),family = tw())
m7<-gam(G~s(x,y) + s(dist2coast), data = data, offset=log(data$effort_km2),family = tw())
m8<-gam(G~s(x,y) + s(dist2shelf), data = data, offset=log(data$effort_km2),family = tw())

diag<-data.frame()
for (idx in 1:8){
  m<-eval(parse(text=paste0('m',idx)))
  diag<-rbind(diag,data.frame(model=paste0('m',idx),aic=m$aic,gcv=m$gcv.ubre))
}
rownames(diag)<-NULL
diag<-diag[order(diag$aic),]
print(diag)
model<-m1

y<-predict(m1,newdata=predGrid,se.fit=T,type='response')

PG<-predGrid_data$PREDGRID
PG$pG<-predGrid$pG<-as.numeric(y$fit)
PG$pG_se<-predGrid$pG_se<-as.numeric(y$se.fit)
PG$pD<-predGrid$pD<-predGrid$pG*ds_data$groups$gs
PG$pD_se<-predGrid$pD<-predGrid$pG*ds_data$groups$gs
PG$pD_lo<-predGrid$pD_lo<-predGrid$pD-1.96*PG$pD_se
PG$pD_hi<-predGrid$pD_hi<-predGrid$pD+1.96*PG$pD_se
PG$pD_lo[PG$pD_lo<0]<-predGrid$pD_lo[predGrid$pD_lo<0]<-0

LL<-sp::SpatialPoints(cbind(predGrid$x,predGrid$y),ANT_POL_STEREO)
LL<-sp::SpatialPointsDataFrame(LL,predGrid)
r1<-sp::SpatialPixels(LL)
r1<-sp::SpatialPixelsDataFrame(LL,predGrid)
cellArea<-prod(r1@grid@cellsize/1000)
r1$pN<-LL$pN<-r1$pD*cellArea
r1$pN_lo<-LL$pN_lo<-r1$pD_lo*cellArea
r1$pN_hi<-LL$pN_hi<-r1$pD_hi*cellArea
r<-raster::raster(r1)
responses<-list('pD','pG','pN','pN_lo','pN_hi')
names(responses)<-c('predicted density','predicted group density','predicted cell abundance','predicted cell abundance (95 low)','predicted cell abundance (95 high)')
island.threshold<-0.9

hotspot_summary<-data.frame()
rasterList<-c()
for (i in 1:length(responses)){
  resp<-responses[[i]]
  fancyName<-names(responses)[i]
  r_exp<-raster::rasterize(LL,r,resp)
  names(r_exp)<-resp
  q <- quantile(r_exp, island.threshold)
  hotspots <- reclassify(r_exp, cbind(-Inf, q, NA))
  v<-hotspots@data@values
  v<-v[!is.na(v)]
  area<-length(v)*cellArea
  hotspot_summary<-rbind(hotspot_summary,data.frame(ref='hotspots',name=fancyName,area=area,sum=sum(v),avg=mean(v)))
  raster::writeRaster(r_exp,file.path(SPRESDIR,paste0('PS112_',fancyName)),format='GTiff',overwrite =T)
  raster::writeRaster(hotspots,file.path(SPRESDIR,paste0('PS112_hotspots_',fancyName)),format='GTiff',overwrite =T)
  
  rasterList<-c(rasterList,file.path(SPRESDIR,paste0('PS112_',fancyName)),file.path(SPRESDIR,paste0('PS112_hotspots_',fancyName)))
  
  png(file.path(GFXRESDIR,paste0('prediction_',fancyName,'.png')),res=300,width=2400,height=2400)
  plot(r_exp,main=fancyName)
  graphics.off()
  
  png(file.path(GFXRESDIR,paste0('hotspot_',fancyName,'.png')),res=300,width=2400,height=2400)
  plot(hotspots,main=fancyName)
  graphics.off()
  
}
hotspot_summary<-rbind(hotspot_summary,data.frame(ref='survey area',name='Abundance',area=cellArea*nrow(r1),sum=sum(r1$pN),avg=mean(r1$pN)))
hotspot_summary<-rbind(hotspot_summary,data.frame(ref='survey area',name='Density',area=cellArea*nrow(r1),sum=sum(r1$pD),avg=mean(r1$pD)))
openxlsx::write.xlsx(hotspot_summary,file.path(RESDIR,'PS112_summary_predictions.xlsx'))

gam_data<-list(predGrid=predGrid,PREDGRID=PG,response.key=responses, summary = hotspot_summary, rasterList = rasterList)
save(gam_data, file=file.path(DATRESDIR,'PS112_gam_data.RData'),compress='gzip')
