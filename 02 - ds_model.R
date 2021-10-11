rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #only works in RStudio!
SCRIPTDIR<-file.path(getwd(),'SCRIPT')
for (f in list.files(SCRIPTDIR,pattern = '*.R',full.names = T)){ source(f)}

#### Distance Sampling detection function
library(Distance)
library(raster)
library(sqldf)

load(file.path(DATRESDIR,'PS112_modified_survey_data.RData'))
load(file.path(DATRESDIR,'PS112_spatialData.RData'))
data.fin$seastate<-as.factor(data.fin$seastate)
data.fin$subj<-as.factor(data.fin$subj)
truncation<-1750

png(file.path(GFXRESDIR,'group_size_regression.png'),res=300,width=2400,height=2400)
size_regression_plot(data.fin, truncation_width = truncation)
graphics.off()

m1<-ds(data.fin, truncation = truncation, adjustment = NULL, key='hn', formula=~1)
m2<-ds(data.fin, truncation = truncation, adjustment = NULL, key='hn', formula=~seastate)
m3<-ds(data.fin, truncation = truncation, adjustment = NULL, key='hn', formula=~subj)

table<-summarize_ds_models(m1,m2,m3,output='plain')
openxlsx::write.xlsx(table, file=file.path(RESDIR,'Table_2_det_function.xlsx'))

ds_model<-m1

png(file.path(GFXRESDIR,'Figure 2.png'),res=300,width=2400,height=2400)
det.fct.plot(ds_model)
graphics.off()

o<-ds_model$ddf$data
o$esw<-predict(ds_model$ddf,esw=T,newdata=o)$fitted
esw<-mean(unique(o$esw))/1000

dsm_fin<-sqldf::sqldf('select seg_label, count(*) as G, sum(best_number) as I from data where species = "bphy" group by seg_label')
dsm_seg<-sqldf::sqldf('select date,transect,seg_label, max(seg_length_km) as effort_km, avg(x) as x, avg(y) as y from data group by seg_label')
dsm_seg$effort_km2<-2*esw*dsm_seg$effort_km

cellsize=2500
predGrid<-sp::spsample(prediction,type='regular',cellsize=cellsize)
crs(predGrid)<-ANT_POL_STEREO
predGrid<-sp::SpatialPointsDataFrame(predGrid,data.frame(ID=1:length(predGrid)))
predGrid$x<-coordinates(predGrid)[,1]
predGrid$y<-coordinates(predGrid)[,2]

LL<-sp::SpatialPoints(cbind(dsm_seg$x,dsm_seg$y),ANT_POL_STEREO)
envFiles<-list.files(ENVDIR,pattern='*.tif',full.names = T)
envVars<-NULL
for (env in envFiles){
  varName<-strsplit(basename(env),'\\.')[[1]][1]
  varName<-strsplit(varName,'env_')[[1]][2]
  
  if (varName == 'shelf_habitat'){next}
  
  r<-raster::raster(env)
  raster::crs(r)<-ANT_POL_STEREO
  print(varName)
  predGrid[[varName]]<-raster::extract(r,predGrid)
  dsm_seg[[varName]]<-raster::extract(r,LL)
  envVars<-c(envVars,env)

}

#pairs(dsm_seg[,c(5,6,8:15)],pch=16,col='grey',cex=.5)

predGrid_data<-list(predGrid=predGrid@data,PREDGRID=predGrid,cellsize=cellsize)
save(predGrid_data,file=file.path(DATRESDIR,'PS112_predGrid.RData'),compress='gzip') 

data.dsm<-sqldf::sqldf('select * from dsm_seg as a left join dsm_fin as b on a.seg_label = b.seg_label')
data.dsm<-data.dsm[,-which(names(data.dsm)=='seg_label')[2]] #remove duplicate column seg_label
data.dsm[is.na(data.dsm)]<-0
gs<-sum(dsm_fin$I)/sum(dsm_fin$G)
gs_se<-sd(dsm_fin$I/dsm_fin$G)/nrow(data.fin)

dsm_data<-list(data=data.dsm,sigs=dsm_fin)
ds_data<-list(data=data.fin,model=ds_model, esw = esw, groups=list(gs=gs,gs_se=gs_se))

save(dsm_data, file=file.path(DATRESDIR,'PS112_dsm_data.RData'),compress='gzip')
save(ds_data, file=file.path(DATRESDIR,'PS112_ds_data.RData'),compress='gzip')
