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

ds_covar_table<-sqldf::sqldf('select count(*) as sightings, seastate,subj from "data.fin" where distance < 1750 group by seastate, subj')
openxlsx::write.xlsx(ds_covar_table, file=file.path(AUXDIR,'ds_covars_table.xlsx'))

#### mock flatfile format - only needed for group size estimator
data.fin$Region.Label='A'
data.fin$size<-data.fin$best_number
data.fin$Sample.Label<-data.fin$transect
data.fin$Effort<-data.fin$transect_length_km
data.fin$Area<-1
####

m1<-ds(data.fin, truncation = truncation, adjustment = NULL, key='hn', formula=~1)
m2<-ds(data.fin, truncation = truncation, adjustment = NULL, key='hn', formula=~seastate)
m3<-ds(data.fin, truncation = truncation, adjustment = NULL, key='hn', formula=~subj)
c1<-ds(data.fin, truncation = truncation, adjustment = 'cos', order = 2, key='hn', formula=~1)
c2<-ds(data.fin, truncation = truncation, adjustment = 'cos', order = 2, key='hn', formula=~seastate)
c3<-ds(data.fin, truncation = truncation, adjustment = 'cos', order = 2, key='hn', formula=~subj)
h1<-ds(data.fin, truncation = truncation, adjustment = NULL, key='hr', formula=~1)
h2<-ds(data.fin, truncation = truncation, adjustment = NULL, key='hr', formula=~seastate)
h3<-ds(data.fin, truncation = truncation, adjustment = NULL, key='hr', formula=~subj)

ds_table<-summarize_ds_models(m1,m2,m3,c1,c2,c3,h1,h2,h3,output='plain', delta_only=F)

for (idx in 1:3){
  modelname<-paste0('m',idx)
  model<-eval(parse(text=modelname))
  
  png(file.path(AUXDIR,paste0('ds_det_fct_plot_',modelname,'.png')),res=300,width=2400,height=2400)
  det.fct.plot(model)
  graphics.off()
  
  modelname<-paste0('h',idx)
  model<-eval(parse(text=modelname))
  
  png(file.path(AUXDIR,paste0('ds_det_fct_plot_',modelname,'.png')),res=300,width=2400,height=2400)
  det.fct.plot(model)
  graphics.off()
  
  modelname<-paste0('c',idx)
  model<-eval(parse(text=modelname))
  
  png(file.path(AUXDIR,paste0('ds_det_fct_plot_',modelname,'.png')),res=300,width=2400,height=2400)
  det.fct.plot(model)
  graphics.off()
}

ds_model<-c1

o<-ds_model$ddf$data
o$esw<-predict(ds_model$ddf,esw=T,newdata=o)$fitted
esw<-mean(unique(o$esw))/1000

dsm_fin<-sqldf::sqldf('select seg_label, 2*count(*) as G, sum(best_number) as I from data where species = "bphy" group by seg_label')
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

predGrid_data<-list(predGrid=predGrid@data,PREDGRID=predGrid,cellsize=cellsize)
save(predGrid_data,file=file.path(DATRESDIR,'PS112_predGrid.RData'),compress='gzip') 

data.dsm<-sqldf::sqldf('select * from dsm_seg as a left join dsm_fin as b on a.seg_label = b.seg_label')
data.dsm<-data.dsm[,-which(names(data.dsm)=='seg_label')[2]] #remove duplicate column seg_label
data.dsm[is.na(data.dsm)]<-0

#### Group size estimation ####
# we use the average group size from the distance sampling model (using the abundance estimate from mock data, see lines 21ff)
gs<-ds_model$dht$Expected.S$Expected.S[1]
gs_se<-ds_model$dht$Expected.S$se.Expected.S[1]

png(file.path(AUXDIR,'ds_group_size_regression.png'),res=300,width=2400,height=2400)
size_regression_plot(data.fin, truncation_width = truncation)
graphics.off()

dsm_data<-list(data=data.dsm,sigs=dsm_fin)
ds_data<-list(data=data.fin,model=ds_model, esw = esw, groups=list(gs=gs,se=gs_se, total_I = sum(data.fin$best_number), total_groups = length(data.fin$best_number)))

save(dsm_data, file=file.path(DATRESDIR,'PS112_dsm_data.RData'),compress='gzip')
save(ds_data, ds_table, ds_model, file=file.path(DATRESDIR,'PS112_ds_data.RData'),compress='gzip')
