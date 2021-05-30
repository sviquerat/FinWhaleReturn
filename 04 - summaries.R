rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #only works in RStudio!
SCRIPTDIR<-file.path(getwd(),'SCRIPT')
for (f in list.files(SCRIPTDIR,full.names = T)){ source(f)}

#### SUMMARIES ####
library(mgcv)
library(raster)

load(file.path(DATRESDIR,'PS112_dsm_data.RData'))
load(file.path(DATRESDIR,'PS112_ds_data.RData'))
load(file.path(DATRESDIR,'PS112_gam_data.RData'))

#load(file.path(DATRESDIR,'PS112_spatialData.RData'))
data<-dsm_data$data
data$date<-as.Date(data$date,format="%Y-%M-%D")
dsm_sigs<-sqldf::sqldf('select count(*) as N_sig, min(I) as I_min, max(I) as I_max, avg(I) as I_avg, sum(I) as I_sum, 
                       min(G) as G_min, max(G) as G_max, sum(G) as G_sum, avg(G) as G_avg, avg(I/G) as gs from data where G>0')
dsm_eff<-sqldf::sqldf('select "" as date_min, "" as date_max,count(*) as N_seg, min(effort_km) as effort_km_min, max(effort_km) as effort_km_max, 
             avg(effort_km) as effort_km_avg, sum(effort_km) as effort_km_sum from data')
seg_sum<-cbind(dsm_eff,dsm_sigs)
seg_sum$date_min<-min(data$date)
seg_sum$date_max<-max(data$date)

load(file.path(DATRESDIR,'PS112_modified_survey_data.RData'))
ds_eff<-sqldf::sqldf('select "" as date_min, "" as date_max, "" as transects, 0 as effort_km_avg, sum(track_length_km ) as effort_km_sum from data')
transect_lengths<-sqldf::sqldf('select avg(transect_length_km) as effort_km_avg from data group by transect')[,1]
ds_eff$effort_km_avg<-mean(transect_lengths)
ds_eff$transects<-length(transect_lengths)
ds_eff$date_min<-min(data$date)
ds_eff$date_max<-max(data$date)

data<-ds_data$data
ds_sigs<-sqldf::sqldf('select count(*) as G, min(best_number) as I_min, max(best_number) as I_max, avg(best_number) as I_avg, sum(best_number) as I_sum, 0 as nL, 0 as gs, 0 as gs_se from data')
ds_sigs$gs<-ds_data$groups$gs
ds_sigs$gs_se<-ds_data$groups$gs_se
ds_sum<-cbind(ds_eff,ds_sigs)
ds_sum$nL<-ds_sum$G/ds_sum$effort_km_sum

openxlsx::write.xlsx(seg_sum,file.path(RESDIR,'PS112_summary_segments.xlsx'))
openxlsx::write.xlsx(ds_sum,file.path(RESDIR,'PS112_summary_effort.xlsx'))

getData<-function(N,Dg,D,areaname='survey area'){
  out<-data.frame(name=areaname)
  
  N<-raster(N)
  Dg<-raster(Dg)
  D<-raster(D)
  
  for (ras in c('N','Dg','D')){
    r<-eval(parse(text=ras))
    r<-values(r)
    r<-r[!is.na(r)]
    cells<-length(r)
    out[[ras]]<-round(mean(r),4)
    out[[paste0(ras,'95CI')]]<-paste0(round(quantile(r,c(.05,.95)),4),collapse=' - ')
    
    if (ras=='N'){
      out[[ras]]<-round(sum(r),0)
      out[[paste0(ras,'95CI')]]<-paste0(round(quantile(r,c(.05,.95))*cells,0),collapse=' - ')
    }
  }
  return(out)
}

N<-  file.path(SPRESDIR,'PS112_predicted cell abundance.tif')
Dg<-  file.path(SPRESDIR,'PS112_predicted group density.tif')
D<-  file.path(SPRESDIR,'PS112_predicted density.tif')
summaries<-getData(N,Dg,D,'survey area')

N<-  file.path(SPRESDIR,'PS112_hotspots_predicted cell abundance.tif')
Dg<-  file.path(SPRESDIR,'PS112_hotspots_predicted group density.tif')
D<-  file.path(SPRESDIR,'PS112_hotspots_predicted density.tif')
summaries<-rbind(summaries,getData(N,Dg,D,'hotspots'))
