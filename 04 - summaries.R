rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #only works in RStudio!
SCRIPTDIR<-file.path(getwd(),'SCRIPT')
for (f in list.files(SCRIPTDIR,pattern = '*.R',full.names = T)){ source(f)}

#### SUMMARIES ####
library(mgcv)
library(raster)

load(file.path(DATRESDIR,'PS112_dsm_data.RData'))
load(file.path(DATRESDIR,'PS112_ds_data.RData'))
load(file.path(DATRESDIR,'PS112_gam_data.RData'))

#### SIGHTING DATA SUMMARIES (AUXILIARY) ####
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

openxlsx::write.xlsx(seg_sum,file.path(AUXDIR,'PS112_summary_segments.xlsx'))
openxlsx::write.xlsx(ds_sum,file.path(AUXDIR,'PS112_summary_effort.xlsx'))

#### FIGURES AND TABLES FROM DETECTION FUNCTION MODELLING ####
pretty_ds_table<-ds_table[,1:4]
names(pretty_ds_table)<-c('model','key function', 'covariates', 'Cramér von Mises p')
pretty_ds_table[,4]<-round(pretty_ds_table[,4],4)
pretty_ds_table$covariates<-gsub('~','',pretty_ds_table$covariates)
pretty_ds_table$covariates<-gsub('1','-',pretty_ds_table$covariates)
pretty_ds_table$`p0 ± SE`<-paste0(round(ds_table$`Average detectability`,4), ' ± ', round(ds_table$`se(Average detectability)`,4))
pretty_ds_table$AIC<-round(ds_table$AIC,2)
pretty_ds_table$delta_AIC<-round(ds_table$`Delta AIC`,2)
pretty_ds<-pretty_ds_table[order(pretty_ds_table$model),]
openxlsx::write.xlsx(pretty_ds_table, file=file.path(RESDIR,'Table_2_det_function.xlsx'))

png(file.path(RESDIR,'Figure 2.png'),res=600,width=4000,height=4000)
det.fct.plot(ds_model)
graphics.off()

#### FIGURES AND TABLES FROM ADDITIVE MODELLING ####
summaries<-gam_data$summary
pretty_summary<-summaries[,1:2]
names(pretty_summary)<-c('name','area_km2')
for (resp in c('N','Dg','D')){
  sig.digits<-4
  if (length(grep('N', resp))!=0){
    sig.digits=0
  }
  lo<-summaries[,paste0(resp,'_lo')]
  hi<-summaries[,paste0(resp,'_hi')]
  CI<-paste0(round(lo,sig.digits),' - ',round(hi,sig.digits))
  PRED<-round(summaries[[resp]],sig.digits)
  pretty_summary[[paste0(resp,'_95CI')]]<-paste0(PRED, ' (',CI,')')
}
pretty_summary<-pretty_summary[rev(order(pretty_summary$name)),]
openxlsx::write.xlsx(pretty_summary,file.path(RESDIR,'Table_4_summary_abundance.xlsx'))

diags<-gam_data$gam_diags
diags$aic<-round(diags$aic,2)
diags$gcv<-round(diags$gcv,2)
diags$r_squared<-round(diags$r_squared,2)
diags$dev_expl<-paste0(100*round(diags$dev_expl,4),'%')
diags<-diags[order(diags$model),]
openxlsx::write.xlsx(diags, file=file.path(RESDIR,'Table_3_gam_diagnostics.xlsx'))
