rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #only works in RStudio!
SCRIPTDIR<-file.path(getwd(),'SCRIPT')
for (f in list.files(SCRIPTDIR,pattern = '*.R',full.names = T)){ source(f)}

#### Density surface model #####
library(mgcv)
library(sp)
library(raster)

load(file.path(DATRESDIR,'PS112_dsm_data.RData'))
load(file.path(DATRESDIR,'PS112_ds_data.RData'))

data<-dsm_data$data

gam_m1<-gam(G~s(x,y), data = data, offset=log(data$effort_km2),family = tw())
gam_m2<-gam(G~s(x,y) + s(depth), data = data, offset=log(data$effort_km2),family = tw())
gam_m3<-gam(G~s(x,y) + s(TPI), data = data, offset=log(data$effort_km2),family = tw())
gam_m4<-gam(G~s(x,y) + s(TRI), data = data, offset=log(data$effort_km2),family = tw())
gam_m5<-gam(G~s(x,y) + s(aspect), data = data, offset=log(data$effort_km2),family = tw())
gam_m6<-gam(G~s(x,y) + s(roughness), data = data, offset=log(data$effort_km2),family = tw())
gam_m7<-gam(G~s(x,y) + s(dist2coast), data = data, offset=log(data$effort_km2),family = tw())
gam_m8<-gam(G~s(x,y) + s(dist2shelf), data = data, offset=log(data$effort_km2),family = tw())

for (model_index in 1:8){
  model <- eval(parse(text=paste0('gam_m',model_index)))
  b<-mgcViz::getViz(model)
  png(file.path(AUXDIR,paste0('gam_check_gam_0',model_index,'.png')),2400,2400,res=300)  
  p<-mgcViz::check.gamViz(b, type='auto', a.qq=list(CI='quantile', level=0.75))
  print(p)
  graphics.off()
}

diag<-data.frame()
for (idx in 1:8){
  m<-eval(parse(text=paste0('gam_m',idx)))
  summ<-summary(m)
  form<-as.character(m$formula)[3]
  diag<-rbind(diag,data.frame(model=paste0('g',idx),covariates=form, aic=m$aic,gcv=m$gcv.ubre, r_squared = summ$r.sq, dev_expl = summ$dev.expl))
}
rownames(diag)<-NULL
diag<-diag[order(diag$aic),]
gam_diags<-diag

### PREDICTION
model<-gam_m1
load(file.path(DATRESDIR,'PS112_predGrid.RData'))
predGrid<-predGrid_data$predGrid

cellsize<-predGrid_data$cellsize # cellsize in one dimension of the prediction grid cells [meters]
cellArea<-(cellsize/1000)^2 # area per prediction grid cell [km²]

y<-predict(model,newdata=predGrid,se.fit=T,type='response')

PG<-predGrid
PG$pG<-as.numeric(y$fit)
PG$pG_se<-as.numeric(y$se.fit)
PG$pG_lo<-PG$pG-1.96*PG$pG_se
PG$pG_lo[PG$pG_lo<0]<-0
PG$pG_hi<-PG$pG+1.96*PG$pG_se

PG$pD<-PG$pG*ds_data$groups$gs

#se of group size is very small (0.06), therefore we assume se of group size is 0
PG$pD_se<-PG$pG_se

PG$pD_lo<-PG$pG_lo*ds_data$groups$gs
PG$pD_hi<-PG$pG_hi*ds_data$groups$gs

PG$pN<-PG$pD*cellArea
PG$pN_lo<-PG$pD_lo*cellArea
PG$pN_hi<-PG$pD_hi*cellArea

responses<-list('pD','pD_lo','pD_hi','pG','pG_lo','pG_hi','pN','pN_lo','pN_hi','pG_se')
names(responses)<-c('predicted density','predicted density (95 low)','predicted density (95 high)',
                    'predicted group density','predicted group abundance (95 low)',
                    'predicted group abundance (95 high)','predicted cell abundance','predicted cell abundance (95 low)',
                    'predicted cell abundance (95 high)', 'predicted standard error')

# create raster stack of predictions
stack<-list()
for (i in 1:length(responses)){
  resp<-responses[[i]]
  fancyName<-names(responses)[i]
  r_exp<-raster::rasterFromXYZ(cbind(PG$x,PG$y,PG[[resp]]),
                               res=c(predGrid_data$cellsize,predGrid_data$cellsize),crs=ANT_POL_STEREO)
  names(r_exp)<-resp
  stack[[resp]]<-r_exp
}

stack<-raster::stack(stack)
stack$p_CV<-(stack$pG_se*100)/stack$pG
stack$CV_mask<-stack$p_CV
values(stack$CV_mask)[values(stack$p_CV)>100] <- 0
values(stack$CV_mask)[values(stack$p_CV)<=100] <- 1

hotspot.threshold<-0.9 # we want to identify the density /abundance that only 10% of all cells pass
q <- quantile(stack$pD, hotspot.threshold)
stack$hotspot_mask <- reclassify(stack$pD, cbind(-Inf, q, NA))

values(stack$hotspot_mask)[!is.na(values(stack$hotspot_mask))]<-1
values(stack$hotspot_mask)[is.na(values(stack$hotspot_mask))]<-0

for (i in 1:length(responses)){
  resp<-responses[[i]]
  fancyName<-names(responses)[i]
  r<-stack[[resp]]
  raster::writeRaster(r,file.path(SPRESDIR,paste0('PS112_',fancyName)),format='GTiff',overwrite =T)
  r_hot<-r*stack$hotspot_mask
  raster::writeRaster(r_hot,file.path(SPRESDIR,paste0('PS112_hotspots_',fancyName)),format='GTiff',overwrite =T)
}

raster::writeRaster(stack$p_CV,file.path(SPRESDIR,'PS112_CV'),format='GTiff',overwrite =T)
raster::writeRaster(stack$CV_mask,file.path(SPRESDIR,'PS112_CV_mask'),format='GTiff',overwrite =T)
raster::writeRaster(stack$hotspot_mask,file.path(SPRESDIR,'PS112_hotspot_mask'),format='GTiff',overwrite =T)

abundance_summary<-data.frame(ref='hotspots',area=sum(values(stack$hotspot_mask)*cellArea,na.rm=T), 
                              Ng = sum(values(stack$pG*stack$hotspot_mask),na.rm=T),
                              Ng_lo = sum(values(stack$pG_lo*stack$hotspot_mask),na.rm=T),
                              Ng_hi = sum(values(stack$pG_hi*stack$hotspot_mask),na.rm=T),
                              N = sum(values(stack$pN*stack$hotspot_mask),na.rm=T),
                              N_lo = sum(values(stack$pN_lo*stack$hotspot_mask),na.rm=T),
                              N_hi = sum(values(stack$pN_hi*stack$hotspot_mask),na.rm=T))
abundance_summary<-rbind(abundance_summary,
                         data.frame(ref='survey area',area=cellArea*ncell(stack),
                                    Ng = sum(values(stack$pG),na.rm=T),
                                    Ng_lo = sum(values(stack$pG_lo),na.rm=T),
                                    Ng_hi = sum(values(stack$pG_hi),na.rm=T),
                                    N = sum(values(stack$pN),na.rm=T),
                                    N_lo = sum(values(stack$pN_lo),na.rm=T),
                                    N_hi = sum(values(stack$pN_hi),na.rm=T)))
abundance_summary<-rbind(abundance_summary,
                         data.frame(ref='cv area',area=cellArea*sum(values(stack$CV_mask),na.rm=T),
                                    Ng = sum(values(stack$pG*stack$CV_mask),na.rm=T),
                                    Ng_lo = sum(values(stack$pG_lo*stack$CV_mask),na.rm=T),
                                    Ng_hi = sum(values(stack$pG_hi*stack$CV_mask),na.rm=T),
                                    N = sum(values(stack$pN*stack$CV_mask),na.rm=T),
                                    N_lo = sum(values(stack$pN_lo*stack$CV_mask),na.rm=T),
                                    N_hi = sum(values(stack$pN_hi*stack$CV_mask),na.rm=T)))

abundance_summary$Dg<-abundance_summary$Ng/abundance_summary$area
abundance_summary$Dg_lo<-abundance_summary$Ng_lo/abundance_summary$area
abundance_summary$Dg_hi<-abundance_summary$Ng_hi/abundance_summary$area
abundance_summary$D<-abundance_summary$N/abundance_summary$area
abundance_summary$D_lo<-abundance_summary$N_lo/abundance_summary$area
abundance_summary$D_hi<-abundance_summary$N_hi/abundance_summary$area

v<-values(stack$pD)
abundance_summary$dmax<-max(v,na.rm=T)
p_max<-which(v==max(v,na.rm=T))
abundance_summary$dmax_lo<-values(stack$pD_lo)[p_max]
abundance_summary$dmax_hi<-values(stack$pD_hi)[p_max]

####create proper raster stack with hard drive files (previously, data was only in memory)
stack<-list()
for (i in 1:length(responses)){
  resp<-responses[[i]]
  fancyName<-names(responses)[i]
  stack[[resp]]<-raster::raster(file.path(SPRESDIR,paste0('PS112_',fancyName,'.tif')))
}
stack<-raster::stack(stack)  
stack$p_CV<-raster::raster(file.path(SPRESDIR,'PS112_CV.tif'))
stack$hotspot_mask<-raster::raster(file.path(SPRESDIR,'PS112_hotspot_mask.tif'))
stack$CV_mask<-raster::raster(file.path(SPRESDIR,'PS112_CV_mask.tif'))
  
gam_data<-list(predGrid=PG,response.key=responses, summary = abundance_summary, model = model, prediction_stack = stack, gam_diags = gam_diags)
save(gam_data, file=file.path(DATRESDIR,'PS112_gam_data.RData'),compress='gzip')
