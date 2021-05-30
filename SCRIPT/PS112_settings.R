cat("\014")  #clear console
set.seed(666)

.requiredPackages<-c('raster','mgcv','sp','rgdal','sqldf','rgeos')
packageInstaller<-function(){ #this works only on windows
  install.packages(setdiff(.requiredPackages, rownames(installed.packages())))
  print('All required packages are installed on System.')
}

packageInstaller()

#### projections ####
WGS84<-sp::CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0')
ANT_POL_STEREO<-sp::CRS('+init=epsg:3031')

#### directories ####
DATADIR<-file.path(getwd(),'DATA')
SPDIR<-file.path(DATADIR,'SPATIAL')
ENVDIR<-file.path(DATADIR,'ENV')

RESDIR<-file.path(getwd(),'RESULTS')
DATRESDIR<-file.path(RESDIR,'DATA')
SPRESDIR<-file.path(RESDIR,'SPATIAL')
GFXRESDIR<-file.path(RESDIR,'GFX')

dir.create(DATRESDIR,recursive=T,showWarnings = F)
dir.create(GFXRESDIR,recursive=T,showWarnings = F)
dir.create(SPRESDIR,recursive=T,showWarnings = F)
dir.create(SPDIR,recursive=T,showWarnings = F)
dir.create(DATADIR,recursive=T,showWarnings = F)
dir.create(ENVDIR,recursive=T,showWarnings = F)