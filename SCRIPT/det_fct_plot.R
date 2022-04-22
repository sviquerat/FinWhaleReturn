##### detection function plot for distance sampling models

det.fct.plot<-function(ds.model,COLOURED=F,esw.col='black',...){
  options(warn=-1)
  covar<-gsub('~','',ds.model$ddf$ds$aux$ddfobj$scale$formula)
  o<-ds.model$ddf$data
  o$esw<-predict(ds.model$ddf,esw=T,newdata=o)$fitted
  plot(ds.model,showpoints=F,type='n',pl.col='#F1F1F1',border='darkgrey')
  h<-add_df_covar_line(ds.model,data=o,lty=0)
  h<-apply(h,2,'mean')
  left <- ds.model$ddf$meta.data$left
  width <- ds.model$ddf$meta.data$width
  xx <- seq(left, width, length.out=250)
  lines(h~xx,col=esw.col,lty=1,lwd=1.5)
  if(covar!='1'){
    txt<-NULL
    lvls<-sort(levels(as.factor(ds.model$ddf$data[[covar]])))
    cols<-hcl.colors(length(lvls),palette='Set 2')
    for (i in 1:length(lvls)){
      lvl<-lvls[i]
      col<-cols[i]
      df<-data.frame(covar=lvl)
      names(df)<-covar
      add_df_covar_line(ds.model,data=df,col=col,lwd=1)
      esw<-unique(o$esw[o[[covar]]==lvl])
      N<-sum(ds.model$ddf$data[[covar]]==lvl)
      txt<-c(txt,paste0(covar,': ',lvl,' (N: ',N,', esw: ',round(esw,0),' m)'))
    }
    legend('topright',legend=txt,lty=2,col=cols,lwd=2)
  }
  esw<-mean(unique(o$esw))
  abline(v=esw,col=esw.col,lwd=1.5)
  text(esw,.9,paste(round(esw,0),' m',sep=' '),srt=90,pos=4,col=esw.col)
  options(warn=0)
}

#### size regression plot (group size vs. distance)
size_regression_plot<-function(ds_data, best_col='best_number',distance_col='distance',distance_unit='m',
                          truncation_width=max(ds_data[[distance_col]],na.omit=T)){
 ds_data$distance<-ds_data[[distance_col]]
 ds_data$best_number<-ds_data[[best_col]]
 ds_data<-subset(ds_data, best_number>0)
 ds_data<-subset(ds_data, !is.na(distance))
 m1<-glm(best_number~distance,data=subset(ds_data,distance<=truncation_width),family='poisson')
 m0<-glm(best_number~1,data=subset(ds_data,distance<=truncation_width),family='poisson')
 gs_0<-1+coefficients(m0)[1]
 gs_1<-1+coefficients(m1)[1]
 anv<-anova(m0,m1,test="Chisq")
 p_dissimilarity=1-pchisq( abs(anv$Deviance[2]), abs(anv$Df[2]))

 f0<-paste0('H0: group size ~ ',round(gs_0,4))
 f1<-paste0('H1: group size ~ ',round(gs_1,4),' + ',round(coefficients(m1)[2],4), ' x distance')
 f2<-paste0('Probability of identity: ',round(1-p_dissimilarity,4)*100, ' %')
 
 plot(best_number~distance,data=ds_data,type='n',axes=F,xlab= 'perpendicular distance [m]',ylab='observed group size')
 grid()
 points(best_number~distance,data=subset(ds_data,distance<=truncation_width),pch=16,col=adjustcolor('purple',.6))
 points(best_number~distance,data=subset(ds_data,distance>truncation_width),pch=17,col=adjustcolor('purple',.3))
 p<-data.frame(distance=0:truncation_width)
 p$y<-predict(m1,newdata=p,type='response')
 lines(y~distance,data=p,col=adjustcolor('red',.7))
 abline(h=1+round(coefficients(m0)[1],4),col=adjustcolor('green',.7))
 abline(v=truncation_width,col=adjustcolor('blue',.4))
 legend('topright',legend=c(f0,f1,f2))
 axis(2,at = 1:max(ds_data$best_number))
 axis(4,at = seq(1.25,max(ds_data$best_number)+.25,.5))
 axis(1)
 box()
 return(list(gs_0=gs_0,gs_1=gs_1,p_dissimilarity=p_dissimilarity,anv=anv))
}
