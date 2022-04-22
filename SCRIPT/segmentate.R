#### segmentation script for line transect data
#### splits continuous stretches of effort into discrete segments of user defined length

segmentate<-function(DATA,TRANSECTCOL=NULL,DISTANCECOL=NULL, LIMIT=NULL, EXTRA=NULL, MINROW=1, MIN_SEG_LENGTH=0, MIN_POINT_DISTANCE=NULL,
                     MULTIPLIER=5,COORD_COL=c('lat','lon'),MIDPOINT_METHOD='AVG',QUIET=F){
  DATA[[TRANSECTCOL]]<-factor(DATA[[TRANSECTCOL]])

  if( !is.null( MIN_POINT_DISTANCE ) ){
    MULTIPLIER <- NULL
  }else{
    print('Using automatic estimation of minimum distance between Points')
  }

  if(!is.null(EXTRA)){
    DATA[[EXTRA]]<-factor( DATA[[EXTRA]] )
    out<-.segmentate_extra(DATA,TRANSECTCOL,DISTANCECOL,LIMIT,EXTRA,MINROW,MIN_SEG_LENGTH,MIN_POINT_DISTANCE,MULTIPLIER,QUIET)
  }else{
    out<-.segmentate_normal(DATA,TRANSECTCOL,DISTANCECOL,LIMIT,MINROW,MIN_SEG_LENGTH,MIN_POINT_DISTANCE,MULTIPLIER,QUIET)
  }

  midpoints<-.create_Midpoints(out,COORD_COL,MIDPOINT_METHOD)
  OUT<-list(data=out,centre = midpoints)
  return(OUT)

}

.create_Midpoints<-function(DATA,COORD_COL,MIDPOINT_METHOD){
  midpoints<-data.frame()
  midpoints$seg_label<-NULL
  midpoints[[ COORD_COL[1] ]] <-NULL
  midpoints[[ COORD_COL[2] ]] <-NULL
  outcoords<-paste('seg_',COORD_COL,sep='')

  if( length( which( names(DATA) %in%  outcoords)) == 2 & MIDPOINT_METHOD %in% c('MED','NEAR','AVG')){
    print(paste('calculating ',MIDPOINT_METHOD,' midpoints',sep=''))
    for (seg in unique(DATA$seg_label)){
      d<-subset(DATA, seg_label==seg)
      d<-d[,which(names(d) %in% outcoords)]
      switch(MIDPOINT_METHOD,
             'MED' = X <- .find_Median(d),
             'NEAR' = X <- .find_Nearest(d),
             'AVG' = X <- .find_Average(d)
      )
      midpoints<-rbind(midpoints,cbind(seg,X[1],X[2],MIDPOINT_METHOD))
    }
    names(midpoints) <- c('seg_label',outcoords, 'METHOD')
  }else{midpoints<-NA}
  return(midpoints)
}

.find_Nearest<-function(COORDS){
  X <- c(mean( COORDS[,1], na.rm=T), mean( COORDS[,2], na.rm=T))
  for (r in 1:nrow(COORDS)){
    COORDS$distance[r] <- sqrt(abs(X[1]-COORDS[r,1])^2 + abs(X[1]-COORDS[r,2])^2) #euclidean distance
  }
  X<-COORDS[which(COORDS$distance==min(COORDS$distance)),]
  out<-c(X[,1],X[,2],X[,3])
  return(out)
}

.find_Median<-function(COORDS){
  return( c(median( COORDS[,1], na.rm=T), median( COORDS[,2], na.rm=T)) )
}

.find_Average<-function(COORDS){
  return(c(mean( COORDS[,1], na.rm=T), mean( COORDS[,2], na.rm=T)))
}

.dist_A_B<-function(X0,Y0,X1,Y1){
  DISTANCE<-sqrt(abs(X1-X0)^2 + abs(Y1-Y0)^2) #euclidean distance
  return(DISTANCE) 
}

.segmentate_normal<-function(DATA,TRANSECTCOL,DISTANCECOL,LIMIT,MINROW,MINLENGTH,MINDISTANCE,MULTIPLIER,QUIET){
  shortdata<-data.frame(DATA,label=0,length=0)
  out<-shortdata[-c(1:nrow(shortdata)),]
  for (t in levels(shortdata[[TRANSECTCOL]])){
    d<-shortdata[(shortdata[[TRANSECTCOL]] == t),]
    count<-1
    length<-0
    if(!QUIET){
    print(paste("Transect ",t,sep=''))
    pb<-txtProgressBar(1,nrow(d),initial=0,style=3)
    }
    bar<-0
    for (row in 1:nrow(d)){
      cond1<-(length+d[[DISTANCECOL]][row])>LIMIT #is current segment larger than segment limit?
      cond2<-(nrow(d)-row)>MINROW #is the remainder of the transect more than minimum row numbers?
      cond3<-sum(d[[DISTANCECOL]][row:nrow(d)]) > MINLENGTH #is the remainder of the transect longer than minimum length?
      if (cond1 && cond2 && cond3){
        length<-0
        count <- count+1
      }
      else{
        length<-length+d[[DISTANCECOL]][row]
      }
      d$label[row]<-paste(t,"_S",formatC(count, width = 4, format = "d", flag = "0"),sep='')
      d$length[row]<-length
      bar<-bar+1
      if(!QUIET){setTxtProgressBar(pb,bar)}
    }
    if(!QUIET){close(pb)}
    out<-rbind(out,d)
  }

  out$length<-as.numeric(as.character(out$length))
  names(out)<-paste('seg_',names(out),sep='')
  return(out)
}

.segmentate_extra<-function(DATA,TRANSECTCOL=NULL,DISTANCECOL=NULL,
                            LIMIT=NULL,EXTRA=NULL,MINROW=1,MINLENGTH=0,MINDISTANCE=NULL,MULTIPLIER,QUIET){
  shortdata<-data.frame(DATA,label=0,length=0)
  out<-shortdata[-c(1:nrow(shortdata)),]
  for(E in levels(factor(shortdata[[EXTRA]]))){
    print(paste("Extra: ",E,sep=''))
    shortdata2<-shortdata[(shortdata[[EXTRA]]==E),]
    shortdata2[[TRANSECTCOL]]<-factor(shortdata2[[TRANSECTCOL]])
    for (t in levels(shortdata2[[TRANSECTCOL]])){
      d<-shortdata2[(shortdata2[[TRANSECTCOL]] == t),]
      count<-1
      length<-0
      if(!QUIET){
        print(paste("Transect ",t,sep=''))
      pb<-txtProgressBar(0,nrow(d),initial=0,style=3)
      }
      bar<-0
      for (row in 1:nrow(d)){
        cond1<-(length+d[[DISTANCECOL]][row])>LIMIT #is current segment larger than segment limit?
        cond2<-(nrow(d)-row)>MINROW #is the remainder of the transect more than minimum row numbers?
        cond3<-sum(d[[DISTANCECOL]][row:nrow(d)]) > MINLENGTH #is the remainder of the transect longer than minimum length?
        if (cond1 && cond2 && cond3){#then go to next segment
          length<-0
          count = count+1
        }
        else{
          length<-length+d[[DISTANCECOL]][row]
        }
        d$label[row]<-paste(E,'_',t,"_S",formatC(count, width = 4, format = "d", flag = "0"),sep='')
        d$length[row]<-length
        bar<-bar+1
        if(!QUIET){setTxtProgressBar(pb,bar)}
      }
      if(!QUIET){close(pb)}
      out<-rbind(out,d)
    }
  }
  out$length<-as.numeric(as.character(out$length))
  names(out)<-paste('seg_',names(out),sep='')
  return(out)
}