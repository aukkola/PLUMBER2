
# Plotting.R
#
# author: Anna Ukkola UNSW 2017
# Function comes from FluxnetLSM R package

#-----------------------------------------------------------------------------

#' Plots a smoothed and non-smoothed time seris
Timeseries <- function(obslabel,tsdata,varname,ytext,legendtext,
                       plotcex,timing,smoothed=FALSE,winsize=1,plotcolours,modlabel='no',
                       vqcdata=matrix(-1,nrow=1,ncol=1),na.rm=FALSE){
  #
  errtext = 'ok'
  metrics = list()
  ncurves = length(tsdata[1,]) # Number of curves in final plot:
  ntsteps = length(tsdata[,1]) # Number of timesteps in data:
  tstepinday=86400/timing$tstepsize # number of time steps in a day
  ndays = ntsteps/tstepinday # number of days in data set
  nyears=as.integer(ndays/365) # find # years in data set
  # x-axis labels:
  xxat=c()
  xxlab=c()
  data_smooth = matrix(NA,(ndays-winsize-1),ncurves) # init
  if(smoothed){
    for(p in 1:ncurves){
      # Reshape into column timesteps, row days:
      data_days=matrix(tsdata[,p],ncol=tstepinday,byrow=TRUE) 
      for(i in 1:(ndays-winsize-1)){
        # Find evaporative fraction using averaging window:
        data_smooth[i,p] = mean(data_days[i:(i+winsize-1),],na.rm=na.rm)
      }
      if(p==1){
        yvalmin = as.character(signif(min(tsdata[,p],na.rm=na.rm),3))
        yvalmax = as.character(signif(max(tsdata[,p],na.rm=na.rm),3))
        datamean = as.character(signif(mean(tsdata[,p],na.rm=na.rm),3))
        datasd = as.character(signif(sd(tsdata[,p],na.rm=na.rm),3))
        
      }else{
        yvalmin = paste(yvalmin,', ',as.character(signif(min(tsdata[,p],na.rm=na.rm),3)),sep='')
        yvalmax = paste(yvalmax,', ',as.character(signif(max(tsdata[,p],na.rm=na.rm),3)),sep='')
        datamean = paste(datamean,', ',as.character(signif(mean(tsdata[,p],na.rm=na.rm),3)),sep='')
        datasd = paste(datasd,', ',as.character(signif(sd(tsdata[,p]),3),na.rm=na.rm),sep='')
      }
    }
    ymin = signif(min(data_smooth,na.rm=na.rm),3)
    ymax = signif(max(data_smooth,na.rm=na.rm),3)
    # If we're adding a gap-filling QC line, make space for it:
    if(vqcdata[1,1] != -1) {
      ymin = ymin - (ymax-ymin)*0.06
    }
    #If ignoring NA, make space for printing % missing
    #Also shift other labels and legend down in this case
    y_adj=1
    if(na.rm){
      ymax=ymax*1.1
      y_adj = 0.94
    }
    xmin = 1
    xmax = length(data_smooth[,1])
    xloc=c(1:xmax)
    # Draw plot:
    #All missing, plot empty
    if(all(is.na(data_smooth[,1]))){
      plot(xloc,xloc,type="n",ylab=ytext, xaxt='n',cex.lab=plotcex,cex.axis=plotcex,xlab='',mgp = c(2.5+plotcex*0.9,0.8,0))
      mtext(side=3, "All values missing", col="red", line=-4)
    } else {
      plot(xloc,data_smooth[,1],type="l",ylab=ytext,lwd=3,
           col=plotcolours[1],ylim=c((ymin),(ymin + (ymax-ymin)*1.2)),
           xaxt='n',cex.lab=plotcex,cex.axis=plotcex,xlab='',mgp = c(2.5+plotcex*0.9,0.8,0))
    }
    # Calculate NME scores:
    if(ncurves>1){
      smoothscore = c()
      allscore = c()
      for(p in 2:ncurves){ # for each additional curve
        lines(data_smooth[,p],lwd=3,col=plotcolours[p])
        smoothscore[p-1] = sum(abs(data_smooth[,1] - data_smooth[,p]))/
          sum(abs(mean(data_smooth[,1]) - data_smooth[,1]))
        allscore[p-1] = sum(abs(tsdata[,1] - tsdata[,p]))/
          sum(abs(mean(tsdata[,1]) - tsdata[,1]))
      }
      # Report NME metric:
      metricname = paste('NME',winsize,'day',sep='')
      if(ncurves==2){ # model only
        metrics[[1]] = list(name='Bias',model_value=mean(tsdata[,2]-tsdata[,1],na.rm=TRUE)) 
        metrics[[2]] = list(name='NME',model_value=allscore[1]) 
        metrics[[3]] = list(name=metricname,model_value=smoothscore[1]) 
      }else if(ncurves==3){
        metrics[[1]] = list(name='Bias',model_value=mean(tsdata[,2]-tsdata[,1],na.rm=TRUE),
                            bench_value=list(bench1=mean(tsdata[,3]-tsdata[,1],na.rm=TRUE) ))
        metrics[[2]] = list(name='NME',model_value=allscore[1],bench_value=list(bench1=allscore[2]))	
        metrics[[3]] = list(name=metricname,model_value=smoothscore[1],
                            bench_value=list(bench1=smoothscore[2]))	
      }else if(ncurves==4){
        metrics[[1]] = list(name='Bias',model_value=mean(tsdata[,2]-tsdata[,1],na.rm=TRUE),
                            bench_value=list(bench1=mean(tsdata[,3]-tsdata[,1],na.rm=TRUE),
                                             bench2=mean(tsdata[,4]-tsdata[,1],na.rm=TRUE) ))
        metrics[[2]] = list(name='NME',model_value=allscore[1],
                            bench_value=list(bench1=allscore[2],bench2=allscore[3]))
        metrics[[3]] = list(name=metricname,model_value=smoothscore[1],
                            bench_value=list(bench1=smoothscore[2],bench2=smoothscore[3]))	
      }else if(ncurves==5){
        metrics[[1]] = list(name='Bias',model_value=mean(tsdata[,2]-tsdata[,1],na.rm=TRUE),
                            bench_value=list(bench1=mean(tsdata[,3]-tsdata[,1],na.rm=TRUE),
                                             bench2=mean(tsdata[,4]-tsdata[,1],na.rm=TRUE),
                                             bench3=mean(tsdata[,5]-tsdata[,1],na.rm=TRUE) ))
        metrics[[2]] = list(name='NME',model_value=allscore[1],
                            bench_value=list(bench1=allscore[2],bench2=allscore[3],bench3=allscore[4]))
        metrics[[3]] = list(name=metricname,model_value=smoothscore[1],
                            bench_value=list(bench1=smoothscore[2],bench2=smoothscore[3],bench3=smoothscore[4]))	
      }
    }
    for(l in 1:nyears){
      xxat[(2*l-1)] = (l-1)*365 + 1
      xxat[(2*l)] = (l-1)*365 + 152	
      xxlab[(2*l-1)]=paste('1 Jan',substr(as.character(timing$syear+l-1),3,4))
      xxlab[(2*l)]=paste('1 Jun',substr(as.character(timing$syear+l-1),3,4))
    }
    # place legend:
    legend(xmin-(xmax-xmin)*0.03,(ymin + (ymax-ymin)*(y_adj+0.24)),legend=legendtext[1:ncurves],lty=1,
           col=plotcolours[1:ncurves],lwd=3,bty="n",cex=max((plotcex*0.75),1))
    # Add title:
    if(modlabel=='no'){
      title(paste('Smoothed ',varname[1],': ',winsize,'-day running mean.   Obs - ',
                  obslabel,sep=''),cex.main=plotcex)
    }else{
      title(paste('Smoothed ',varname[1],': ',winsize,'-day running mean.   Obs - ',
                  obslabel,'	 Model - ',modlabel,sep=''),cex.main=plotcex)
    }
    # Add Max/Min/Mean/SD numbers:
    text(x=(xmin+(xmax-xmin)*0.25),y=c(ymin + (ymax-ymin)*(y_adj+0.19),ymin + (ymax-ymin)*(y_adj+0.14)),
         labels=c(paste('Min = (',yvalmin,')',sep=''),
                  paste('Max = (',yvalmax,')',sep='')),
         cex=max((plotcex*0.75),1),pos=4)
    text(x=(xmin+(xmax-xmin)*0.25),y=c(ymin + (ymax-ymin)*(y_adj+0.09),ymin + (ymax-ymin)*(y_adj+0.04)),
         labels=c(paste('Mean = (',datamean,')',sep=''),paste('SD = (',datasd,')',sep='')),
         cex=max((plotcex*0.75),1),pos=4)
    # Add NME scores to plot (if there is at least one model output):
    if(ncurves>1){
      sscorestring = paste(signif(smoothscore,digits=3),collapse=', ')
      ascorestring = paste(signif(allscore,digits=3),collapse=', ')
      text(x=c(xmin+(xmax-xmin)*0.65),y=(ymin + (ymax-ymin)*(y_adj+0.18)),
           labels=paste('Score_smooth: ',sscorestring,sep=''),,pos=4)
      text(x=c(xmin+(xmax-xmin)*0.65),y=(ymin + (ymax-ymin)*(y_adj+0.12)),
           labels= paste('Score_all: ',ascorestring,sep=''),pos=4)
      text(x=c(xmin+(xmax-xmin)*0.65),y=(ymin + (ymax-ymin)*(y_adj+0.06)),
           labels=' (NME)',pos=4)
    }
    #Print percentage of data missing if na.rm=TRUE and some data missing
    if(na.rm){
      perc_missing = signif(sapply(1:ncol(tsdata), function(x) 
        sum(is.na(tsdata[,x]))/length(tsdata[,x])), digits=3)	   
      if(any(perc_missing > 0)){
        text(xmin-(xmax-xmin)*0.03,y=(ymin + (ymax-ymin)*(y_adj+0.24)),
             paste("(",paste(perc_missing,collapse=", "), ")% data missing", sep=""),
             pos=4,offset=1, col="red", cex=plotcex)
      }
    }
    # Calculate QC time series information, if it exists:
    if(vqcdata[1,1] != -1){
      qcliney = ymin - (ymax-ymin)*0.015# y-location of qc line
      qctexty = ymin + (ymax-ymin)*0.02 # y-location of qc text
      qcpc = signif((1-mean(vqcdata[,1]))*100,2) # % of data that's gapfilled
      # Construct line-plottable version of qc timeseries:
      origline =	qcliney/(vqcdata[,1]) # 0s will become 'Inf'
      gapline = (qcliney/(vqcdata[,1]-1))*-1 # 1s will become 'Inf'
      # Plot qc time series line:
      xloc_qc = c(1:length(origline))/length(origline) * length(xloc)
      lines(xloc_qc,origline,lwd=20,col='gray80', lend=1)
      lines(xloc_qc,gapline,lwd=10,col='indianred', lend=1)
      text(x=xmin,y=qctexty,cex=max((plotcex*0.75),0.85),pos=4,
           labels=paste(qcpc,'% of observed ',varname[1],' is gap-filled:',sep=''))
    }	
    
    #Not smoothed
  }else{
    xmin = 1
    xmax = ntsteps
    xloc=c(1:xmax)
    y_adj=1
    
    #All missing
    if(all(is.na(tsdata[,1]))){
      plot(xloc,xloc,type="n",ylab=ytext,lwd=3,
           yaxt="n", xaxt='n',cex.lab=plotcex,cex.axis=plotcex,xlab='')
      mtext(side=3, "All values missing", col="red", line=-4,mgp = c(2.5+plotcex*0.9,0.8,0))
      #Else plot
    } else {
      # this code not functioning but kept for future modification:
      yvalmin = signif(min(tsdata, na.rm=na.rm),3)
      yvalmax = signif(max(tsdata, na.rm=na.rm),3)
      datamean = signif(mean(tsdata[,1], na.rm=na.rm),3)
      datasd = signif(sd(tsdata[,1], na.rm=na.rm),3)
      ymin = yvalmin
      ymax = yvalmax
      
      #If ignoring NA, make space for printing % missing
      #Also shift other labels and legend down in this case
      if(na.rm){
        ymax=ymax*1.1
        y_adj = 0.94
      }
      
      plot(xloc,tsdata[,1],type="l",ylab=ytext,lwd=3,
           col=plotcolours[1],ylim=c(ymin,(ymin + (ymax-ymin)*1.3)),
           xaxt='n',cex.lab=plotcex,cex.axis=plotcex,xlab='',mgp = c(2.5+plotcex*0.9,0.8,0))
      # Add smoothed curve over whole timeseries:
      data_days=matrix(tsdata[,1],ncol=tstepinday,byrow=TRUE) 
      data_smooth = c()
      dayssmooth = 30
      for(i in 1:(ndays-dayssmooth-1)){
        # Find evaporative fraction using averaging window:
        data_smooth[i] = mean(data_days[i:(i+dayssmooth-1),], na.rm=na.rm)
      }
      xct = c(1:(ndays-dayssmooth-1))
      xsmooth = xct*tstepinday + (tstepinday*dayssmooth / 2 - tstepinday)
      lines(xsmooth,data_smooth,lwd=3,col='gray')
      
      if(ncurves>1){
        for(p in 2:ncurves){ # for each additional curve
          lines(tsdata[,p],lwd=3,col=plotcolours[p])
        }  
      }
      
      legend(0-(xmax-xmin)*0.05,(ymin + (ymax-ymin)*(y_adj+0.42)),legend=legendtext[1:ncurves],lty=1,
             col=plotcolours[1:ncurves],lwd=3,bty="n",cex=max((plotcex*0.75),1))
      # Locations of max,min,mean,sd text:
      stattextx = c(xmin,xmin+(xmax-xmin)*0.5)
      stattexty = c(ymin + (ymax-ymin)*(y_adj+0.18),ymin + (ymax-ymin)*(y_adj+0.24))
      # Write max,min,mean,sd to plot in two lines:
      text(x=stattextx,y=stattexty[2],
           labels=c(paste('Min = ',ymin,sep=''),paste('Max = ',ymax,sep='')),
           cex=max((plotcex*0.75),1),pos=4)
      text(x=stattextx,y=stattexty[1],
           labels=c(paste('Mean = ',datamean,sep=''),paste('SD = ',datasd,sep='')),
           cex=max((plotcex*0.75),1),pos=4)
      #Print percentage of data missing if na.rm=TRUE and some data missing
      if(na.rm){
        perc_missing = signif(sapply(1:ncol(tsdata), function(x) 
          sum(is.na(tsdata[,x]))/length(tsdata[,x])), digits=3)     
        if(any(perc_missing > 0)){
          text((xmax-xmin)*0.5,y=(ymin + (ymax-ymin)*(y_adj+0.42)),
               paste("(",paste(perc_missing,collapse=", "), ")% data missing", sep=""),
               pos=1,offset=1, col="red",cex=plotcex)
        }
        # Calculate QC time series information, if it exists:
        if(vqcdata[1,1] != -1){
          qcliney = ymin + (ymax-ymin)*(y_adj+0.04) # y-location of qc line
          qctexty = ymin + (ymax-ymin)*(y_adj+0.09) # y-location of qc text
          qcpc = signif((1-mean(vqcdata[,1], na.rm=TRUE))*100,2) # % of data that's gapfilled
          # Construct line-plottable version of qc timeseries:
          origline =	qcliney/(vqcdata[,1]) # 0s will become 'Inf'
          gapline = (qcliney/(vqcdata[,1]-1))*-1 # 1s will become 'Inf'
          # Plot qc time series line:
          lines(origline,lwd=20,col='gray80', lend=1)
          lines(gapline,lwd=10,col='red', lend=1)
          text(x=stattextx[1],y=qctexty,cex=max((plotcex*0.75),1),pos=4,
               labels=paste(qcpc,'% of time series is gap-filled:',sep=''))
        }
      }
    } #all NA?
    for(l in 1:nyears){
      xxat[(2*l-1)] = (l-1)*365*tstepinday + 1
      xxat[(2*l)] = (l-1)*365*tstepinday + 183*tstepinday
      xxlab[(2*l-1)]=paste('1 Jan',substr(as.character(timing$syear+l-1),3,4))
      xxlab[(2*l)]=paste('1 Jul',substr(as.character(timing$syear+l-1),3,4))
    }
    title(paste(obslabel,varname[1]),cex.main=plotcex)
    axis(1,at=xxat,labels=xxlab,cex.axis=plotcex,mgp = c(2.3,plotcex*0.7,0))
    result = list(err=FALSE,errtext = errtext,metrics=metrics)
    return(result)
  }
}

#---------------------------
GetTimingNcfile = function(fid){
  # This function gets the time step size, number of timesteps
  # and start date and time details from a netcdf file.
  errtext='ok'
  # Get the name of the time variable in this file
  timevar = FindTimeVarName(fid)
  if(timevar$err){ # Report fatal error
    timing=list(errtext=timevar$errtext,err=TRUE)
    return(timing)
  }
  # Get time units:
  tunits = GetTimeUnits(fid,timevar$name)
  if(tunits$err){ # Report fatal error
    timing=list(errtext=tunits$errtext,err=TRUE)
    return(timing)
  }
  # Get number of time steps:
  ntsteps = GetNumberTimesteps(fid,timevar)
  
  # Get the time step size:
  tstep = GetTimestepSize(fid,timevar$name,tunits,ntsteps)
  if(tstep$err){ # Report fatal error
    timing=list(errtext=tstep$errtext,err=TRUE)
    return(timing)
  }
  
  # Create return list:
  timing = list(err=FALSE,errtext=errtext,tstepsize=tstep$size,tsteps=ntsteps,
                syear=tunits$syear,smonth=tunits$smonth,sdoy=tunits$sdoy,whole=tstep$wholeyear,interval=tstep$interval)
  return(timing)
}

#----

#' Finds time variable name
FindTimeVarName = function(fid){
  # Finds the name of the time variable in an open netcdf file.
  errtext='ok' 
  exists = FALSE # initialise
  nvars = length(fid$var) # number of variables in netcdf file
  ndims = length(fid$dim) # number of dimensions in netcdf file
  for (v in 1:nvars){ # Search through all variables in netcdf file
    if((substr(fid$var[[v]]$name,1,6)=='t_ave_') | # i.e. ORCHIDEE file
       (substr(fid$var[[v]]$name,1,5)=='mscur') | 
       (fid$var[[v]]$name == 'time')){
      # i.e. 	ORCHIDEE file, CLM file or non-dimension variable named time, respectively
      exists = TRUE
      dimvar = FALSE # i.e. time variable is not a dimension variable
      timevarname = fid$var[[v]]$name
      timedimid = fid$var[[v]]$dimids
      break # leave for loop for variables
    }
  }
  if(!exists){ # i.e. none of the above time variables were found
    # Search for time as a dimension variable:
    for (d in 1:ndims){ # Search through all dimensions in netcdf file
      if(fid$dim[[d]]$name=='time' | fid$dim[[d]]$name=='t' |
         fid$dim[[d]]$name=='time_counter' | fid$dim[[d]]$name=='Time'){
        # Now check for time dimension variable:
        if(fid$dim[[d]]$dimvarid$id != -1){ # i.e. dim var exists	
          exists = TRUE # i.e. found time dimension variable
          dimvar = TRUE # time variable is a dimension variable
          timevarname = fid$dim[[d]]$name
          timedimid = fid$dim[[d]]$id
        }else{ # i.e. time dim exists but no dim var
          errtext = paste('T1: Cannot interpret timing in ',stripFilename(fid$filename),
                          ': time dimension exists but no dimension variable.',sep='')
          timevar = list(err=TRUE,errtext=errtext)
          return(timevar)
        }
        break	
      }
    }
    if(!exists){ # Still cannot identify time variable
      # Return to parent function with error:
      errtext = paste('T1: Unable to ascertain name of time variable in', stripFilename(fid$filename))
      timevar = list(err=TRUE,errtext=errtext)
      return(timevar)
    }
  }
  # Return result:
  timevar = list(err=FALSE, errtext=errtext,name=timevarname,dimid=timedimid,dimvar=dimvar)
  return(timevar)
}

#----

#' Gets time units
GetTimeUnits = function(fid,timevarname){
  # Fetches and processes time units from a netcdf file.
  errtext = 'ok'
  if(substr(timevarname,1,5)=='mscur'){ # i.e. CLM file
    # Read date variable:
    date1=as.character(ncvar_get(fid,'mcdate',start=1,count=1))
    syear = as.numeric(substr(date1,1,4))
    smonth = as.numeric(substr(date1,5,6))
    sdoy = as.numeric(substr(date1,7,8))
    units = 'clm_blah'
  }else{
    units = ncatt_get(fid,timevarname,'units')
    if(! units$hasatt){
      errtext = paste('T1: Unable to find time units in', stripFilename(fid$filename))
      tunits = list(err=TRUE,errtext=errtext)
      return(tunits)	
    }
    if(substr(units$value,1,4)=='seco'){ # time units are seconds
      syear = as.numeric(substr(units$value,15,18))
      smonth = as.numeric(substr(units$value,20,21))
      sdoy = as.numeric(substr(units$value,23,24))
      units = 'seconds'
    }else if(substr(units$value,1,4)=='days'){ # time units are days
      syear = as.numeric(substr(units$value,11,14))
      smonth = as.numeric(substr(units$value,16,17))
      sdoy = as.numeric(substr(units$value,19,20))
      units = 'days'
    }else{
      errtext = paste('T1: Unable to interpret time units in', stripFilename(fid$filename))
      tunits = list(err=TRUE,errtext=errtext)
      return(tunits)	
    }
  }
  tunits = list(err=FALSE,errtext=errtext,syear=syear,smonth=smonth,sdoy=sdoy,units=units)
  return(tunits)	
}

#----

#' Gets time step size
GetTimestepSize = function(fid,timevarname,tunits,ntsteps){
  # Fetches time step size
  errtext = 'ok'
  # Read first 2 timesteps of time variable:
  time=ncvar_get(fid,timevarname,start=1,count=2)
  time_end=ncvar_get(fid,timevarname,start=ntsteps,count=1)
  # Define time step size:
  tsize=time[2]-time[1]
  tperiod=time_end - time[1] + tsize
  
  if(tunits$units == 'days'){
    if((tsize >359) && (tsize < 367)){
      interval = 'annual'
      wholeyear = TRUE
    }else if((tsize>27) && (tsize<32)){
      interval = 'monthly'
      if((ntsteps %% 12) == 0){
        wholeyear = TRUE	
      }else{
        wholesyear==FALSE
      }	
    }else if(tsize==1){
      interval = 'daily'
      intyear = Yeardays(tunits$syear,tperiod)
      wholeyear = intyear$whole	
    }else{
      errtext = paste('T1: Unable to interpret time step size in', stripFilename(fid$filename))
      tstep = list(err=TRUE,errtext=errtext)
      return(tstep)	
    }
  }else if(tunits$units == 'seconds'){
    if(tsize <= (3600*3)){ # i.e. less than 3-hourly
      interval = 'timestep'
      tstepinday=86400/tsize # time steps in a day
      ndays = ntsteps/tstepinday # number of days in file
      intyear = Yeardays(tunits$syear,ndays)
      wholeyear = intyear$whole
    }else if(tsize == (3600*24)){
      interval = 'daily'
      intyear = Yeardays(tunits$syear,(tperiod/3600/24))
      wholeyear = intyear$whole
    }else if( (tsize > (27*24*3600)) && (tsize < (32*24*3600)) ){
      interval = 'monthly'
      if((ntsteps %% 12) == 0){
        wholeyear = TRUE	
      }else{
        wholesyear==FALSE
      }	
    }else if( (tsize > (359*24*3600)) && (tsize < (367*24*3600)) ){
      interval = 'annual'
      wholeyear = TRUE
    }else{
      errtext = paste('T1: Unable to interpret time step size in', stripFilename(fid$filename))
      tstep = list(err=TRUE,errtext=errtext)
      return(tstep)	
    }
    
    
    
  }else{
    errtext = paste('T1: Unable to interpret time units in', stripFilename(fid$filename))
    tstep = list(err=TRUE,errtext=errtext)
    return(tstep)
  }
  tstep = list(err=FALSE,errtext=errtext,size=tsize,interval=interval,wholeyear=wholeyear)
  return(tstep)
}

#----

#' Gets the number of time steps
GetNumberTimesteps = function(fid,timevar){		
  # Gets the number of time steps in a netcdf file
  ndims = length(fid$dim)
  # Find out how many time steps there are in the model file:
  for (d in 1:ndims){ # Search through all dimensions in netcdf file
    if(fid$dim[[d]]$id == timevar$dimid){ # i.e. this is the dim of the time variable
      ntsteps = fid$dim[[d]]$len
      break # stop searching for unlim dim
    }
  }
  return(ntsteps)
}

#----

#' Gets the day of month and month, given day of year
doydate = function(doy,leap=FALSE){
  # Doydate returns the day of month and month, given day of year:
  month=getMonthDays(leap)
  # Find month of this doy
  for(m in 1:12){
    if(doy >= month$start[m] && doy < month$start[m+1]){
      doymonth = m
      doyday = doy - month$start[m] + 1
    }
  }
  date = list(month = doymonth, day = doyday)
  return(date)
}

#----

#' Finds number of days per year
Yeardays <- function(startyear,ndays) {
  # Returns: an integer vector of possible number of days in each year of a 
  # dataset, and whether it contains a whole number of years
  if(ndays<365){
    whole=FALSE
    daysperyear = ndays
  }
  daysperyear = c()
  ctr=0 # initialise
  year = startyear # initialise
  days=ndays # initialise
  lpyrs = 0
  # Incrementally remove year of days from total number of days:
  repeat {
    ctr = ctr + 1
    if(is.leap(year)){	
      days = days - 366
      daysperyear[ctr] = 366
      lpyrs = lpyrs + 1
    }else{
      days = days - 365
      daysperyear[ctr] = 365
    }
    year = year + 1
    if(days<365){
      if(days>0 && days!=(365-lpyrs)){ # ie. after removing whole years, days are left over
        daysperyear[ctr+1] = days
        whole=FALSE
      }else if(days==(365-lpyrs)){ # i.e. non leap year data set
        daysperyear[ctr+1] = days
        whole=TRUE
      }else{ # =0
        whole=TRUE
      }
      break
    }
  }
  # Create return list:
  yeardays = list(daysperyear=daysperyear,whole=whole)
  return(yeardays)
}

#----

#' Finds leap years
is.leap = function(year){
  if((((year %% 4)==0) & ((year %% 100)!=0)) || 
     (((year %% 4)==0) & ((year %% 400)==0))){
    leap=TRUE	
  }else{
    leap=FALSE
  }
  return(leap)
}
