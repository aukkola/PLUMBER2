diagnostic_plot <- function(site_code, outdir,  met_file, flux_file)
{
  
  
  ### Variables to plot ###
  
  vars <- c("SWdown", "LWdown", "Precip", "Tair", "Qair", "Psurf", "Wind", "CO2air", 
            "Rnet", "Qle", "Qh", "Qg", "Ebal", 
            "LAI", "LAI_alternative")
  
  #Variable types
  var_type <- c("Met", "Met", "Met", "Met", "Met", "Met", "Met", "Met",
                "Flux", "Flux", "Flux", "Flux", "Flux", 
                "Met", "Met") 
  
  
  ### Get file information ###
  
  #Open file handles
  nc_met  <- nc_open(met_file)
  nc_flux <- nc_open(flux_file)
  
  
  #Get time (used to var lengths etc)
  time <- ncvar_get(nc_flux, "time")
  
  #Get timing information
  timing <- GetTimingNcfile(nc_flux)  
  
  
  #QC information
  QC_measured <- 0
  if(substr(site_code, 1, 3) == "AU-") QC_measured <- append(QC_measured, 10)
  
    
    
  
  ### Set up file ###
  
  outfile <- paste0(outdir, "/", site_code, "_timeseries_plot.png")
  png(outfile, height=1300, width=3500, res=50, pointsize=20)
  
  par(mfcol=c(3,5))
  
  
  
  ### Loop through variables ###
  
  for (v in 1:length(vars)) {
    
    #Get data
    
    #Energy balance
    if (vars[v] == "Ebal") {
      
      rnet <- tryCatch(ncvar_get(nc_flux, "Rnet"), error=function(e) rep(0, length(time)))
      qle  <- tryCatch(ncvar_get(nc_flux, "Qle"), error=function(e) rep(0, length(time)))
      qh   <- tryCatch(ncvar_get(nc_flux, "Qh"), error=function(e) rep(0, length(time)))
      
      #Look for ground heat flux, set to zero if not found
      qg <- tryCatch(ncvar_get(nc_flux, "Qg"), error=function(e) rep(0, length(time)))
      
      #Calculate ratio (Rnet-G) / (Qle + Qh)
      var_data <- (rnet - qg) / (qle + qh)
      
      qc_data <- NA
      
      data_unit <- "-"
      
      
      #Cap values to 0-5
      var_data[var_data > 1.5] <- 1.5
      var_data[var_data < 0]  <- 0
      var_data[which(is.infinite(var_data))] <- 1.5
      
      
      
    #All other variables
    } else {
      
      #Set file for met variable
      if (var_type[v] == "Met") {
        nc <- nc_met
        
        #Set file for  flux variable  
      } else if (var_type[v] == "Flux") { 
        nc <- nc_flux
      }
      
      #Get variable data
      var_data <- tryCatch(ncvar_get(nc, vars[v]), error=function(e) rep(NA, length(time)))
      
      #Get QC data
      qc_data <- tryCatch(ncvar_get(nc, paste0(vars[v], "_qc")), error=function(e) NA)
      
      #Get data units
      data_unit <- tryCatch(ncatt_get(nc, varid=vars[v], "units")$value, error=function(e) NA)
      
    }
    
    #Convert precip and Tair units
    if (vars[v] == "Precip") {
      
      var_data <- var_data * (time[2] - time[1])
      data_unit <- paste("mm/", (time[2] - time[1])/60, "min", sep="")
    }
    
    if (vars[v] == "Tair") {
      
      var_data <- var_data - 273.15
      data_unit <- "deg C"
    }
    
    
    
    #Extract QC data and replace all gap-filled values with 0
    # and measured with 1 (opposite to Fluxnet but what PALS expects)
    if (any(!is.na(qc_data))) {
      
      
      qc_data[!(qc_data %in% QC_measured)]  <- 2 #replace gap-filled values with a temporary value
      qc_data[qc_data %in% QC_measured]     <- 1 #set measured to 1
      qc_data[qc_data == 2] <- 0 #set gap-filled to 0
      
      #If first value missing, set to measured (to avoid an error when PALS
      #checks if first value -1)
      if(is.na(qc_data[1])){ qc_data[1] <- 0}
      
      
      #Else set to PALS option corresponding to no QC data
    } else {
      qc_data <- matrix(-1, nrow = 1, ncol = 1)
    }
    
    
    
    #Plot file
    Timeseries(obslabel=site_code, tsdata=as.matrix(var_data),
               varname=vars[v],
               ytext=paste(vars[v], " (", data_unit, ")", sep=""), 
               legendtext=vars[v],
               plotcex=2 , timing=timing, 
               smoothed = FALSE, winsize = 1, 
               plotcolours="black",
               vqcdata = as.matrix(qc_data),
               na.rm=TRUE)
    
    
    
  } #vars
    
  dev.off()
  

  #Close nc handles
  nc_close(nc_met)
  nc_close(nc_flux)
  
  
} #function



#------------------------------------------

#Plotting function, from FluxnetLSM/PALS

# Plots a smoothed and non-smoothed time series
Timeseries <- function(obslabel, tsdata, varname, ytext, legendtext,
                       plotcex, timing, smoothed=FALSE, winsize=1, plotcolours, modlabel='no',
                       vqcdata=matrix(-1,nrow=1,ncol=1), na.rm=FALSE){
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
           labels=paste('Score_smooth: ',sscorestring,sep=''),pos=4)
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







