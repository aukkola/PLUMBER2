#All functions are from FluxnetLSM but reproduced here to allow post-processing
#of PLUMBER2 inputs

#' Synthesises downward longwave radiation based on Tair and rel humidity
SynthesizeLWdown <- function(TairK, RH, technique='Abramowitz_2012') {
  
  #Three techniques available, see Abramowitz et al. (2012),
  #Geophysical Research Letters, 39, L04808 for details
  
  zeroC <- 273.15
  
  #Inputs missing, set lwdown missing
  if (is.na(TairK) | is.na(RH)) {
    lwdown <- NA
    
    #Else synthesise value
  } else {
    
    if (technique=='Swinbank_1963') {
      # Synthesise LW down from air temperature only:
      lwdown <- 0.0000094*0.0000000567*TairK^6
      
    } else if(technique=='Brutsaert_1975') {
      satvapres <- 611.2*exp(17.67*((TairK-zeroC)/(TairK-29.65)))
      vapres    <- pmax(5,RH)/100*satvapres
      emiss     <- 0.642*(vapres/TairK)^(1/7)
      lwdown    <- emiss*0.0000000567*TairK^4
      
    } else if(technique=='Abramowitz_2012') {
      satvapres <- 611.2*exp(17.67*((TairK-zeroC)/(TairK-29.65)))
      vapres    <- pmax(5,RH)/100*satvapres
      lwdown    <- 2.648*TairK + 0.0346*vapres - 474
      
    } else {
      CheckError('S4: Unknown requested LWdown synthesis technique.')
    }
  }
  
  
  return(lwdown)
}

#-----------------------------------------------------------------------------
