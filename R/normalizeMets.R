normalizeMets <- function(dat, metabolites=NULL){
     df <- dat
     if (!is.null(metabolites)) {
          comp.id <- metabolites
     } else {
          comp.id <- names(df)[-1]
     }

     if (length(comp.id)>1){
     df[,comp.id] <- apply(df[,comp.id],2,as.numeric)
     df[,comp.id] <- apply(df[,comp.id],2,FitAR::glog)
     df[,comp.id] <- apply(df[,comp.id],2, function(x){
          (x - mean(x,na.rm=T) ) / sd(x,na.rm=T)
     })} else {
     df[[comp.id]] <- as.numeric(df[[comp.id]])
     df[[comp.id]] <- FitAR::glog(df[[comp.id]])
     df[[comp.id]] <- (df[[comp.id]] - mean(df[[comp.id]],na.rm=T))/
             sd(df[[comp.id]], na.rm=T)
     }
     return(df)
}
