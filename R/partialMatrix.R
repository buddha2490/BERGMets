
partialMatrix <- function(dat,source,target,covariates){


if (anyNA(dat[,source])) {
        warning("Source variables include missing data and will be dropped")
}
if (anyNA(dat[,target])) {
        warning("Target variables include missing data and will be dropped")
}
if (anyNA(dat[,covariates])) {
        warning("Covariates variables include missing data and will be dropped")
}
check <- data.frame(var=NULL,type=NULL)
for (i in 1:length(covariates)){
        foo <- data.frame(var=covariates[i],
                          type=class(dat[[covariates[i]]]),
                          stringsAsFactors=F)
        check <- rbind(check,foo)
}
for (i in 1:nrow(check)){
        if (check$type[i] != "numeric"){
                stop(paste0("Covariate ", check$var[i]," must be numeric.",
                            "If categorical, you must dummy code the variable."))
        }
}


mydata <<- dat
mysource <<- source
mytarget <<- target
mycovars <<- covariates


cl <- parallel::makeCluster(parallel::detectCores(logical=T))
parallel::clusterExport(cl, c("mydata","mysource","mytarget","mycovars"))
     one <- parallel::parLapply(cl, mysource, function(x) {
          foo <- sapply(mytarget, function(y) {
               x <- mydata[[x]]
               y <- mydata[[y]]
               z <- as.matrix(mydata[,mycovars])
               if (x==y) { return(1)} else {
                    return(ppcor::pcor.test(x,y,z)$estimate)
               }
               })
          foo <- data.frame(foo=foo)
          names(foo) <- x
          return(foo)
     })
     parallel::stopCluster(cl)
     mat <- Reduce(function(x,y) cbind(x,y), one)
     mat <- data.frame(COMP_ID=mysource,mat,stringsAsFactors=F)
     names(mat) <- c("COMP_ID",mytarget)
     row.names(mat) <- NULL
     rm(mydata, mysource, mytarget, mycovars, envir=globalenv())
     return(mat)
}
