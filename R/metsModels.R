metsModels <- function(dat, biochem, outcome, compid,intvar=NULL,
                       covariates=NULL, normalize=T){
     if (is.null(intvar)){
          foo <- metsAssoc(dat,biochem,outcome,compid,covariates,normalize)
     } else {
          foo <- metsInt(dat,biochem,outcome,compid,intvar,
                          covariates,normalize)
     }
     return(foo)
}
