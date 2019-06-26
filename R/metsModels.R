metsModels <- function(dat, biochem, outcome, metabolites,intvar=NULL,
                       covariates=NULL, normalize=F){
     if (is.null(intvar)){
          foo <- metsAssoc(dat,biochem,outcome,metabolites,covariates,normalize)
     } else {
          foo <- metsInt(dat,biochem,outcome,metabolites,intvar,
                          covariates,normalize)
     }
     return(foo)
}
