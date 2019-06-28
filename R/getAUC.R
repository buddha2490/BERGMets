
getAUC <- function(dat, outcome, compid, covariates=NULL, normalize=T){

# Data checks
checkoutcome <- lapply(outcome, function(x) {
     f <- data.frame(outcome=x, stringsAsFactors=F)
     f$coding <- ifelse(sum(unique(dat[[x]]),na.rm=T) != 1,"Bad","Good")
     f$Class <- ifelse(class(dat[[x]]) != "numeric", "Bad", "Good")
     f
     }) %>% do.call("rbind",.)

checkoutcome <- checkoutcome[checkoutcome$coding != "Good" | checkoutcome$Class != "Good",]
if (nrow(checkoutcome)>0){
     stop("Check coding on your outcome variables.  They must be NUMERIC and coded [0,1] for a logistic regresison")
}
if (normalize==T){
   message("Metabolite values will be log transformed, mean centered, and normalized prior to modeling")
}

mydf <<- dat
paroutcomes <<- outcome
parcompid <<- compid
parcovars <<- covariates
parnormalize <<- normalize

 cl <- parallel::makeCluster(parallel::detectCores(logical=T))
 parallel::clusterExport(cl, c("mydf","paroutcomes","parcompid","parcovars","parnormalize"))

foo <- parallel::parLapply(cl, paroutcomes, function(y) {
     mydf$OUTCOME <- mydf[[y]]
     b <- lapply(parcompid, function(x){
          tmp <- mydf[!is.na(mydf[[x]]),]
          if (parnormalize==T){
               tmp <- BERGMets::normalizeMets(tmp,x)
          }
          form <- formula(paste0("OUTCOME ~",paste0(c(x,parcovars), collapse="+")))
          fit <- glm(form,data=tmp,family="binomial")
          auc <- pROC::roc(response=fit$model$OUTCOME,
                 predictor=fit$fitted.values,
                 plot=F,
                 auc=T)$auc
          data.frame(COMP_ID=x, Outcome=y, AUC=auc, stringsAsFactors=F)
     })
     m <- Reduce(function(x,y) rbind(x,y), b)
 })
   stopCluster(cl)
   names(foo) <- outcome
rm(mydf, paroutcomes, parcompid, parcovars, parnormalize, envir=globalenv())
if (length(outcome)==1) foo <- foo[[1]]
   return(foo)
}


