summary.MetsResults <- function(dat){
     dat <- Reduce(function(x,y) rbind(x,y),dat)
     dat <- dplyr::right_join(mybiochem, dat, "COMP_ID")
        dat$Estimate <- format(round(dat$Estimate,2), nsmall=2, scientific=F)
        dat$LL <- format(round(dat$LL,2), nsmall=2, scientific=F)
        dat$UL <- format(round(dat$UL,2), nsmall=2, scientific=F)
        dat$Association <- paste0(dat$Estimate, " (",
                                  dat$LL,", ",
                                  dat$UL,")")
     dat$P.val <- format(dat$P,scientific=T, digits=3)
     dat <- dplyr::select(dat, -Estimate, -LL, -UL, -P)
     dat$Association <- sub("^\\s+", "", dat$Association)
     row.names(dat) <- NULL
     return(dat)
}

summary.MetsStrat <- function(results){
        new <- results
        P.Int <- unique(new$P_Int)
        foo <- lapply(split(new,new$columns),function(x){
        dat <- x
        strata <- dat$columns
        dat$Estimate <- format(round(dat$Estimate,2), nsmall=2, scientific=F)
        dat$LL <- format(round(dat$LL,2), nsmall=2, scientific=F)
        dat$UL <- format(round(dat$UL,2), nsmall=2, scientific=F)
        dat$Association <- paste0(dat$Estimate, " (",
                                  dat$LL,", ",
                                  dat$UL,")")
        dat$P.val <- format(dat$P, scientific=T, digits=3)
        dat$P_Int <- format(dat$P, scientific=T, digits=3)
        dat <- dat[,c("COMP_ID","Association","P.val")]
        dat$Association <- sub("^\\s+", "",dat$Association)
        names(dat) <- c("COMP_ID", paste0(strata,":",c("Association","P")))
        row.names(dat) <- NULL
        return(dat)
})
foo <- Reduce(function(x,y) merge(x,y,"COMP_ID",all=T), foo)
foo$P.Interaction <- P.Int
return(foo)
}



summary.MetsInt <- function(results){
        new <- results
        P.Int <- unique(new$P_Int)
        dat <- new[new$columns=="Interaction",]
        dat$Estimate <- format(round(dat$Estimate,2), nsmall=2, scientific=F)
        dat$LL <- format(round(dat$LL,2), nsmall=2, scientific=F)
        dat$UL <- format(round(dat$UL,2), nsmall=2, scientific=F)
        dat$Association <- paste0(dat$Estimate, " (",
                                  dat$LL,", ",
                                  dat$UL,")")
        dat$P.val <- format(dat$P, scientific=T, digits=3)
        dat <- dat[,c("COMP_ID","Association","P.val")]
        dat$Association <- sub("^\\s+", "",dat$Association)
        names(dat) <- c("COMP_ID", c("IntTerm:Association","IntTerm:P"))
        row.names(dat) <- NULL
dat$P.Interaction <- P.Int
return(dat)
}

