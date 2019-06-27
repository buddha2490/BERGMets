## Pull metabolite data functions



# dat:
# 1.  breast_metabolomics
# 2.  Parkinsons_C8
# 3.  Parkinsons_HILICNEG
# 4.  Parkinsons_HILICPOS
# 5.  prostate_metabolomics
# 6.  Validation_Blood
# 7.  Validation_Urine
getMetabolites <- function(dat, IDlist=NULL, compid=NULL){

cps2 <- file.path("S:/CPS/BIOSPECIMENS/Metabolomics")
park <- file.path(cps2,"Parkinsons disease")
cps3 <- file.path("s:/CPS3/biospecimens/metabolomics")

if (toupper(dat)=="BREAST_METABOLOMICS"){
     load(file.path(cps2,"Breast cancer - NCIA-03-16ml/breast_metabolomics.rdata"))
}
if (toupper(dat)=="PARKINSONS_C8"){
     load(file.path(park,"Parkinsons_C8.rdata"))
}
if (toupper(dat)=="PARKINSONS_HILICNEG"){
     load(file.path(park,"PARKINSONS_HILICNEG.rdata"))
}
if (toupper(dat)=="PARKINSONS_HILICPOS"){
     load(file.path(park,"PARKINSONS_HILICPOS.rdata"))
}
if (toupper(dat)=="PROSTATE_METABOLOMICS"){
     load(file.path(cps2,"Prostate Cancer - ACS0-01-16MD+/prostate_metabolomics.rdata"))
}
if (toupper(dat)=="VALIDATION_BLOOD"){
     load(file.path(cps3,"Validation study blood - ACS0-01-18ML/Validation_Blood.rdata"))
}
if (toupper(dat)=="VALIDATION_URINE"){
     load(file.path(cps3,"Validation study urine - ACS0-02-18ML/Validation_Urine.rdata"))
}
if (!is.null(compid) &
    toupper(dat) %in% c("PARKINSONS_C8","PARKINSONS_HILICNEG","PARKINSONS_HILICPOS")){
     warning("The Parkinsons data are not organized by COMP_ID, all data will be returned")
     compid <- NULL
}
if (!is.null(compid)){
     warning(paste0("You've chosen to subset the data to only include ",length(compid),
                    " metabolites."))
}
if (!is.null(IDlist)){
     warning(paste0("You've chosen to subset the data to only include ", length(IDlist),
                    " samples."))
}

f <- paste0("?",dat)
message(paste0("To access the documentation file for this data set, please search help with the following command: " ,
        f))

# return all the data if subsetting not requested
if (is.null(compid) & is.null(IDlist)){
     return(metabolites)
} else { # otherwise, need to subset
     metdata <- metabolites$metabolites
     ID <- metabolites$ID
     biochem <- metabolites$biochem

     if (!is.null(compid)){
          metdata <- metdata[,names(metdata) %in% c("ID","Round",compid)]
          biochem <- biochem[biochem$COMP_ID %in% compid,]
     }
     if (!is.null(IDlist)){
          metdata <- metdata[metdata$ID %in% IDlist,]
          ID <- ID[ID$ID %in% IDlist, ]
     }
     metabolites <- list(metabolites=metdata, ID=ID, biochem=biochem)
     return(metabolites)
}
}
