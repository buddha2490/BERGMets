# getID

# 1.  breast_metabolomics
# 2.  Parkinsons_C8
# 3.  Parkinsons_HILICNEG
# 4.  Parkinsons_HILICPOS
# 5.  prostate_metabolomics
# 6.  Validation_Blood
# 7.  Validation_Urine

getID <- function(dat,filepath){

     cps2 <- file.path("S:/CPS/BIOSPECIMENS/Metabolomics")
     park <- file.path(cps2,"Parkinsons disease")
     cps3 <- file.path("s:/CPS3/biospecimens/metabolomics")

     if (toupper(dat)=="BREAST_METABOLOMICS"){
          load(file.path(cps2,"Breast cancer - NCIA-03-16ml/breast_metabolomics.rdata"))
          cohort <- "CPS2"
     }
     if (toupper(dat)=="PARKINSONS_C8"){
          load(file.path(park,"Parkinsons_C8.rdata"))
          cohort <- "CPS2"
     }
     if (toupper(dat)=="PARKINSONS_HILICNEG"){
          load(file.path(park,"PARKINSONS_HILICNEG.rdata"))
          cohort <- "CPS2"
     }
     if (toupper(dat)=="PARKINSONS_HILICPOS"){
          load(file.path(park,"PARKINSONS_HILICPOS.rdata"))
          cohort <- "CPS2"
     }
     if (toupper(dat)=="PROSTATE_METABOLOMICS"){
          load(file.path(cps2,"Prostate Cancer - ACS0-01-16MD+/prostate_metabolomics.rdata"))
          cohort <- "CPS2"
     }
     if (toupper(dat)=="VALIDATION_BLOOD"){
          load(file.path(cps3,"Validation study blood - ACS0-01-18ML/Validation_Blood.rdata"))
          cohort <- "CPS3"
     }
     if (toupper(dat)=="VALIDATION_URINE"){
          load(file.path(cps3,"Validation study urine - ACS0-02-18ML/Validation_Urine.rdata"))
          cohort <- "CPS3"
     }

     # Pull ID file
     ID <- metabolites$ID

     # set filename
     filename <- paste0("ID file from ",dat," - ",Sys.Date(),".sas7bdat")

     # Write a SAS file of IDs
     haven::write_sas(ID,file.path(filepath,filename))
     message(paste0("ID file ", filename," has been written to ",filepath))

}
