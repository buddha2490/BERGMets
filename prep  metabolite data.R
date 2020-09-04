# Get_Metabolites
# This will pull the metabolite files

# What files do we have

# Breast cancer - NCIa-03-16ML
# Prostate case-cohort - ACS0-01-16md+
# Parkinsons
# Validation blood
# Validation urine
# QC analysis

"CPS2_Breast"
"CPS2_Prostate"
"Park_C8"
"Park_HILIC_NEG"
"Park_HILIC_POS"
"Validation_blood"
"Validation_urine"



cps2 <- file.path("S:/CPS/BIOSPECIMENS/Metabolomics")
load(file.path(cps2,"breast cancer - ncia-03-16ml/original data/Cleaned breast metabolomics file - 10Nov2017.rdata"))
breast_missing <- readRDS(file.path(cps2,"breast cancer - ncia-03-16ml/original data/N missing for each metabolite - 17June2019.rds"))
breast <- metabolites; rm(metabolites)

load(file.path(cps2,"prostate cancer - acs0-01-16md+/original data/Cleaned prostate metabolomics file - 2018-09-24.rdata"))
prostate_missing <- readRDS(file.path(cps2,"prostate cancer - acs0-01-16md+/original data/N missing for each metabolite - 17June2019.rds"))
prostate <- metabolites; rm(metabolites)


# Breast cancer file ------------------------------------------------------
scaledmet <- breast$scaledmet_noexclusions
ICC <- breast$ICC[,c("COMP_ID","ICC_TECH","CV")]
ICC$ICC <- ICC$ICC_TECH
ICC <- ICC[,names(ICC) != "ICC_TECH"]
ID <- breast$ID.conv[,c("ID","subjid","origID")]
caco <- breast$caco
biochem <- breast$biochem
missing <- breast_missing
missing$COMP_ID <- paste0("X",missing$COMP_ID)

biochem <- Reduce(function(x,y) full_join(x,y,"COMP_ID"), list(biochem,ICC,missing))
biochem <- dplyr::rename(biochem,N_Missing=missing, Prop_Missing=prop.miss)
rm(ICC,missing,breast_missing)

names(caco) <- c("subjid","matchID","caco")
ID <- right_join(ID,caco,"subjid")
rm(caco)
breast_metabolomics <- list(metabolites=scaledmet,ID=ID,biochem=biochem)
rm(breast,ID,scaledmet,biochem)
format(object.size(breast_metabolomics),"Mb")  # 17.8Mb


# Prostate file -----------------------------------------------------------
scaledmet <- prostate$scaledmet_noexclusions
biochem <- prostate$biochem
ICC <- prostate$ICC
caco <- prostate$caco[,c("ID","Overall_group","PrCa")]
names(caco) <- c("ID","CACO","Subgroup")


biochem <- Reduce(function(x,y) full_join(x,y,"COMP_ID"), list(biochem,ICC,prostate_missing))
biochem <- dplyr::rename(biochem,N_Missing=missing, Prop_Missing=prop.miss)

ID <- prostate$ID.conv[,c("ID","subjid","CURRENT_LABEL")]
     names(ID) <- c("ID","subjid","origID")
ID <- ID[!duplicated(ID),]
ID <- full_join(ID,caco,"ID")
rm(caco)

prostate_metabolomics <- list(metabolites=scaledmet,ID=ID,biochem=biochem)
rm(prostate,ID,scaledmet,biochem,prostate_missing,ICC)
format(object.size(prostate_metabolomics),"Mb")  # 6.5Mb



# Fix the metabolite names for breast and prostate cancer -----------------
load(file.path(cps2,"Alternative biochemical names and classes - Dec2018.Rdata"))

alt.breast <- altnames_Dec2018$breast_cancer[,c("COMP_ID","Alt_Name","Metabolite_Class")]
alt.prostate <- altnames_Dec2018$prostate_cancer[,c("COMP_ID","Alt_Name","Metabolite_Class")]

breast_metabolomics$biochem <- left_join(breast_metabolomics$biochem,
                                         alt.breast,
                                         "COMP_ID")
breast_metabolomics$biochem$BIOCHEMICAL <- with(breast_metabolomics$biochem,
                                                ifelse(!is.na(Alt_Name),Alt_Name,BIOCHEMICAL))


prostate_metabolomics$biochem <- left_join(prostate_metabolomics$biochem,
                                         alt.prostate,
                                         "COMP_ID")
prostate_metabolomics$biochem$BIOCHEMICAL <- with(prostate_metabolomics$biochem,
                                                ifelse(!is.na(Alt_Name),Alt_Name,BIOCHEMICAL))


breast_metabolomics$biochem <- dplyr::select(breast_metabolomics$biochem,-Alt_Name)
prostate_metabolomics$biochem <- dplyr::select(prostate_metabolomics$biochem,-Alt_Name)
table(is.na(breast_metabolomics$biochem$BIOCHEMICAL))
table(is.na(prostate_metabolomics$biochem$BIOCHEMICAL))


# Parkinsons data ---------------------------------------------------------
park <- file.path(cps2,"Parkinsons disease")
park.orig <- file.path(park,"Original data")
x <- list.files(park.orig,"Rdata")
for (i in 1:length(x)) load(file.path(park.orig,x[i]))
rm(x,i)

c8 <- C8_Pos; rm(C8_Pos)
hneg <- HILIC_Neg; rm(HILIC_Neg)
hpos <- HILIC_Pos; rm(HILIC_Pos)

## c8
met <- c8$NONQC
ID <- c8$samples[is.na(c8$samples$QC),c("ID","subjid","Sample_ID")]
biochem <- c8$biochem
biochem <- dplyr::rename(biochem,N_Missing=n.miss)
biochem$Prop_Missing <- biochem$N_Missing / nrow(biochem)
biochem <- biochem[,c("Compound","BIOCHEMICAL","Metabolite","COMP_ID","Method","M_Z","RT",
                      "HMDB_ID","Order","N_Missing","Prop_Missing","CV")]
Parkinsons_C8 <- list(metabolites=met,ID=ID,biochem=biochem)
rm(met,ID,biochem,c8)
format(object.size(Parkinsons_C8),"Mb")  # 75.8Mb


# HILIC_Pos
met <- hpos$NONQC
ID <- hpos$samples[is.na(hpos$samples$QC),c("ID","Sample_ID","subjid")]
biochem <- hpos$biochem
biochem <- dplyr::rename(biochem,N_Missing=n.miss)
biochem$Prop_Missing <- biochem$N_Missing / nrow(biochem)
biochem <- biochem[,c("Compound","BIOCHEMICAL","Metabolite","COMP_ID","Method","M_Z","RT",
                      "HMDB_ID","Order","N_Missing","Prop_Missing","CV")]
Parkinsons_HILICPOS <- list(metabolites=met, ID=ID, biochem=biochem)
rm(met,ID,biochem,hpos)
format(object.size(Parkinsons_HILICPOS),"Mb")  # 60.3Mb


# HILIC_Neg
met <- hneg$NONQC
ID <- hneg$samples[is.na(hneg$samples$QC),c("ID","Sample_ID","subjid")]
biochem <- hneg$biochem
biochem <- dplyr::rename(biochem,N_Missing=n.miss)
biochem$Prop_Missing <- biochem$N_Missing / nrow(biochem)
biochem <- biochem[,c("Compound","BIOCHEMICAL","Metabolite","COMP_ID","Method","M_Z","RT",
                      "HMDB_ID","Order","N_Missing","Prop_Missing","CV")]
Parkinsons_HILICNEG <- list(metabolites=met, ID=ID, biochem=biochem)
rm(met,ID,biochem,hneg)
format(object.size(Parkinsons_HILICNEG),"Mb")  # 0.7Mb



# CPS3 files --------------------------------------------------------------
cps3 <- file.path("s:/CPS3/biospecimens/metabolomics")


# Blood first -------------------------------------------------------------
load(file.path(cps3,"Validation study blood - ACS0-01-18ML/original data/CPS3 Validation study - Blood metabolomics - 2019-02-11.Rdata"))

met <- CPS3_Validation_Blood$volnormimpdata
biochem <- CPS3_Validation_Blood$biochem
ID <- CPS3_Validation_Blood$samples
QC <- CPS3_Validation_Blood$QC_samples

biochem <- dplyr::rename(biochem,N_Missing=N.Missing, Prop_Missing=Prop.Missing)
ID <- ID[,c("ID","Current_Label","SAMPLE_NAME","Round")]
Validation_Blood <- list(metabolites=met,ID=ID,biochem=biochem)
rm(met,biochem,ID,QC,CPS3_Validation_Blood)
format(object.size(Validation_Blood),"Mb")  # 16.5Mb




# Urine -------------------------------------------------------------------
load(file.path(cps3,"Validation study urine - ACS0-02-18ML/original data/CPS3 Validation study - Urine metabolomics - 2019-02-12.Rdata"))

met <- CPS3_Validation_Urine$OsmoNormImpData
biochem <- CPS3_Validation_Urine$biochem
ID <- CPS3_Validation_Urine$samples[,c("ID","Current_Label","SAMPLE_NAME","Round")]
biochem <- dplyr::rename(biochem,N_Missing=N.Missing, Prop_Missing=Prop.Missing)
Validation_Urine <- list(metabolites=met,ID=ID,biochem=biochem)
rm(met,biochem,ID,CPS3_Validation_Urine)
format(object.size(Validation_Urine),"Mb")  # 18.7Mb


save(breast_metabolomics,file=file.path(cps2,"Breast cancer - NCIA-03-16ML/breast_metabolomics.rdata"), compress=F)
save(prostate_metabolomics,file=file.path(cps2,"Prostate Cancer - ACS0-01-16MD+/prostate_metabolomics.rdata"), compress=F)

save(Parkinsons_C8,file=file.path(park,"Parkinsons_C8.rdata"), compress=F)
save(Parkinsons_HILICNEG,file=file.path(park,"Parkinsons_HILICNEG.rdata"), compress=F)
save(Parkinsons_HILICPOS,file=file.path(park,"Parkinsons_HILICPOS.rdata"), compress=F)


save(Validation_Blood,file=file.path(cps3,"Validation study blood - ACS0-01-18ML/Validation_Blood.rdata"), compress=F)
save(Validation_Urine,file=file.path(cps3,"Validation study urine - ACS0-02-18ML/Validation_Urine.rdata"), compress=F)
rm(list=ls())



# I've changed my mine about data delivery, need to rename the fil --------
# Rather than deliver the data within the function, I'm going to write a functin to read the data
# since the data are formatted identically, may as well name the objects the same
# to work with them more easily.

urine <- file.path("S:/CPS3/Biospecimens/Metabolomics/Validation study urine - ACS0-02-18ML")
blood <- file.path("S:/CPS3/Biospecimens/Metabolomics/Validation study blood - ACS0-01-18ML")
park <- file.path("S:/CPS/BIOSPECIMENS/Metabolomics/Parkinsons disease")
breast <- file.path("S:/CPS/BIOSPECIMENS/Metabolomics/Breast cancer - NCIA-03-16ML")
prostate <- file.path("S:/CPS/BIOSPECIMENS/Metabolomics/Prostate Cancer - ACS0-01-16MD+")




# urine
load(file.path(urine,list.files(urine,"rdata")))
metabolites <- Validation_Urine
save(metabolites, file=file.path(urine,list.files(urine,"rdata")))
rm(Validation_Urine, metabolites, urine)

# blood
load(file.path(blood,list.files(blood,"rdata")))
metabolites <- Validation_Blood
save(metabolites, file=file.path(blood,list.files(blood,"rdata")))
rm(Validation_Blood, metabolites, blood)

# breast
load(file.path(breast,list.files(breast,"rdata")))
metabolites <- breast_metabolomics
save(metabolites, file=file.path(breast,list.files(breast,"rdata")))
rm(breast_metabolomics, metabolites, breast)

# prostate
load(file.path(prostate,list.files(prostate,"rdata")))
metabolites <- prostate_metabolomics
save(metabolites, file=file.path(prostate,list.files(prostate,"rdata")))
rm(prostate_metabolomics, metabolites, prostate)

# parkinsons
files <- list.files(park,"rdata")
for (i in 1:3){
     load(file.path(park,files[i]))
}

metabolites <- Parkinsons_C8
save(metabolites, file=file.path(park,"Parkinsons_C8.Rdata"))
rm(Parkinsons_C8, metabolites)

metabolites <- Parkinsons_HILICPOS
save(metabolites, file=file.path(park,"Parkinsons_HILICPOS.Rdata"))
rm(Parkinsons_HILICPOS, metabolites)

metabolites <- Parkinsons_HILICNEG
save(metabolites, file=file.path(park,"Parkinsons_HILICNEG.Rdata"))
rm(Parkinsons_HILICNEG, metabolites)



