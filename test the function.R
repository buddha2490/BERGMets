rm(list=ls())
library(BERGMets)



# getMetabolites ----------------------------------------------------------

c8 <- getMetabolites("parkinsons_c8")
pos <- getMetabolites("parkinsons_HILICPOS")
neg <- getMetabolites("parkinsons_hilicneg")
blood <- getMetabolites("Validation_blood")
urine <- getMetabolites("validation_urine")
breast <- getMetabolites("breast_metabolomics")
prostate <- getMetabolites("prostate_metabolomics")

summary(breast$biochem$ICC_TECH)
summary(breast$biochem$Prop_Missing)
summary(prostate$biochem$ICC)
summary(prostate$biochem$Prop_Missing)

summary(blood$biochem$ICC_TECH)
summary(blood$biochem$Prop_Missing)


table(duplicated(urine$metabolites$ID))
table(duplicated(urine$ID$ID))
ncol(urine$metabolites)

metabolites <- getMetabolites("prostate_metabolomics")
nrow(metabolites$metabolites)
table(metabolites$ID$caco,useNA="ifany")

id <- metabolites$metabolites$ID[1:10]
comp <- metabolites$biochem$COMP_ID[1:10]

# subset by ID
metabolites2 <- getMetabolites("breast_metabolomics",id)
     dim(metabolites2$metabolites)
     dim(metabolites2$ID)
     dim(metabolites2$biochem)
# subset by COMP_ID
metabolites3 <- getMetabolites("breast_metabolomics",compid=comp)
     dim(metabolites3$metabolites)
     dim(metabolites3$ID)
     dim(metabolites3$biochem)
# subset by both
metabolites4 <- getMetabolites("breast_metabolomics",id,comp)
     dim(metabolites4$metabolites)
     dim(metabolites4$ID)
     dim(metabolites4$biochem)
rm(metabolites,metabolites2,metabolites3,metabolites4)

# Prepare data for testing ------------------------------------------------
src <- file.path("/Users/briancarter/OneDrive - American Cancer Society/Rdata/Metabolomics/BMI Metabolomics - Program Review/data")
load(file.path(src,"final analytic women list - 2018-12-13.rdata"))

breast_metabolomics <- breast$scaledmet
survey <- breast$survey

df <- left_join(survey,breast_metabolomics,"ID")
comp.id <- names(breast_metabolomics)[-1]
biochem <- breast$biochem[,c("COMP_ID","BIOCHEMICAL","SUPER_PATHWAY","SUB_PATHWAY")]

df$agecat <- gtools::quantcut(df$AGE_INT,4)
levels(df$agecat) <- c("Q1","Q2","Q3","Q4")
table(df$agecat, useNA="ifany")
df$bmicat <- with(df,ifelse(BMI<25,0,ifelse(BMI>30,1,NA)))
df$BMI_SD <- df$BMI/sd(df$BMI, na.rm=T)
table(df$bmicat,useNA="ifany")
df$WAISTCAT <- as.numeric(gtools::quantcut(df$WAIST,2))-1
covars <- c("AGE_INT","DIABETES","CHOLESTEROL","LASTATE","HRTKIND","PHYSACT","ALCOHOL","BPDRUG")



# Association analyses ----------------------------------------------------


# Simple association - linear/logistic
one_lin <- metsAssoc(dat=df, biochem, outcome="BMI_SD", compid=comp.id[1:2], covariates=covars)
one_log <- metsAssoc(dat=df, biochem, outcome="bmicat", compid=comp.id[1:2], covariates=covars)

# stratified models - linear/logistic
two_lin <- metsInt(dat=df,biochem,outcome="BMI",compid=comp.id[1:2],intvar="agecat",
                covariates=covars)
two_log <- metsInt(dat=df,biochem,outcome="bmicat",compid=comp.id[1:2],intvar="agecat",
                covariates=covars)

# interaction models - linear/logistic
three_lin <- metsInt(dat=df,biochem,outcome="BMI",compid=comp.id[1:2],intvar="AGE_INT",
                covariates=covars)
three_log <- metsInt(dat=df,biochem,outcome="bmicat",compid=comp.id[1:2],intvar="AGE_INT",
                covariates=covars)



# metsModels() ------------------------------------------------------------
# Check to see if metsModels() works
# Simple association - linear/logistic
one_lin2 <- metsModels(dat=df, biochem, outcome="BMI_SD", compid=comp.id[1:2], covariates=covars)
one_log2 <- metsModels(dat=df, biochem, outcome="bmicat", compid=comp.id[1:5], covariates=covars)

 # stratified models - linear/logistic
two_lin2 <- metsModels(dat=df,biochem,outcome="BMI",compid=comp.id[1:2],intvar="agecat",
                   covariates=covars)
two_log2 <- metsModels(dat=df,biochem,outcome="bmicat",compid=comp.id[1:2],intvar="agecat",
                   covariates=covars)

# interaction models - linear/logistic
three_lin2 <- metsModels(dat=df,biochem,outcome="BMI",compid=comp.id[1:2],intvar="AGE_INT",
                     covariates=covars)
three_log2 <- metsModels(dat=df,biochem,outcome="bmicat",compid=comp.id[1:2],intvar="AGE_INT",
                     covariates=covars)

identical(one_lin,one_lin2)
identical(one_log,one_log2)
identical(two_lin,two_lin2)
identical(two_log,two_log2)
identical(three_lin,three_lin2)
identical(three_log,three_log2)



# getID() -----------------------------------------------------------------
out <- file.path("/users/briancarter/onedrive - american cancer society/rdata/my r packages/SAS")

# Create an ID file in the global environment
getID("breast_metabolomics")

# SAS xport
getID("breast_metabolomics",out)

# makeSAS
makeSAS("breast_metabolomics",filepath=out)

?makeSAS


# partialMatrix -----------------------------------------------------------

# Full matrix - 781 metabolites - 7.5 minutes
system.time(
mat <- corMatrix(df,comp.id,comp.id,"AGE_INT")
)


# Subset matrix - 20 metabolites - 2.78 seconds
system.time(
partialmat <- corMatrix(df,comp.id[1:20],comp.id[1:20],comp.id[1:20])
)


# Subset - source has more than target
system.time(
     mat <- corMatrix(df, comp.id[1:20], comp.id[1:10], "AGE_INT")
)

# Subset - target has more than source - 3.68sec
system.time(
     mat <- corMatrix(df, comp.id[1:10], comp.id[1:20], "AGE_INT")
)

# simple correlation matrix
mat <- corMatrix(df, comp.id[1:10],comp.id[1:10])

# Correlation matrix, source has more than target
mat <- corMatrix(df, comp.id[1:20],comp.id[1:10])

# Correlation matrix, target has more than source
mat <- corMatrix(df, comp.id[1:10],comp.id[1:20])



# getAUC() ----------------------------------------------------------------
rm(list=ls())
metabolites <- getMetabolites("breast_metabolomics")
df <- dplyr::left_join(survey,metabolites$metabolites,"ID")

comp.id <- metabolites$biochem$COMP_ID[1:10]

# code BMI categories as [0,1] - obese vs normal
df$bmi_binary <- ifelse(df$BMI < 25,0, ifelse(df$BMI >= 30,1,NA))


system.time(myAUC <- getAUC(df,
                outcome="bmi_binary",
                compid = comp.id,
                covariates="AGE_INT",
                normalize=T)
)
head(myAUC)


# heatmaps ----------------------------------------------------------------
breast_metabolomics <- getMetabolites("breast_metabolomics")
biochem <- breast_metabolomics$biochem
metdata <- dplyr::left_join(survey,breast_metabolomics$metabolites,"ID")
comp.id <- sample(biochem$COMP_ID,50,replace=T)
load("C:/Users/briancarter/OneDrive - American Cancer Society/Rdata/Metabolomics/BMI Metabolomics - Program Review/output/List of significant metabolites - 2018-12-13.Rdata")
metabolite.subset <- sig.women[1:50]


png(file.path("test mets heatmap.png"),height=10,width=10,units="in",res=400)


metsHeatmap(dat=metdata,
            compid=metabolite.subset,
            covariates=NULL,
            biochem=biochem,
            method="pearson",
            hclust=T,
            title="Title goes here",
            subtitle="This is where a subtitle will be",
            plot_colors=NULL,
            xaxis_size=7,
            yaxis_size=7,
            title_size=12,
            subtitle_size=8,
            plot_theme=NULL)


dev.off()




# ROC curves --------------------------------------------------------------
bmi <- file.path("C:/Users/briancarter/OneDrive - American Cancer Society/Rdata/Metabolomics/BMI Metabolomics - Program Review/output")
load(file.path(bmi,"list of significant metabolites - 2018-12-13.rdata"))
metdata$BMIBIN <- metdata$BMIBIN-1
out <- file.path("/users/briancarter/onedrive - american cancer society/rdata/my r packages/testoutput")

png(file.path(out,"multi ROC.png"),height=6,width=6,units="in",res=400)
metsROC(dat=df,
        outcome="bmi_binary",
        compid=comp.id[1:5],
        covariates="AGE_INT",
        normalize=T,
        xlabel="1 - Specificity",
        ylabel="Sensitivity",
        title=NULL,
        colors=NULL,
        title_size=10,
        #subtitle_size=12,
        #xaxis_size=12,
        #yaxis_size=12,
        plot_theme=NULL)
dev.off()


png(file.path(out,"single ROC.png"),height=6,width=6,units="in",res=400)
metsROC(dat=df,
        outcome="bmi_binary",
        compid=comp.id[1],
        covariates="AGE_INT",
        normalize=T,
        xlabel="1 - Specificity",
        ylabel="Sensitivity",
        title=NULL,
        colors=NULL,
        title_size=10,
        #subtitle_size=12,
        #xaxis_size=12,
        #yaxis_size=12,
        plot_theme=NULL)
dev.off()
# Help files --------------------------------------------------------------

?breast_metabolomics
?Parkinsons_C8
?Parkinsons_HILICNEG
?Parkinsons_HILICPOS
?prostate_metabolomics
?Validation_Blood
?Validation_Urine



rm(list=ls())


src <- file.path("/Users/briancarter/OneDrive - American Cancer Society/Rdata/Metabolomics/BMI Metabolomics - Program Review/data")
load(file.path(src,"final analytic women list - 2018-12-13.rdata"))

cohort <- breast$survey

metdata <- getMetabolites("breast_metabolomics")
biochem <- metdata$biochem
biochem$exclusions <- with(biochem,ifelse(
        ICC_TECH<0.5,1, ifelse(
                Prop_Missing>=0.9, 2, ifelse(
                        SUPER_PATHWAY== "Unknown",3,0))))
biochem <- filter(biochem, exclusions==0)
metabolites <- dplyr::select(metdata$metabolites, ID,biochem$COMP_ID)
biochem <- dplyr::select(biochem,COMP_ID, BIOCHEMICAL, SUPER_PATHWAY)

mydata <- dplyr::left_join(cohort,metabolites,"ID")
comp.id <- biochem$COMP_ID
mydata$BMI_SD <- mydata$BMI/sd(mydata$BMI)
covars.women <- c("AGE_INT","DIABETES","CHOLESTEROL","LASTATE","HRTKIND","PHYSACT","ALCOHOL","BPDRUG")


v <- mydata[,comp.id]

# Initial association analysis
mymodels <- metsModels(mydata,
                       biochem,
                       "BMI_SD",
                       comp.id,
                       covariates=covars.women,
                       normalize=T)

# select based on Bonferonni P-value
sig.mets <- filter(mymodels,as.numeric(P.val) < 0.05/781)$COMP_ID


mymodels$P2 <- as.numeric(mymodels$P.val)

head(mymodels[,c("P2","P.val","P3")])

mymodels$P3 <- as.numeric(format(mymodels$P2,scientific=T))

summary(as.numeric(sig.mets$P.val))


head(mymodels)

summary(mydata$BMI)

mymodels[mymodels$COMP_ID == "X42459","P.val"]



