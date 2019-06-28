
icc9 <- function(y,subject,fixed=NULL,dat,method="ML"){

# Quick install from github
p <-   rownames(installed.packages())
if ("devtools" %in% p == "FALSE") install.packages("devtools")
if ("bootmlm" %in% p == "FALSE") devtools::install_github("marklhc/bootmlm")

data <- dat

# quick function to format my estimates/confidence intervals
formatEstimate <- function(est, lcl, ucl){
    a <- format(round(est,2),nsmall=2)
    b <- format(round(lcl,2),nsmall=2)
    c <- format(round(ucl,2),nsmall=2)
    return(paste0(a, " (",b,",",c,")"))
}

# drop any missings from the dataset
dat <- dat[!is.na(dat[[subject]]),]
dat <- dat[!is.na(dat[[y]]),]

# Define my input variables
SUBJECT <- data[[subject]]
Y <- data[[y]]
if (is.null(fixed)) {FIXED <- NULL} else {FIXED <- as.matrix(data[,fixed])}

# Frequency stuff first
tabx.f <- data.frame(table(SUBJECT))
tabx.p <- data.frame(prop.table(table(SUBJECT)))
tabx <- data.frame(subjid=tabx.f$SUBJECT, COUNT=tabx.f$Freq, PERCENT=tabx.p$Freq*100,
                    stringsAsFactors=F)
rm(tabx.f,tabx.p)


# tab3 lists the number of unique subjects, number of observations, and average obs per subject
tab3 <- data.frame(numw2=nrow(tabx),
                   obsw2=sum(tabx$COUNT,na.rm=T))
tab3$avgrep <- tab3$obsw2/tab3$numw2
rm(tabx)

# Step 1 - run the model
    # Subject ID is always a random effects
    # Y is the response
    # FIXED is optional
    # Method - inputs are "ML" or "REML", default is "ML"

# Sets modelmethod to logical
# the lmer() function has an option "REML=T/F"
# This little bit of code will just fill in that option
if (method == "ML") modelmethod <- F
if (method != "ML") modelmethod <- T

if (!is.null(FIXED)) {
f <- formula(paste0("Y ~ ",paste(c("(1|SUBJECT)","FIXED"),collapse="+")))
    } else {
f <- formula("Y ~ (1|SUBJECT)")
    }
fit <- lme4::lmer(f, REML=modelmethod)  # model object
z <- summary(fit)       # summary object

# Extract the asymptotic covariance matrix
asycov <- data.frame(bootmlm::vcov_vc(fit,F))
names(asycov) <- c("CovP1","CovP2")  # Renaming to match the SAS program
row.names(asycov) <- c("Intercept","Residual")

# Pull out the covariance parameters
estimates <- data.frame(z$varcor)
estimates <- estimates[,c("grp","vcov","sdcor")]
names(estimates) <- c("Parm","Estimate","SD")

# Pull out the fixed effects estimates
fixed_effects <- data.frame(z$coefficients)[1,]
    names(fixed_effects) <- c("Estimate","SE","T")

# Everything I need is in the asycov, random-estimates and fixed-estimates matrices
subject_estimate <- estimates$Estimate[estimates$Parm=="SUBJECT"]
residual_estimate <- estimates$Estimate[estimates$Parm=="Residual"]
residual_covp2 <- asycov["Residual","CovP2"]
intercept_covp1 <- asycov["Intercept","CovP1"]
intercept_covp2 <- asycov["Intercept","CovP2"]
fixed_estimate <- fixed_effects$Estimate
variance_fixed <- fixed_effects$SE^2


# ICC and confidence intervals
ICC_Variance <- (subject_estimate^2 * residual_covp2+residual_estimate^2 * intercept_covp1-2*subject_estimate*residual_estimate*intercept_covp2)/subject_estimate^4
ICC <- estimates$Estimate[estimates$Parm=="SUBJECT"] / sum(estimates$Estimate)
    logICC <- log(ICC / (1-ICC))
    SE_ICC <- sqrt(((ICC^2) * ICC_Variance) / ((1-ICC)^2))
        logLL <- logICC - 1.96*SE_ICC
        logUL <- logICC + 1.96*SE_ICC
    ICCLL <- exp(logLL) / (1+exp(logLL))
    ICCUL <- exp(logUL) / (1+exp(logUL))
finalICC <- formatEstimate(ICC, ICCLL, ICCUL)
rm(ICC_Variance,logICC,ICC,SE_ICC,logLL,logUL,ICCLL,ICCUL)



# Overall CV
Variance_CV <- residual_covp2 / (4*fixed_estimate^2 * residual_estimate) +
    (residual_estimate / fixed_estimate^4) * variance_fixed

CV <- sqrt(residual_estimate) / z$coefficients["(Intercept)","Estimate"]
    logCV <- log(CV)
    varCV <- (1/CV^2) * Variance_CV
    logLL <- logCV - 1.96*sqrt(varCV)
    logUL <- logCV + 1.96*sqrt(varCV)
    CVLL <- exp(logLL)
    CVUL <- exp(logUL)
finalCV <- formatEstimate(CV, CVLL, CVUL)

# Create a final data frame with:
    # n-unique subject
    # n-observations
    # avg-N per subject
final <- data.frame(Response=y,tab3,ICC=finalICC,CV=finalCV,stringsAsFactors=F)
names(final) <- c("Response_Variable","N_subjects","N_Observations","Average_Measurements","ICC","CV")
return(final)
}
