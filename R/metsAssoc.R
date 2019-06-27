metsAssoc <- function(dat,biochem,outcome,compid,
                            covariates=NULL,normalize=F) {

whatmodel <- ifelse(length(unique(dat[[outcome]][!is.na(dat[[outcome]])])) < 2, NA, ifelse(
     length(unique(dat[[outcome]][!is.na(dat[[outcome]])]))==2,"Logistic","Linear"))

mybiochem <<- biochem
# Check the data
if (is.na(whatmodel)) stop("Outcome must have exactly 2 levels for logistic regression or more than 2 for linear models")

if (whatmodel=="Logistic"){
     if (!class(dat[[outcome]]) %in% c("numeric","double","integer")){
          stop("Outcome must be numeric for Logistic models")}
     if (sum(unique(dat[[outcome]]),na.rm=T) != 1){
          stop("Outcome must be coded 0=control, 1=case for Logistic models")
     }
     if (anyNA(dat[[outcome]])) {
      n <- nrow(dat[is.na(dat[[outcome]]),])
          warning(paste0("Outcome variable contains ",n," missing values that will be dropped from models"))
    }
}
if (whatmodel=="Linear"){
     if (class(dat[[outcome]]) != "numeric"){
          stop("Outcome must be a numeric variable for linear models")
     }
     if (anyNA(dat[[outcome]])) {
      n <- nrow(dat[is.na(dat[[outcome]]),])
          warning(paste0("Outcome variable contains ",n," missing values that will be dropped from models"))
     }
}
if (anyNA(dat[,covars])) {
     warning("Covariates include missing values and will be dropped from models")
}

# Normalize the data if requested
if (normalize==T) {
     dat <- normalizeMets(dat,compid)
}

# Run the models
little.models <- function(x){
     y <- formula(paste(outcome, "~",paste0(c(x, covariates),collapse="+")))
     if (whatmodel=="Linear"){
          fit <- lm(y,data=dat)
             b <- coef(fit)
             ci <- confint.default(fit)
             se <- summary(fit)$coef[,"Std. Error"]
             p <- summary(fit)$coefficients[,"Pr(>|t|)"]
          f <- data.frame(COMP_ID=x,
                          Estimate=b[x],
                          LL=ci[x,"2.5 %"],
                          UL=ci[x,"97.5 %"],
                          P=p[x],
                          stringsAsFactors=F)
          } else {
          fit <- glm(y,family="binomial",data=dat)
             b <- coef(fit)
             ci <- confint.default(fit)
             se <- summary(fit)$coef[,"Std. Error"]
             p <- summary(fit)$coefficients[,"Pr(>|z|)"]
          f <- data.frame(COMP_ID=x,
                          Estimate=exp(b[x]),
                          LL=exp(ci[x,"2.5 %"]),
                          UL=exp(ci[x, "97.5 %"]),
                          P=p[x],
                          stringsAsFactors=F)
          }
     return(f)
}

foo <- lapply(compid, little.models)
class(foo) <- "MetsResults"
foo <- summary(foo)
rm(mybiochem, envir=globalenv())
return(foo)
}


