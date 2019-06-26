

metsInt <- function(dat,biochem,outcome,metabolites,intvar,
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
if (class(dat[[intvar]]) == "factor"){
     message(paste0(intvar,
                    " is a factor variable.  Function will return models of metabolites stratified by ",
                    intvar))
}
if (class(dat[[intvar]]) == "numeric"){
     message(paste0(intvar,
                    " is a numeric variable.  Function will return models of (metabolites * ",intvar,")"))
}
if (anyNA(dat[[intvar]])) {
     warning("Stratification variable includes missing values that will be dropped from models")
}

# Normalize the data if requested
if (normalize==T) {
     dat <- normalizeMets(dat,metabolites)
}

if (class(dat[[intvar]])=="numeric") type="interaction"
if (class(dat[[intvar]])=="factor")  type="stratified"


little.models <- function(x){

if (type=="interaction"){
  y1 <- formula(paste(outcome, "~",paste0(c(paste0(intvar,"*",x),
                                           covariates), collapse="+")))
  tmp <- c(x,paste0(intvar,":",x))
  columns <- c(x,"Interaction")
}
if (type=="stratified"){
  y1 <- formula(paste(outcome, "~",paste0(c(paste0(intvar,"+",x,":",intvar),
                                           covariates), collapse="+")))
  tmp <- paste0(intvar,levels(df[[intvar]]),":",x)
  columns <- paste0(intvar,levels(df[[intvar]]))

}
# base model for p-int
  y2 <- formula(paste(outcome, "~", paste0(c(x,intvar,covariates),collapse="+")))


# Need a vector of names to pull out the correct estimates


if (whatmodel=="Linear"){
        fit1 <- lm(y1,data=dat)
        fit2 <- lm(y2, data=dat)
        pint <- anova(fit1,fit2)[2,"Pr(>F)"]
             b <- coef(fit1)
             ci <- confint.default(fit1)
             se <- summary(fit1)$coef[,"Std. Error"]
             p <- summary(fit1)$coefficients[,"Pr(>|t|)"]
             results <- data.frame(COMP_ID=x,
                                   columns=columns,
                                   Estimate=b[tmp],
                                   LL=ci[tmp,"2.5 %"],
                                   UL=ci[tmp,"97.5 %"],
                                   P=p[tmp],
                                   P_Int=pint,
                                   stringsAsFactors=F)
             results$columns <- factor(results$columns,results$columns)
             if (type=="stratified")   class(results) <- c("MetsStrat","data.frame")
             if (type=="interaction")  class(results) <- c("MetsInt","data.frame")
             results <- summary(results)
              }
if (whatmodel=="Logistic"){
        fit1 <- glm(y1,family="binomial",data=dat)
        fit2 <- glm(y2,family="binomial",data=dat)
        pint <- anova(fit1, fit2,test="Chisq")[2,"Pr(>Chi)"]
             b <- coef(fit1)
             ci <- confint.default(fit1)
             se <- summary(fit1)$coef[,"Std. Error"]
             p <- summary(fit1)$coefficients[,"Pr(>|z|)"]

        results <- data.frame(COMP_ID=x,
                          columns=columns,
                          Estimate=exp(b[tmp]),
                          LL=exp(ci[tmp,"2.5 %"]),
                          UL=exp(ci[tmp, "97.5 %"]),
                          P=p[tmp],
                          P_Int=pint,
                          stringsAsFactors=F)
        results$columns <- factor(results$columns,results$columns)
             if (type=="stratified")   class(results) <- c("MetsStrat","data.frame")
             if (type=="interaction")  class(results) <- c("MetsInt","data.frame")
        results <- summary(results)

}
return(results)
}
foo <- lapply(metabolites,little.models)
foo <- Reduce(function(x,y) rbind(x,y), foo)
foo <- dplyr::right_join(mybiochem,foo,"COMP_ID")
rm(mybiochem, envir=globalenv())
return(foo)
}





