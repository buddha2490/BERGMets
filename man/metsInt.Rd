\name{metsInt}
\alias{metsInt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stratified and interaction metabolite association models
}
\description{
This function will run stratified or interaction models for metabolites with a continuous or binary outcome.  The choice of stratified vs simple interaction depends on the type of interaction variable.  See details.
}
\usage{
metsInt(dat,biochem,outcome,metabolites,intvar,
                covariates=NULL,normalize=F)
}

\arguments{
  \item{dat}{
data frame containing your metabolites and covariates
}
  \item{biochem}{
The biochemical metadata file.  Typically this is in the $biochem element of the metabolites list.  You can subset this data frame any way you wish, some variables are probably unneeded.  Recommend that you subset to only include COMP_ID, BIOCHEMICAL, SUPER_PATHWAY, SUB_PATHWAY.
}
  \item{outcome}{
Outcome for the analysis.  Logistic regression must be coded [0,1]; linear regression must be continuous numeric.  Function will identify which is which and will provide an error message if you code it wrong.
}
  \item{metabolites}{
Vector of COMP_IDs used for the models
}
  \item{intvar)}{
Name of interaction variable.  If class(intvar)=="numeric", simple interaction models will be fit coded as metabolite*intvar.  If class(intvar)=="factor", the function will return estimates for the metabolites stratified by the intvar.
  }
  \item{covariates}{
Optional.  Vector of covariate names from the dat data frame.
}
  \item{normalize}{
Logical. Set to TRUE if you haven't already mean-centered and glog transformed your metabolite data.
}
}
\details{
In order to use metsInt(), users need to take care to code their variables appropriately.  The outcome variable determines the type of model fit: a continuous outcome variable will fit a linear model; a binary outcome will fit a logistic model.  If the variables are coded incorrectly, an error/warning will be returned.

Users also must code the intvar correctly.  If intvar is a factor variable, the function will return metabolite associations stratified across each level of of intvar.  If intvar is coded as a numeric variable, the function will return associations for the interaction term.  Both models will return a p-value for the interaction based on a likelihood ratio test comparing the interaction model with a reduced model.
}
\value{
A data frame containing the results.  Results are merged with the biochem-metadata file and estimates/pvalues are nicely formatted.

Stratified models will return estimates across each level of the interaction variables.

Interaction models will return only the estimate and p-values associated with the interaction term.

Both will return an overall p-value for the interaction.

}

\author{
Brian Carter
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

# Data
breast_metabolomics <- getMetabolites("breast_metabolomics")
data(survey)
df <- left_join(survey,breast_metabolomics$metabolites,"ID")
covars <- c("AGE_INT","LASTATE")
comp.id <- names(breast_metabolomics$metabolites)[-1]
biochem <- breast_metabolomics$biochem[,c("COMP_ID","BIOCHEMICAL")]

#### stratified models

     # categorize my strata variables (quartiles):
     df$agecat <- gtools::quantcut(df$AGE_INT,4)
     metsInt(dat=df,
                     biochem=biochem,
                     outcome="BMI",
                     metabolites=comp.id,
                     intvar="agecat",
                     covariates=covars,
                     normalize=T)

#### Interaction models
     metsInt(dat=df,
               biochem=biochem,
               outcome="BMI",
               metabolites=comp.id,
               intvar="AGE_INT",
               covariates=covars,
               normalize=T)

}
