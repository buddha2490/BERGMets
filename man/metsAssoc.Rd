\name{metsAssoc}
\alias{metsAssoc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Basic metabolomics association analysis.
}
\description{
This function will run a basic association analysis adjusted for covariates for a given list of metabolites.  Depending on the outcome, the function will choose a linear or logistic regression model.  It will calculate the associations between outcome and each metabolite, merge with the metadata files, and return an organized table.
}
\usage{
metsAssoc(dat, biochem, outcome, metabolites, covariates=NULL, normalize=F)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{compid}{
Vector of COMP_IDs used for the models
}
  \item{covariates}{
Optional.  Vector of covariate names from the dat data frame.
}
  \item{normalize}{
Logical. Set to TRUE if you haven't already mean-centered and glog transformed your metabolite data.
}
}
\details{

}
\value{
Data frame organizing the associations and p-values for all metabolites.
}
\author{
Brian Carter
}

\examples{
breast_metabolomics <- getMetabolites("breast_metabolomics")
data(survey)
df <- left_join(survey,breast_metabolomics$metabolites,"ID")
covars<- c("AGE_INT","LASTATE")
comp.id <- names(breast_metabolomics$metabolites)[-1]
biochem <- breast_metabolomics$biochem[,c("COMP_ID","BIOCHEMICAL")]

# Linear regression models
out <- metsAssoc(dat=df,
                       biochem=biochem,
                       outcome="BMI",
                       compid=comp.id[1:10],
                       covariates=covars,
                       normalize=T)  # normalizing data

# Logistic regression models
out <- metsAssoc(dat=df,
                       biochem=biochem,
                       outcome="BMIBIN",
                       compid=comp.id[1:10],
                       covariates=covars,
                       normalize=T)  # normalizing data

}

