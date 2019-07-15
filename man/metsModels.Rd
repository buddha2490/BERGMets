\name{metsModels}
\alias{metsModels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wrapper function that will run association, interaction, or stratified models
}
\description{
metsModels should do everything you need it to do.  If you include a intvar= variable, the function will run interaction or stratified models (depending on intvar class, see details).  If the intvar= argument is left NULL, simple association models will be returned.
}
\usage{
metsModels(dat, biochem, outcome, compid,intvar=NULL,
           covariates=NULL, normalize=T)
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
  \item{intvar)}{
Optional Name of interaction variable.  If left NULL, a simple association analysis will be returned.  If class(intvar)=="numeric", simple interaction models will be fit coded as metabolite*intvar.  If class(intvar)=="factor", the function will return estimates for the metabolites stratified by the intvar.
  }
  \item{covariates}{
Optional.  Vector of covariate names from the dat data frame.
}
  \item{normalize}{
Logical. Set to TRUE if you haven't already mean-centered and glog transformed your metabolite data.  Default is TRUE
}
}
\details{
the metsModels() function is a wrapper function that runs either metsAssocation() for simple association models, or metsInteraction() for interaction/stratified models.

A simple association model includes either linear or logistic regression (depending on the outcome), adjusted for covariates.  If the user includes an intvar= variable, the interaction/stratified models will be returned, depending on the intvar class.  See documentation files for metsAssociation() and metsInteraction() for details of either sets of models
}
\value{
A data frame containing the results.  Results are merged with the biochem-metadata file and estimates/pvalues are nicely formatted.

Simple association models will return estimates, p-values and metabolite metadata

Stratified models will return estimates across each level of the interaction variables.

Interaction models will return only the estimate and p-values associated with the interaction term.

Both will return an overall p-value for the interaction.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Brian Carter
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{metsAssoc}}, \code{\link{metsInt}}, \code{\link{normalizeMets}}
}
\examples{
breast_metabolomics <- getMetabolites("breast_metabolomics")
data(survey)
df <- left_join(survey,breast_metabolomics$metabolites,"ID")
covars<- c("AGE_INT","LASTATE")
comp.id <- names(breast_metabolomics$metabolites)[-1]
biochem <- breast_metabolomics$biochem[,c("COMP_ID","BIOCHEMICAL")]

# Association analysis (logistic regression)
out <- metsAssoc(dat=df,
                       biochem=biochem,
                       outcome="BMIBIN",
                       compid=comp.id[1:10],
                       covariates=covars,
                       normalize=T)  # normalizing data



# Stratified analysis (linear regression)
     metsInt(dat=df,
               biochem=biochem,
               outcome="BMI",
               compid=comp.id,
               intvar="AGE_INT",
               covariates=covars,
               normalize=T)

{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
