\name{getAUC}
\alias{getAUC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates AUC for a logistic regression.
}
\description{
getAUC() takes a series of outcomes, metabolites, and covariates and calculates AUC for each model.  You can input as many outcomes and metabolites as you want.
}
\usage{
getAUC(dat, outcome, compid, covariates=NULL, normalize=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{Data frame containing your metabolites, outcome, and covariates
}
  \item{outcome}{Numeric vector coded [0,1] for your outcome.
}
  \item{compid}{Vector of COMP_IDs
}
  \item{covariates}{Vector of covariate variable names (optional)
}
  \item{normalize}{Normalize the data prior to modeling?  (default is FALSE)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
If a single outcome is provided, getAUC() will return a data frame with the following variables
  \item{COMP_ID }{Each COMP_ID listed in compid argument}
  \item{OUTCOME }{Outcome}
  \item{AUC } {AUC calculation}

If multiple outcomes are requested, getAUC() will return a list.  Each data frame of the list corresponds to each of the outcomes.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Brian Carter
}

\seealso{
 \code{\link{pROC::roc}},
 \code{\link{normalizeMets}},
}
\examples{
# load metabolite data

metabolites <- getMetabolites("breast_metabolomics")
df <- left_join(survey,metabolites$metabolites,"ID")

comp.id <- metabolites$biochem$COMP_ID[1:10]

# code BMI categories as [0,1] - obese vs normal
df$bmi_binary <- ifelse(df$BMI < 25,0, ifelse(df$BMI >= 30,1,NA))

myAUC <- getAUC(df,
                outcome="bmi_binary",
                compid = comp.id,
                covariates="AGE_INT",
                normalize=T)
}
