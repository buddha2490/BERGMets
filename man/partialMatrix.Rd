\name{partialMatrix}
\alias{partialMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates a correlation matrix adjusting for covariates
}
\description{
The partialMatrix() function will calculate partial correlations between x and y while adjusting for z.  A full correlation matrix of 100 metabolites will result in 100^2 individual calculations, 500 metabolites will require 500^2 calculations.  This can take a very long time.  This function will processes these calculations in parallel and speed up the task considerably.

A 500 x 500 partial correlation matrix wih 1500 observations will take ~150seconds.
}
\usage{
partialMatrix(dat, target, source, covariates)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
A data frame containing your metabolites and covariates
}
  \item{source}{
A vector of variable names for your metabolites (see details)
}
  \item{target}{
A vector of variable names for your metabolites (see details)
}
  \item{covariates}{
A vector of names for your adjustment variables (see details).  All these variables must be numeric.
}
}
\details{
`source` and `target` variables are the two objects you wish to correlate, adjusted for `covariates`.

To calculate a full correlation matrix, the `source` and `target` vectors will be identical, but the function will work with any combination of variables.  See examples below.


}
\value{
Data frame partial correlation matrix:  Rows=source;  Columns=target.
}

\author{
Brian Carter
}


\seealso{
 \code{\link{ppcor::pcor.test}}
}


\examples{

# covariate data

data(survey)

# Metabolite data

breast_metabolomics <- getMetabolites("breast_metabolomics")

# vector of metabolites

comp.id <- names(breast_metabolomics$metabolites)[-1]

# Merge metabolites with survey data

df <- dplyr::left_join(survey, breast_metabolomics$metabolites, "ID")

# Run the full correlation matrix

mat <- partialMatrix(dat=df,
                     source=comp.id,
                     target=comp.id,
                     covariates="AGE_INT")
}

