\name{corMatrix}
\alias{corMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates a correlation matrix adjusting for covariates
}
\description{
The corMatrix() will calculate correlation matrices between source and target variables.  If users include a vector of covariates, corMatrix() will calculate partial correlations adjusted for these covariates.  See details.
}
\usage{
corMatrix(dat, target, source, covariates = NULL, method="pearson")
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
  \item{covariates}{Optional vector of names for your adjustment variables (see details). All these variables must be numeric.
}
  \item{method}{a character string indicating which correlation coefficient (or covariance) is to be computed:  "pearson" (default), "kendall", or "spearman"}
}
}
\details{
The function will produce simple pairwise correlations or adjusted partial correlations.

Full correlation matrices are computationally intensive because they must be calculated iteratively.  For example, computing pairwise correlations for 100 metabolites equals 100^2 calculations; 500 metabolites equals 500^2 calculations.

The corMatrix() function utilizes the parallel package to split these calculations across logical CPU cores and run them in parallel.  This cuts down the time it takes to calculate these correlations considerably.  For example, including 800 metabolites on a dataset of 1500 samples will take about 8 minutes of processing time.

}
\value{
Correlation matrix:  Rows=target;  Columns=source.  Diagonal==1.0
}

\author{
Brian Carter and Becky Hodge
}


\seealso{
 \code{\link{ppcor::pcor.test}}
 \code{\link{cor}}
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

# Run a simple correlation matrix - no covariates

mat <- corMatrix(dat=df, source=comp.id, target=comp.id, method="pearson")

# Run the full partial correlation matrix

mat <- corMatrix(dat=df,
                     source=comp.id,
                     target=comp.id,
                     covariates="AGE_INT")
}

