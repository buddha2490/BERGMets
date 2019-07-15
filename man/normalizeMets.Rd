\name{normalizeMets}
\alias{normalizeMets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normalize metabolite data prior to analyses
}
\description{
This function normalizes metabolite values by creating a log-transformed z-score.  This should be done prior to any analysis, although the function can also be done on the fly in the modeling functions by setting the normalize=TRUE argument.

Users have the choice of normalizing an entire dataset or selected metabolites. It is important that the data be subset to your analytic sample prior to running the function.  Either merge with your covariate data or subset the IDs.
}
\usage{
normalizeMets(dat, metabolites=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame containing metabolite data.
}
  \item{metabolites}{
Optional.  Vector of COMP_IDs of the metabolites in your data.
}
}
\details{
Users have several options for running this function, but the data must be subset to your analytic sample prior to running.  This can be done by merging the metabolites with your survey data or just subsetting the metabolite dataset to some set of IDs.

If you choose to merge the metabolites with your survey data, you MUST include a vector of COMP_IDs in the metabolites argument.  See examples below how you might want to do this.

If you choose to normalize the entire dataset, you only need to include the dataset in the `dat` argument, no vector of metabolites is required.  See examples below.
}
\value{
Data frame of metabolites mean centered and log-transformed
}

\author{
Brian Carter
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

### Example 1 - normalize an entire dataset
     data(breast_metabolomics)
     normalized_data <- normalizeMets(breast_metabolomics$metabolites)

### Example 2 - merge with your survey data first
     data(breast_metabolomics)
     mydata <- left_join(survey_data, breast_metabolomics$metbolites, "ID")

     # create a vector of COMP_IDs
     mets <- c("X101","X102","X103")

     # normalize the data
     mydata <- normalizeMets(mydata,metabolites=mets)
}

