\name{qc}
\alias{qc}
\docType{data}
\title{
Example dataset for documenting the icc9() function
}
\description{
Sample dataset using scrambled and deidentified metabolomics data from the American Cancer Society's Cancer Prevention Study II.
}
\usage{data("qc")}
\format{
  A data frame with 138 observations on the following 6 variables.
  \describe{
    \item{\code{subjid}}{a character vector of subject IDs, triplicates of 46 unique subjects}
    \item{\code{BATCH1}}{First numeric batch variable for fixed effects adjustment}
    \item{\code{BATCH2}}{Second numeric batch variable for fixed effects adjustment}
    \item{\code{X48761}}{Normalized metabolite measurement}
    \item{\code{X19130}}{Normalized metabolite measurement}
    \item{\code{X53174}}{Normalized metabolite measurement}
  }
}
\examples{
head(qc)

}

