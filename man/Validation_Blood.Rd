\name{Validation_Blood}
\alias{Validation_Blood}
\docType{data}
\title{
Metabolomics in the CPS3 validation study - Blood samples
}
\description{
These data are a subset of our original metabolomics data collected for the CPS3 validation study.  It includes the cleaned analytic metabolite file, a data set of metadata describing each metabolite, and an ID conversion file.

The data were collected at two time points, so the metabolite dataset includes a "Round" variable indicating the timing of data collection.

}
\usage{getMetabolites("Validation_Blood")}
\format{
  A list with 3 data frames
  \describe{
    \item{\code{metabolites}}{A data frame of 1,472 observations and 1,370 variables.  These metabolites were scaled to the daily median by Metabolon and the missings were imputed to the minimum.  Samples were further normalized by sample volume.  QC duplicates were averaged to provide only a single sample.  No exclusions have been made, but there are duplicates.  Duplicate samples were collected over two period, and there is a "Round" indicating the timing of each collection.  In total there are 757 unique CPS3 participants in the file.}
    \item{\code{ID}}{A data frame of 1,533 observations and 3 variables.  The original data from Metabolon had several IDs, this file will help you merge if you need to go back to the original.}
    \item{\code{biochem}}{A data frame of 1,368 observations and 17 variables.  Each metabolite has a set of metadata from Metabolon, including a biochemical name, pathway information, and various platform information.  This dataset also includes derived metadata including the number/proportion of samples missing a given metabolite and ICC_TECH / CV quality scores.}
  }
}
\details{
The Validation_Blood file included in the package should include everything an analyst requires to complete a metabolomics project with these data.  A link to the original Rdata file is listed in the references.  If you need more information on these data, please see the documentation file at that file path.
}
\references{
S:\CPS3\Biospecimens\Metabolomics\Validation study blood - ACS0-01-18ML\
}
\examples{
Validation_Blood <- getMetabolites("Validation_Blood")
metabolites <- Validation_Blood$metabolites
ID <- Validation_Blood$ID
biochem <- Validation_Blood$biochem
}
