\name{breast_metabolomics}
\alias{breast_metabolomics}
\docType{data}
\title{
Breast cancer metabolomics file
}
\description{
These data are a subset of our original breast cancer metabolomics NCIA-03-16-16ML+.  It includes the cleaned analytic metabolite file, a data set of metadata describing each metabolite, and an ID conversion file.
}
\usage{getMetabolites("breast_metabolomics")}
\format{
  A list with 3 data frames
  \describe{
    \item{\code{metabolites}}{A data frame of 1,547 observations and 1,385 variables.  These metabolites were scaled to the daily median by Metabolon and the missings were imputed to the minimum.  QC duplicates were averaged to provide only a single sample.  No exclusions have been made:  there are 1,547 unique samples in the data, and 1,384 metabolites (plus an ID variable)}
    \item{\code{ID}}{A data frame of 1,564 observations and 5 variables.  This is an ID conversion file for the breast cancer metabolomics and includes duplicate IDs indicating the QC samples.  The original data from Metabolon had several IDs, this file will help you merge if you need to go back to the original.  Also includes the case-control and Matching variables.}
    \item{\code{biochem}}{A data frame of 1,385 observations and 18 variables.  Each metabolite has a set of metadata from Metabolon, including a biochemical name, pathway information, and various platform information.  This dataset also includes derived metadata including the number/proportion of samples missing a given metabolite and ICC_TECH / CV quality scores.}
  }
}
\details{
The breast_metabolomics file included in the package should include everything an analyst requires to complete a metabolomics project with these data.  A link to the original Rdata file is listed in the references.  If you need more information on the breast cancer metabolomics data, please see the documentation file at that file path.
}
\references{
S:\CPS\BIOSPECIMENS\Metabolomics\Breast cancer - NCIA-03-16ML\
}

\examples{
breast_metabolomics <- getMetabolites("breast_metabolomics")
metabolites <- breast_metabolomics$metabolites
ID <- breast_metabolomics$ID
biochem <- breast_metabolomics$biochem

}

