\name{prostate_metabolomics}
\alias{prostate_metabolomics}
\docType{data}
\title{
Prostate cancer metabolomics file
}
\description{
These data are a subset of our original prostate cancer metabolomics ACS0-01-16MD+.  It includes the cleaned analytic metabolite file, a data set of metadata describing each metabolite, and an ID conversion file.  The original data is a case-cohort design.
}
\usage{getMetabolites("prostate_metabolomics")}
\format{
  A list with 3 data frames
  \describe{
    \item{\code{metabolites}}{A data frame of 556 observations and 1,266 variables.  These metabolites were scaled to the daily median by Metabolon and the missings were imputed to the minimum.  QC duplicates were averaged to provide only a single sample.  No exclusions have been made:  there are 556 unique samples in the data, and 1,265 metabolites (plus an ID variable)}
    \item{\code{ID}}{A data frame of 556 observations and 5 variables.  This is an ID conversion file for the prostate cancer metabolomics.  The original data from Metabolon had several IDs, this file will help you merge if you need to go back to the original.  Also includes the case-control and a subtype variable.}
    \item{\code{biochem}}{A data frame of 1,265 observations and 17 variables.  Each metabolite has a set of metadata from Metabolon, including a biochemical name, pathway information, and various platform information.  This dataset also includes derived metadata including the number/proportion of samples missing a given metabolite and ICC_TECH / CV quality scores.}
  }
}
\details{
The prostate_metabolomics file included in the package should include everything an analyst requires to complete a metabolomics project with these data.  A link to the original Rdata file is listed in the references.  If you need more information on the prostate cancer metabolomics data, please see the documentation file at that file path.
}
\references{
S:\CPS\BIOSPECIMENS\Metabolomics\Prostate Cancer - ACS0-01-16MD+\
}

\examples{
prostate_metabolomics <- getMetabolites("prostate_metabolomics")
metabolites <- prostate_metabolomics$metabolites
ID <- prostate_metabolomics$ID
biochem <- prostate_metabolomics$biochem

}

