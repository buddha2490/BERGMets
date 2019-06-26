\name{Parkinsons_HILICPOS}
\alias{Parkinsons_HILICPOS}
\docType{data}
\title{
CPS2 Parkinsons Disease Metabolomics - HILIC-Pos platform
}
\description{
These data were collected as part of a collaboration with Harvard on Parkinsons Disease.  These data are unlike the other metabolomics data we've collected using Metabolon and are not organized around compound IDs and without the extensive metadata and QC analyses.

Users should check the other Parkinsons data we have collected from C8 and HILIC_Neg metabolomics platforms.  There are repeats in these data, some metabolites have been collected more than once across the three platforms.


}
\usage{getMetabolites("Parkinsons_HILICPOS")}
\format{
  A list with 3 data frames
  \describe{
    \item{\code{metabolites}}{A data frame of 636 observations and 11,626 variables.  QC samples have been removed.  No scaling has been done to the metabolite data and missing values have NOT been imputed}
    \item{\code{ID}}{A data frame of 636 observations and 3 variables.  The original data had several IDs, this file will help you merge if you need to go back to the original.}
    \item{\code{biochem}}{A data frame of 11,625 observations and 12 variables.  Each metabolite has a set of metadata from the lab, including a biochemical name, pathway information, and various platform information.  This dataset also includes derived metadata including the number/proportion of samples missing a given metabolite and CV quality scores.}
  }
}
\details{
The Parkinsons_HILICPOS file included in the package should include everything an analyst requires to complete a metabolomics project with these data.  A link to the original Rdata file is listed in the references.  If you need more information on these data, please see the documentation file at that file path.
}
\references{
S:\CPS\BIOSPECIMENS\Metabolomics\Parkinsons disease\
}
\examples{
Parkinsons_HILICPOS <- getMetabolites("Parkinsons_HILICPOS")
metabolites <- Parkinsons_HILICPOS$metabolites
ID <- Parkinsons_HILICPOS$ID
biochem <- Parkinsons_HILICPOS$biochem
}
