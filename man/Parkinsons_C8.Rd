\name{Parkinsons_C8}
\alias{Parkinsons_C8}
\docType{data}
\title{
CPS2 Parkinsons Disease Metabolomics - C8 platform
}
\description{
These data were collected as part of a collaboration with Harvard on Parkinsons Disease.  These data are unlike the other metabolomics data we've collected using Metabolon and are not organized around compound IDs and without the extensive metadata and QC analyses.

Users should check the other Parkinsons data we have collected from HILIC_Pos and HILIC_Neg metabolomics platforms.  There are repeats in these data, some metabolites have been collected more than once across the three platforms.


}
\usage{getMetabolites("Parkinsons_C8")}
\format{
  A list with 3 data frames
  \describe{
    \item{\code{metabolites}}{A data frame of 635 observations and 14,624 variables.  QC samples have been removed.  No scaling has been done to the metabolite data and missing values have NOT been imputed}
    \item{\code{ID}}{A data frame of 635 observations and 3 variables.  The original data had several IDs, this file will help you merge if you need to go back to the original.}
    \item{\code{biochem}}{A data frame of 14,623 observations and 12 variables.  Each metabolite has a set of metadata from the lab, including a biochemical name, pathway information, and various platform information.  This dataset also includes derived metadata including the number/proportion of samples missing a given metabolite and CV quality scores.}
  }
}
\details{
The Parkinsons_C8 file included in the package should include everything an analyst requires to complete a metabolomics project with these data.  A link to the original Rdata file is listed in the references.  If you need more information on these data, please see the documentation file at that file path.
}
\references{
S:\CPS\BIOSPECIMENS\Metabolomics\Parkinsons disease\
}
\examples{
Parkinsons_C8 <- getMetabolites("Parkinsons_C8")
metabolites <- Parkinsons_C8$metabolites
ID <- Parkinsons_C8$ID
biochem <- Parkinsons_C8$biochem
}
