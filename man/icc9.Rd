\name{icc9}
\alias{icc9}
\title{R Version of the ICC9 SAS Macro}
\description{
The icc9() function adapts the ICC9 macro published by Ellen Hertzmark and Donna Spiegelman (2010).  It computes intraclass correlation coefficients (ICC) and coefficients of variation (CV) and their confidence intervals.  These metrics can be calculated with and without adjustment for fixed effects.
}
\usage{
icc9(y,subject,fixed=NULL,dat,method="ML")
}
\arguments{
  \item{y}{
Name of the continuous response vector (REQUIRED)
}
  \item{subject}{
Name of the subject vector (REQUIRED)
}
  \item{fixed}{
Vector of fixed effects one or more variable names (OPTIONAL)
}
  \item{dat}{
Data frame containing variables used in analysis (REQUIRED)
}
  \item{method}{
  Users can indicate "ML" for maximum likelihood method or "REML" for restricted maximum likelihood.

By default, the lmer() function calculates estimates based on restricted maximum likelihood (REML).  The original SAS ICC9 macro defaults to a maximum likelihood (ML) method, and so that is the default for this R function.
}
}
\details{
Users familiar with the icc9 SAS macro will notice a few differences with the R function.  Notably, the macro options for subsetting the data or including BY= and WHERE= options have been removed.  Users should prepare their data prior to running the function.

This function has been tested on a limmited set of data, but has so far replicated the SAS ICC and CV estimates with only an occasional rounding difference.  If the user identifies major differences in their dataset, please contact the author of this function with details.
}
\value{
Data frame containing ICC and CV estimates
\item{Response_Variable}{Reproduces the response variable for the model (the "y" argument in the function)}
\item{N_subjects}{Number of unique "subject" values in the input dataset}
\item{N_Observations}{Number of observations in the input dataset}
\item{Average_Measurements }{Average number of observations per unique subject}
\item{ICC }{Intraclass correlation coefficient with confidence intervals}
\item{CV }{Coefficient of variation with confidence intervals}
}
\references{
Users can refer to the documentation for the original SAS icc9 macro at the following web address:  https://cdn1.sph.harvard.edu/wp-content/uploads/sites/271/2012/09/icc9.pdf

}
\author{
Brian Carter (brian.carter@cancer.org)
}


\examples{

# Example 1 - basic run
icc9(y="X48761", subject="subjid", dat=qc)

# Example 2 - using BATCH1 as a fixed effect
icc9(y="X48761", subject="subjid", fixed=BATCH1, dat=qc)

# Example 3 - Using multiple fixed effects
icc9(y="X48761", subject="subjid", fixed=c(BATCH1,BATCH2), dat=qc)

# Exmample 4 - Using the function for multiple response variables
library(dplyr)
lapply(c("X48761","X19130","X53174"), function(response) {
     icc9(y=response,subject="subjid",dat=qc)
     }) \%>\% do.call("rbind",.)
}


