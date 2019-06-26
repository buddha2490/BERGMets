\name{survey}
\alias{survey}
\docType{data}
\title{
Sample covariate file used in metabolomics analyses
}
\description{
This is real lifelink survey data, subset to only include a few variables.
}
\usage{data("survey")}
\format{
  A data frame with 1534 observations on the following 5 variables.
  \describe{
    \item{\code{ID}}{Subject identifier}
    \item{\code{BMI}}{BMI at blood draw - continuous numeric}
    \item{\code{AGE_INT}}{Age at blood draw - continuous numeric}
    \item{\code{LASTATE}}{a factor for time since last meal}
    \item{\code{BMIBIN}}{Binary BMI outcome (0= <25; 1=>30)}
  }
}

\examples{
data(survey)
str(survey)
}


