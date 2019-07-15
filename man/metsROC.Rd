\name{metsROC}
\alias{metsROC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates ROC curves using ggplot()
}
\description{
Users provide a binary outcome, metabolites, and covariates.  Logistic regression models are fit for each of the metabolites and ROC curves are returned.  Users have the choice of plotting single or multiple curves on each figure.
}
\usage{
metsROC(dat, outcome, compid, covariates, normalize=T, xlabel = "1 - Specificity",
ylabel = "Sensitivity", title = NULL, colors = NULL, title_size = 10,
subtitle_size = 10, xaxis_size = 9, yaxis_size = 9, plot_theme = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame containing outcome, metabolites, and covariates data
}
  \item{outcome}{
Binary outcome variable, must be coded [0,1] for a logistic regression.
}
  \item{compid}{
Vector of metabolite COMP_IDs.  If you include more than one COMP_ID, all will be plotted on a single figure.
}
  \item{covariates}{
Vector of covariate names.
}
  \item{normalize}{
Logical. Set to TRUE if you haven't already mean-centered and glog transformed your metabolite data.  Default is TRUE.
}
  \item{xlabel}{
X-axis title - defaults to "1 - Specificity"
}
  \item{ylabel}{
Y-axis title - defaults to "Sensitivity"
}
  \item{title}{
Optional title.  If you don't include a title, the function will provide one for you.
}
  \item{colors}{
Optional vector of colors for each curve.  length(colors) must equal length(compid)
}
  \item{title_size}{
Optional text size for title.
}
  \item{subtitle_size}{
Optional text size for subtitle.
}
  \item{xaxis_size}{
Optional text size for x-axis.
}
  \item{yaxis_size}{
Optional text size for y-axis.}
  \item{plot_theme}{
Optional ggplot theme object.  Users can create their own custom theme.  This provides the ability to edit parts of the figure as needed.
}
}
\details{
the metsROC() function uses ggplot2 to create ROC curves for any number of metabolites.  Users should be able to edit the most important aspects of the plot using the arguments in the function; others can be adjusted using theme() elements, or accept the defaults.  The function returns a ggplot2 object, so you can edit the plot after running the function.

}
\value{
a ggplot2 object
}

\author{
Brian Carter and Becky Hodge
}
\seealso{
 \code{\link{ggplot2}},
 \code{\link{plotROC}},
 \code{\link{normalizeMets}}
}
\examples{
breast_metabolomics <- getMetabolites("breast_metabolomics")

biochem <- breast_metabolomics$biochem

metdata <- dplyr::left_join(survey,breast_metabolomics$metabolites,"ID")

comp.id <- sample(biochem$COMP_ID,50,replace=T)

# categorize BMI
metdata$BMIBIN <- metdata$BMIBIN-1



png("test ROC.png",height=6,width=6,units="in",res=400)
metsROC(dat=metdata,
        outcome="BMIBIN",
        compid=comp.id[1:5],
        covariates="AGE_INT",
        normalize=T,
        xlabel="1 - Specificity",
        ylabel="Sensitivity",
        title=NULL,
        colors=NULL,
        title_size=10,
        #subtitle_size=12,
        #xaxis_size=12,
        #yaxis_size=12,
        plot_theme=NULL)
dev.off()
}
