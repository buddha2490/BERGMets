\name{metsHeatmap}
\alias{metsHeatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a correlation heatmap using ggplot2 syntax
}
\description{
This function will create a heatmap using ggcorrplot().  The function will produce the correlation (or partial correlation) matrix using corMatrix() and then construct the figure.  Users are given a limited number of options for adjusting the figures and can supply their own theme().
}
\usage{
metsHeatmap(dat, compid, covariates = NULL, biochem = NULL, method = "pearson", hclust = T, title, subtitle = NULL, plot_colors = NULL, xaxis_size = 7, yaxis_size = 7, title_size = 12, subtitle_size = 8, plot_theme = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
data frame including metabolites and (optionally) covariates
}
  \item{compid}{
Character vector of COMP_IDs included in the heatmap
}
  \item{covariates}{
Optional vector of covariates.  If supplied, heatmap will be based on a partial correlation matrix
}
  \item{biochem}{
Data frame of biochemical metadata.  At a minimum it must have two variables: COMP_ID and BIOCHEMICAL.  The BIOCHEMICAL will be used to label the y-axis.  If left NULL, COMP_IDs will be used to label the x- and y-axes.
}
  \item{method}{
Correlation methods: "pearson" (default), "spearman", or "kendall"
}
  \item{hclust}{
Logical: use heirarchical clustering for arranging the heatmap?  Defaults to TRUE
}
  \item{title}{
Title for your plot, required.
}
  \item{subtitle}{
Optional subtitle.  If left NULL, some summary information about the correlations will be added.
}
  \item{plot_colors}{
Adjust the color scheme by provided exactly THREE colors.  Correlations range from -1.0 to 1.0.  The first color wil be -1.0, middle color will be 0, and third color will be 1.0.  The function will create gradations between these colors.
}
  \item{xaxis_size}{
Text size for the x-axis
}
  \item{yaxis_size}{
Text size for the y-axis
}
  \item{title_size}{
Text size for the title
}
  \item{subtitle_size}{
Text size for the subtitle
}
  \item{plot_theme}{
Optional ggplot2 theme.  metsHeatmap() will create a theme based on some of your options.  You can create a custom theme using theme() and include the name of this theme object in this argument to apply it.
}
}
\details{
Building heatmaps from ggcorrplot() is a tedious process, hopefully this makes it a bit easier.  Users should be able to edit most of the plot using theme() elements, or they can accept the defaults.  The function returns a ggplot2 object, so you can edit the plot after running the function.
}
\value{
ggplot object
}

\author{
Brian Carter
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ggcorrplot}},
\code{\link{corMatrix}},
\code{\link{theme}}
}
\examples{

breast_metabolomics <- getMetabolites("breast_metabolomics")
biochem <- b$biochem
metdata <- dplyr::left_join(survey,breast_metabolomics$metabolites,"ID")
comp.id <- sample(biochem$COMP_ID,50,replace=T)


png(file.path("test mets heatmap.png"),height=10,width=10,units="in",res=400)
metsHeatmap(metdata,
            compid=comp.id,
            covariates=NULL,
            biochem=NULL,
            method="pearson",
            hclust=T,
            title="Title to test my new metsHeatmap() function",
            subtitle="here's a subtitle describing something or another",
            plot_colors=NULL,
            xaxis_size=7,
            yaxis_size=7,
            title_size=12,
            subtitle_size=8,
            plot_theme=NULL)
dev.off()

}
