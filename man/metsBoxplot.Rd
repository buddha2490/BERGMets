\name{metsBoxplot}
\alias{metsBoxplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates layered boxplots of metabolite values across a class variable.
}
\description{
Boxplots are often useful for illustrating the distribution of metabolite values across a categorical outcome variable.  metsBoxplot() will create these boxplots using ggplot2.  The plots are somewhat customizable, and can be further customizable by including a custom theme object.
}
\usage{
metsBoxplot(dat, classvar, metabolite, ylog = F, boxcolors = NULL, ref_line = NULL, title, subtitle = NULL, xlab = NULL, ylab = NULL, xaxis_size = 12, yaxis_size = 12, titlesize = 14, plot_theme = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Data frame used for the analysis
}
  \item{classvar}{
Factor variable used for the x-axis
}
  \item{metabolite}{
Character vector (length==1) of the COMP_ID used for the y-axis.  May be raw or log transformed
}
  \item{ylog}{
  Logical, Plot y-axis on a log10 scale.
}
  \item{boxcolors}{
Vector of colors equal to the length of levels(classvar).  If you leave this NULL, a black and white image will be returned, but you can color each of the boxes if needed.
}
  \item{ref_line}{
Numeric vector indicating a horizontal reference line.  You may want a line to indicate overall minimum or median, or any other reference you might like to plot.
}
  \item{title}{
Character vector of a plot title, optional.
}
  \item{subtitle}{
Character vector of a subtitle, optional
}
  \item{xlab}{
X-axis label.  Defaults to the classvar.
}
  \item{ylab}{
Y-axis label.  Defaults to metabolite COMP_ID
}
  \item{xaxis_size}{
Font size for x-axis
}
  \item{yaxis_size}{
Font size for y-axis
}
  \item{titlesize}{
Font size for title
}
  \item{plot_theme}{
You can create your own ggplot2 theme object and include it here.  If you leave it as NULL, the plot will return with Brian's preferred aesthetic tastes, for better or worse.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
ggplot2 object.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Brian Carter
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

breast <- getMetabolites("breast_metabolomics")

df <- left_join(survey,
                breast$metabolites,
                "ID")

df$BMICAT <- with(df,
                  ifelse(BMI < 18.5, 1, ifelse(
                       BMI < 25, 2, ifelse(
                            BMI < 30, 3, ifelse(
                                 BMI < 35, 4, 5)
                            )))) %>%
     factor(.,1:5, c("<18.5", "18.5-24.9",
                     "25-29.9", "30-34.9",
                     "35+"))

myplot <- metsBoxplot(dat=df,
            classvar="BMICAT",
            metabolite="X42459",
            ylog = T,
            boxcolors = c("royalblue","darkred","yellow",
                          "orange","forestgreen"),
            ref_line = min(df$X42459),
            title = "Title goes here",
            subtitle = "Subtitle goes here",
            xlab = "BMI Categories",
            ylab = "Sphingomyelin",
            xaxis_size = 12,
            yaxis_size = 12,
            titlesize = 14,
            plot_theme = NULL)


png("Boxplot.png", height=6, width=6, units="in", res=400)
myplot
dev.off()

}

