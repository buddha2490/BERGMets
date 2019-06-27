\name{makeSAS}
\alias{makeSAS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a SAS datafile from metabolomics data
}
\description{
makeSAS() will import metabolomics data and convert it to a SAS file.  You can subset the metabolomics files based on IDs or COMP_IDs.  See details how to finish the import in SAS
}
\usage{
makeSAS(dat,IDlist=NULL,compid=NULL,filepath)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{Character vector of the dataset you are requesting (see details below). Options are NOT case sensitive.}
  \item{IDlist}{Optional: Character vector of a subset of IDs. If you have already prepared your covariate dataset, you can include a vector of IDs for subsetting the metabolite data. If you leave it NULL, all data will be returned.}
   \item{compid}{Optional: Character vector of a subset of COMP_IDs. If you know what metabolites you are looking for, you can enter a vector of these COMP_IDs and only pull these data.}
  \item{filepath} {File path to save the files.  Remember that R requires forward slashes "/" in file paths. Example: "S:/user/bcarter/" NOT "S:\user\bcarter"}

}
\details{
Users can pull data from any of 7 metabolomics projects using the following character strings. For more information, you can find documentation by searching the BERGMets documentation:

makeSAS() will create 4 files at the output location.  Two .txt files will be fore the actual metabolites files and a metadata file listing all the relevant data describing each metabolite including biochemical name, pathway information, platform, and some QC metrics.

The other two files are .SAS programs.  In order to get the metabolite data into SAS, simply run these SAS programs.  This will create a dataset in the WORK library that you can work with in whatever way you want.

 }
\value{
Four files will be created in your output folder.  Naming conventions based on the specific metabolomics project you are working on (see examples)

     1.  breast_metabolomics_{date}.txt
     2.  breast_metabolomics_metadata_{date}.txt
     3.  import Breast_metabolomics_{date}.SAS
     4.  import breast_metabolomics_metadata_{date}.SAS

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
\code{\link{getMetabolites}
}}
\examples{

makeSAS(dat="breast_metabolomics",
        filepath="s:/users/bcarter/")


}
