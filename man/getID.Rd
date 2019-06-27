\name{getID}
\alias{getID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pulls IDs from the metabolite files and optionally write to a SAS dataset for later use.
}
\description{
Each metabolite has an ID file.  The getID() function will pull these IDs.  By default, getID() will create a data frame in your global environment.  You can also request export to SAS via the foreign:: package.  See details.
}
\usage{
getID(dat,filepath=NULL)
}


\arguments{
  \item{dat}{Character vector of the dataset you are requesting (see details below).  Options are NOT case sensitive.
}
  \item{filepath}{Optional file path to save the ID file.  Remember that R requires forward slashes "/" in file paths.  Example:  "S:/user/bcarter/" NOT "S:\user\bcarter".  If you leave filepath=NULL, getID() will simply create an ID file within your global environment.
}

}
\details{

Users can pull IDs from any of 7 metabolomics projects using the following character strings. For more information, you can find documentation by searching the BERGMets documentation


SAS foreign files:  The foreign package does two things.  First, it creates a .txt file of the data you are exporting from R.  Second, it writes a SAS program that will import this .txt file into SAS and format everything correctly.  Selecting this option will write a data file and a .SAS program to your filepath.  Simply run the SAS program and your data will be imported.


}
\value{
getID() will only return a data frame named "Metabolomics_IDs" if you leave filepath=NULL.  Otherwise, it will save an .xpt or .txt/.sas files to your filepath.  You can move on to SAS to finish your analysis.

}
\seealso{
\code{\link{getMetabolites}
}}
\author{
Brian Carter
}


\examples{

# Pull the IDs for use in R

# Will return a data frame named "Metabolomics_IDs"

getID("breast_metabolomics")


# Export to SAS using the foreign package

# Creates 2 files:

  # "Import ID_Breast_metabolomics.SAS"

  # "ID_breast_metabolomics"

# Simply run the "Import ID_Breast_metabolomics.SAS" and continue in SAS

  getID(dat="breast_metabolomics",
        filepath="S:/user/bcarter/metabolites")

}

