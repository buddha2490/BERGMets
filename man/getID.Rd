\name{getID}
\alias{getID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pulls IDs from the metabolite files and write to a SAS dataset for later use.
}
\description{
Each metabolite has an ID file.  The getID() function will pull these IDs and write them to a SAS dataset that can be merged with CPS survey files.
}
\usage{
getID(dat, filepath)
}


\arguments{
  \item{dat}{Character vector of the dataset you are requesting (see details below).  Options are NOT case sensitive.
}
  \item{filepath}{File path to save the ID file.  Remember that R requires forward slashes "/" in file paths.  Example:  "S:/user/bcarter/" NOT "S:\user\bcarter"
}
}
\details{

Users can pull IDs from any of 7 metabolomics projects using the following character strings. For more information, you can find documentation by searching the BERGMets documentation:

 1.  "breast_metabolomics"

 2.  "Parkinsons_C8"

 3.  "Parkinsons_HILICNEG"

 4.  "Parkinsons_HILICPOS"

 5.  "prostate_metabolomics"

 6.  "Validation_Blood"

 7.  "Validation_Urine"
}
\value{
getID() will not return any objects into the global environment; rather, a sas7bdat file will be saved in the file path with the ID file for the metabolomics project of interest.
}
\author{
Brian Carter
}
\seealso{

\examples{
getID("breast_metabolomics", "s:/user/bcarter")
}

