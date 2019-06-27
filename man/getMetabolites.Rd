\name{getMetabolites}
\alias{getMetabolites}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Loads and/or subsets metabolite data.
}
\description{
Each metabolomics dataset has multiple versions based on how the lab chose to scale them.  This function will pull the correct data.
}
\usage{
getMetabolites(dat, IDlist=NULL, compid=NULL)
}


\arguments{
  \item{dat}{Character vector of the dataset you are requesting (see details below).  Options are NOT case sensitive.
}
  \item{IDlist}{Optional:  Character vector of a subset of IDs.  If you have already prepared your covariate dataset, you can include a vector of IDs for subsetting the metabolite data.  If you leave it NULL, all data will be returned.
}
  \item{compid}{Optional:  Character vector of a subset of COMP_IDs.  If you know what metabolites you are looking for, you can enter a vector of these COMP_IDs and only pull these data.
}
}
\details{

Users can pull data from any of 7 metabolomics projects using the following character strings. For more information, you can find documentation by searching the BERGMets documentation:

 1.  "breast_metabolomics"

 2.  "Parkinsons_C8"

 3.  "Parkinsons_HILICNEG"

 4.  "Parkinsons_HILICPOS"

 5.  "prostate_metabolomics"

 6.  "Validation_Blood"

 7.  "Validation_Urine"
}
\value{
getMetabolites will return a list of 3 objects:

     $metabolites - Data frame of the metabolite data organized by ID

     $ID - and ID conversion file, useful if you need to pull lifelink data or refer to the original data files

     $biochem - Biochemical metadata:  includes biochemical names, pathway information, platform/batch information, QC metrics (CV, ICC, N-missing, etc)

}
\author{
Brian Carter
}

\examples{

# Access documentation for the breast_metabolomics file:

?breast_metabolomics

# Pull all data

breast_metabolites <- getMetabolites("breast_metabolomics")

names(breast_metabolites)

# Pull just some IDs

ID <- c("42101011001022",
        "55010005007081",
        "43017011009062",
        "56017043004052",
        "06026241020062")

breast_metabolites <- getMetabolites("breast_metabolomics", IDlist=ID)

nrow(breast_metabolites$metabolites)

# Pull just some metabolites

mets <- c("X34404","X21132",
          "X57463","X37183",
          "X39346")

breast_metabolites <- getMetabolites("breast_metabolomics", compid=mets)

ncol(breast$metabolites$metabolites)

# Pull just some metabolites for just some IDs

breast_metabolites <- getMetabolites("breast_metabolomics",
                                     IDlist=ID, =compid=mets)

dim(breast_metabolites$metabolites)
}

