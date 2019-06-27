# getID

# 1.  breast_metabolomics
# 2.  Parkinsons_C8
# 3.  Parkinsons_HILICNEG
# 4.  Parkinsons_HILICPOS
# 5.  prostate_metabolomics
# 6.  Validation_Blood
# 7.  Validation_Urine

getID <- function(dat,filepath=NULL){

# Pull ID file
metabolites <- getMetabolites(dat)
ID <- metabolites$ID
filename <- paste0("ID_",dat)

# Options:

     # Write the ID file to the global environment - user can continue in R
     if (is.null(filepath)){
             Metabolomics_IDs <<- ID
             message("IDs exported into the global environment.")
     }

if (!is.null(filepath)){
             foreign::write.foreign(ID,
                                    datafile=file.path(filepath,paste0(filename,".txt")),
                                    codefile=file.path(filepath,(paste0("Import ",filename,".SAS"))),
                                    package="SAS",
                                    dataname=filename)
             message("IDs saved using the foreign package.  See documentation for details how to import into SAS")
        }
}
