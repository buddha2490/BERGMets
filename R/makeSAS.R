# getID

# 1.  breast_metabolomics
# 2.  Parkinsons_C8
# 3.  Parkinsons_HILICNEG
# 4.  Parkinsons_HILICPOS
# 5.  prostate_metabolomics
# 6.  Validation_Blood
# 7.  Validation_Urine

makeSAS <- function(dat,IDlist=NULL,compid=NULL,filepath){


metabolites <- getMetabolites(dat,IDlist,compid)


drop <- c("CAS","PUBCHEM","CHEMSPIDER","KEGG","HMDB_ID")
metdata <- metabolites$metabolites
biochem <- metabolites$biochem[,!names(metabolites$biochem) %in% drop]

metfile <- paste0(dat,"_",Sys.Date())
biochemfile <- paste0(dat,"_metadata_",Sys.Date())

     foreign::write.foreign(metdata,
                            datafile=file.path(filepath,paste0(metfile,".txt")),
                            codefile=file.path(filepath, paste0("Import ",metfile,".SAS")),
                            package="SAS",
                            dataname="metabolite_data")
     foreign::write.foreign(biochem,
                            datafile=file.path(filepath,paste0(biochemfile,".txt")),
                            codefile=file.path(filepath, paste0("Import ",biochemfile,".SAS")),
                            package="SAS",
                            dataname="metabolite_metadata")
     message("Data saved using the foreign package.  See documentation for details how to import into SAS.")

}
