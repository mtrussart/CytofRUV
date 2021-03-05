library(readxl)
library(writexl)

## Define parameters to load the data
wd_data="/Users/trussart.m/WEHI/CytofRUV/extdata/"
metadata_filename="Metadata_ex.xlsx"

## Load the metadata file
md <- readxl::read_excel(paste(wd_data, metadata_filename, sep = ""))

# Save example metadataset from CytofRUV paper
writexl::write_xlsx(md,"Metadata_ex.xlsx")
usethis::use_data(md, overwrite = TRUE, compress = 'xz')
