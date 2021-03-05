library(readxl)
library(writexl)

## Define parameters to load  the data
wd_data="/Users/trussart.m/WEHI/CytofRUV/extdata/"
panel_filename="Panel_ex.xlsx"

## Load the metadata file
panel<- readxl::read_excel(paste(wd_data, panel_filename, sep = ""))

# Save example metadataset from CytofRUV paper
writexl::write_xlsx(panel,"Panel_ex.xlsx")
usethis::use_data(panel, overwrite = TRUE, compress = 'xz')
