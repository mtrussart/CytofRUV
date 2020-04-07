library(readxl)
library(flowCore)

## Define parameters to load the data
wd_data="/Users/trussart.m/WEHI/CytofRUV/extdata/"
metadata_filename="Metadata_ex.xlsx"

## Load the metadata file
md <- readxl::read_excel(paste(wd_data, metadata_filename, sep = ""))

## Load the fcs file
setwd(wd_data)
fcs_raw <- flowCore::read.flowSet(paste(md$file_name[2],sep=""), transformation = FALSE,
                                  truncate_max_range = FALSE)

# Save fcs file from CytofRUV paper
A2=fcs_raw[[1]]
usethis::use_data(A2, overwrite = TRUE, compress = 'xz')
