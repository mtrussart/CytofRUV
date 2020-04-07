
# CytofRUV

<!-- badges: start -->
<!-- badges: end -->

CytofRUV: Removing unwanted variation to integrate multiple CyTOF datasets.

## Installation

This R library can be installed by running the following lines:

``` r
library(devtools)
install_github('mtrussart/CytofRUV')
```

## Using RUV to removing batch effects in CyTOF data

CytofRUV is a computational algorithm which permits the integration of data across CyTOF batches. We provided a vignette Introduction_to_CytofRUV.Rmd that explain step by step how to load and normalise datasets and also how to visualise the diagnostic plots using the R-Shiny application before and after CytofRUV normalisation.

Users are required to provide: the path to the fcs files from all the samples in the study, a metadata file containing the details of each sample and their respective .fcs files, and a panel file containing the details of all proteins used in the study. We also provided in this vignette an example of a set of data that included all .fcs files, a metadata file and a panel file.
The metadata file is an excel file with the following column names: "file_name", "sample_id", "condition", "patient_id", "batch".
The panel file is an excel file with the following column names: "fcs_colname", "antigen", "marker_class".



``` r
Introduction_to_CytofRUV.Rmd
```

-** R-Shiny interface for the identification of batch effects before and after normalisation** 
To examine the batch effects found when comparing CyTOF data from samples replicated across batches, we built an R-Shiny application that exhibits any batch effects present in samples replicated across batches using four different diagnostics plots: Median Protein Expression, Protein Expression Distributions, Clustering Results and Cluster Proportions. 

## R-Shiny interface for the identification of batch effects using samples replicated across batches

``` r
launch_Shiny()
```

-** Normalisation procedure**
The normalize_data function allow the user to adjust for batch effects with parameter settings for the CytofRUV algorithm, such as the replicated samples to use, the clusters to be nornmalised and the value of k. 

## CytofRUV procedure to remove the batch effects

``` r
normalise_data()
```
