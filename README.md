# Repository for the paper titled " Mutualistic network architecture and eco-evolutionary feedbacks modulate the occurrence of transitions and stability in response to rising temperature"
Authors: Gaurav Baruah, Tim Lakaemper

### This repository contains the following Rscripts, and data used for producing the figures

The codes and scripts are prepared by Gaurav Baruah and Tim Lakaemper.

#### Data and Rscripts:
`figure 1.R` R script can be used to reproduce figure 1 in the main-text.


The folder `datasets_1` contains data for all plant-pollinator networks used in the paper and `maintext` folder contains all the R scripts to produce the figures in the maintext.

Inside the `maintext` folder, `Species_data.RData` and `Network_data.RData` are data files that could be used to reproduce figures 2 to figures 4, and supplementary figures and can be done using the R script `warming_analysis_jacobian.R` script. 

`functions.R` R script contains all the functions used for dynamic eco-evolutionary simulations used for figure 1 and others.
`cluster_run_warming.R` is the R script used to generate the network data using the Theobiota clusters in Bielefeld.
