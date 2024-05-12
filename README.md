# Repository for the paper titled: "Stability, resilience, and eco-evolutionary feedbacks of mutualistic networks to rising temperature"
Authors: Gaurav Baruah, Tim Lakaemper

#The codes and scripts are prepared by Gaurav Baruah and Tim Lakaemper.

#Data and Rscripts:


`figure 1.R` R script can be used to reproduce figure 1 in the main-text.

The folder `datasets_1` contains data for all plant-pollinator networks used in the paper and maintext folder contains all the R scripts to produce the figures in the maintext.

 `all_data_species.RData` and `Network_all_data.RData` are data files that could be used to reproduce figures 2 to figures 4, and supplementary figures and can be done using the R script `analysis_and_plots.R`, `figure_4_jacobian_plots.R` script.

`01_functions_cluster.R` R script contains all the functions used for dynamic eco-evolutionary simulations used for figure 1 and for the bulk simulations used in the cluster. `cluster_run_warming.R` is the R script used to generate the network data using the Theobiota clusters in Bielefeld.
