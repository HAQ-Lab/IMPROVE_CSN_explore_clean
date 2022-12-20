# IMPROVE_CSN_explore_clean
Exploration, flag detection, NA interpolation, clustering of the IMPROVE &amp; CSN data downloaded on CIRA website.

A) CSN Explore.R/IMPROVE Expl.R:
basic exploration of CSN/IMPROVE dataset

B) CSN_IMPROVE_compare.R
1) detect the sampling sites in CSN and IMPROVE that share identical longitude and latitude;
2) compare the correlation of given PM species sampled side-by-side in the two datasets.

C) prepare_data_for_interpolation.R
In the origianl datasets, the data were organized to show the concentration of one PM species sampled on a given day on a given site for one row. With this R script, we convert it to show the concentrations of all PM species sampled on a given date at given site within one row.

