Preparation of observation data for multi-site dispersion normalized PMF analysis at available U.S. monitoring sites.
=====
Ting Zhang 2023-01

Data used in this work
-------
The repository references two datasets about PM2.5 specie observations managed by the [Cooperative Institute for Research in the Atmosphere](http://views.cira.colostate.edu/fed/QueryWizard/Default.aspx), including the Interagency Monitoring of Protected Visual Environments (IMPROVE) and the Chemical Speciation Network (CSN). 

### CSN & IMPROVE data

We selected the "IMPROVE Aerosol" and "EPA CSN" sampled at a frequency of 1-in-3 day, included all outputs of every parameter from all monitoring sites, set start date as 2011-01-01 and end date as 2020-12-31, used the "double quote" and downloaded .txt format data and metadata files separately.  

R scripts in this repository
-------
There are four R scripts in the folder now, they are:

A) `CSN Explore.R` and B) `IMPROVE Expl.R` documents the basic exploration of CSN/IMPROVE dataset, such as the spatiotemporal distribution of PM2.5 species concentrations and comparison of the flagged and unflagged PM2.5 species concentrations.

C) `CSN_IMPROVE_compare.R` documents the detection and comparison of side-by-side observations generated from CSN and IMRPOVE datasets.

D) `prepare_data_for_interpolation.R`   
In the origianl datasets, the data were organized to show the concentration of one PM species sampled on a given day on a given site for one row. With this R script, we convert it to show the concentrations of all PM species sampled on a given date at given site within one row.
