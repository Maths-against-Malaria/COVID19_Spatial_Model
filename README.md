# COVID19_Spatial_Model
This repository contains the implementation of the model described in the paper "Predicting the impact of COVID-19 vaccination campaigns - a flexible age-dependent spatially-stratified predictive model, accounting for multiple viral variants and vaccines". 

The script Spatial_Model_Germany.jl illustrates the implementation of the model assuming one location (here the case of Germany). The code can be adapted to any country. The model parameters must be chosen accordingly. The comments in the scripts and the model description in the paper are good guides to help adapt the code to the considered country.

The script General_Spatial_Model.jl illustrates the implementation of the model assuming multiple locations (here two locations are assumed, i.e., Schleswig-Holstein and Saxony in Germany). The code can be adapted to any number of locations or countries. The model parameters must be chosen accordingly. The comments in the scripts and the model description in the paper are good guides to help adapt the code to the considered locations or countries.
