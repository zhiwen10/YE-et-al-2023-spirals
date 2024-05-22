# YE-et-al-2023-spirals

## Description
This repository contains code for all figures in [Ye et al, 2023b](https://doi.org/10.1101/2023.12.07.570517).
The core of code is a spiral detection algorithm. 

## Code Organization 
The preprocessing workflows are described in `pipline1/2/3/4.m` files in the root folder; 
while the figure plotting functions are organized in the `figure1/2/3/4.m` files.

The preprocessing workflows (`pipline1/2/3/4.m`) sometimes take a long time to run, with results saved in data folder.
The figure ploting functions (`figure1/2/3/4.m`) use preprocessed data to generate final plots.

The main functions for spiral detection are listed in the beginning of `pipeline1_spirals.m`, 
which are `getSpiralDetection.m` and `getSpiralGrouping.m`.


## Data download
