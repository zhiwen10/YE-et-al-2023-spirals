# YE-et-al-2023-spirals

## Description
This repository contains code for all figures in [Ye et al, 2023b](https://doi.org/10.1101/2023.12.07.570517).
The core of code is a spiral detection algorithm. 


![Spirals](https://github.com/zhiwen10/YE-et-al-2023-spirals/blob/main/images/spirals.png)

## Code Organization 
The preprocessing workflows are described in `pipline1/2/3/4.m` files in the root folder; 
while the figure plotting functions are organized in the `figure1/2/3/4.m` files.

The preprocessing workflows (`pipline1/2/3/4.m`) sometimes take a long time to run. Therefore, the results of the preprocessing steps are saved in data folder, ahead of time.
The figure ploting functions (`figure1/2/3/4.m`) use preprocessed data to generate final plots within minutes, after specifying the root data folder.

## Spiral detection algorithm
The main functions for spiral detection are listed in the beginning of `pipeline1_spirals.m`, 
which are `getSpiralDetection.m` and `getSpiralGrouping.m`.



![Pipeline](https://github.com/zhiwen10/YE-et-al-2023-spirals/blob/main/images/pipeline.png)


## Data download
The raw and preprocessed data are shared on figshare.

Part1/2: https://doi.org/10.6084/m9.figshare.25884259

Part2/2: https://doi.org/10.6084/m9.figshare.25884280

In total, 2 parts together require 39.2 GB disk space after unzipping.
After downloading and unzipping both repositories, place all 6 subfolders under a root data folder:
`axons`, `ephys`, `spirals`, `spirals_mirror`, `tables`. 

The companion code below will look for data inputs from these 5 subfolders, after specifying the root data folder.
