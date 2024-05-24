# YE-et-al-2023-spirals

## Description
This repository contains code for all figures in [Ye et al, 2023b](https://doi.org/10.1101/2023.12.07.570517).



![Spirals](https://github.com/zhiwen10/YE-et-al-2023-spirals/blob/main/images/spirals.png)

## Code Organization 
The preprocessing workflows are described in `pipline1/2/3/4.m` files in the root folder; 
while the figure plotting functions are organized in the `figure1/2/3/4.m` files.

The preprocessing workflows (`pipline1/2/3/4.m`) sometimes take a long time to run. Therefore, the results of the preprocessing steps are saved in data folder, ahead of time.
The figure ploting functions (`figure1/2/3/4.m`) use preprocessed data to generate final plots within minutes if not seconds, after specifying the root data folder.

These `pipline1/2/3/4.m` and `figure1/2/3/4.m` files also provide an overview of the code organizing structure. Key algorithms are described in details on top of each function, such as spiral detection algorithm. For extra inforamtion, please refer to the manuscript.


## Data download
The raw and preprocessed data are shared on figshare.

Part1/2: https://doi.org/10.6084/m9.figshare.25884259

Part2/2: https://doi.org/10.6084/m9.figshare.25884280

In total, 2 parts together require 39.2 GB disk space after unzipping.
After downloading and unzipping both repositories, place all 5 subfolders under a root data folder:
`axons`, `ephys`, `spirals`, `spirals_mirror`, `tables`. 

The code here will look for data inputs from these subfolders, after specifying the root data folder in the beginning of each main script.

## Data overview
We prepared a short script 'data_overview.m' in the root, to help quickly navigate the raw data in the data folder.


## Spiral detection algorithm
The main functions for spiral detection are listed in the beginning of `pipeline1_spirals.m`, 
which are `getSpiralDetection.m` and `getSpiralGrouping.m`.



![Pipeline](https://github.com/zhiwen10/YE-et-al-2023-spirals/blob/main/images/pipeline.png)

