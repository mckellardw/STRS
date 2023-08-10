# STRS scripts

## `download_scripts`
Scripts to download data used in [STRS preprint](https://www.biorxiv.org/content/10.1101/2022.04.20.488964v1)
#TODO:
- VASAseq data
- SST data

## `rmd_notebooks`
R markdown notebooks containing analysis (and figures) from the [STRS preprint](https://www.biorxiv.org/content/10.1101/2022.04.20.488964v1)

Tools/resources used:
- Spot deconvolution was done with BayesPrism
  - [Version that we used](https://github.com/Danko-Lab/TED)
  - [New version that is much faster...](https://github.com/Danko-Lab/BayesPrism)
- scMuscle single-cell transcriptomic reference:
  - [paper](https://www.nature.com/articles/s42003-021-02810-x)
  - [data download](https://datadryad.org/stash/dataset/doi:10.5061%2Fdryad.t4b8gtj34)
  - [code](https://github.com/mckellardw/scMuscle)

## `misc_R_scripts`
- `seurat_helpers.R`: utility functions that I commonly use alongside Seurat. Check out the development version(s) [here](https://github.com/mckellardw/DWM_utils/tree/main/sc_utils/seurat_helpers)
- `McKolors_v1.csv`: color palettes I have gathered from other sources; info [here](https://github.com/mckellardw/DWM_utils/tree/main/plotting_utils)
- `scThemes.R`: custom ggplot2 themes for different types of plots; info [here](https://github.com/mckellardw/DWM_utils/tree/main/plotting_utils)

## `figures`
R scripts used to generate the figures for the [STRS preprint](https://www.biorxiv.org/content/10.1101/2022.04.20.488964v1).
