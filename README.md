# igvf-cvfg-biobank
This repository contains scripts and tools for performing biobank validation of variant classifications in the All of Us Workbench.

## Analysis scripts
`scripts` contains the scripts for reproducing the biobank validation presented in [**A scalable approach to resolving variants of uncertain significance**](https://doi.org/10.64898/2026.02.14.705848).

## The `cvfgaou` python package
`cvfgaou` contains the python package. It can be installed by cloning this repository and installing with `pip`:

```
git clone https://github.com/Craven-Biostat-Lab/igvf-cvfg-biobank.git
pip install -e igvf-cvfg-biobank
```

## IGVF Portal submission scripts
`igvf-submit` contains the conda environment definition and python script for submitting processed data and associated metadata.