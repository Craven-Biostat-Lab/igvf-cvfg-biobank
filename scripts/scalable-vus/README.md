# A scalable approach to resolving variants of uncertain significance

This folder contains the scripts for reproducing the biobank validation presented in [**A scalable approach to resolving variants of uncertain significance**](https://doi.org/10.64898/2026.02.14.705848).

The analysis was originally run in an *All of Us* Controlled Tier workspace (CDR version 8). [Link to workspace (not publically accessible)](https://workbench.verily.com/workspaces/biobank-validation-of-variant-group-pathogenicity).

The analysis requires two external input tables:
1. `integrated_variant_effect_dataset_analysis_excalibr_assays_oddspath_models.csv.gz`: A curated dataframe of assay measurements, calibrations, and calibrated combined scores as produced by the [data processing pipeline](https://github.com/bbi-lab/igvf-scalable-vus-resolution), this is termed as one of the [variant reclassification dataset checkpoints](https://github.com/bbi-lab/igvf-scalable-vus-resolution/blob/main/docs/build_variant_reclassification_dataset.md#build-variant-reclassification-dataset).
2. `IGVFFI1328FVZD.csv.gz`: A dataframe of variant effect predictor calibrations, available from the [IGVF portal](https://data.igvf.org/tabular-files/IGVFFI1328FVZD/).

The main deliverable of this analysis is the odds ratio estimate dataframe. The resultant dataframe is available for download from the [IGVF portal](https://data.igvf.org/tabular-files/IGVFFI5542XRSO/).

Code to generate the figures that report odds ratios can be found in the [main manuscript repository](https://github.com/bbi-lab/igvf-scalable-vus-resolution), see "[Figure preparation](https://github.com/bbi-lab/igvf-scalable-vus-resolution/tree/main#3-figure-preparation)".

## Instructions for reproducing the odds ratio estimate dataframe

### Setting up your workbench environment

These instructions are specific to the Verily workbench, as of this writing, that is the only computational environment available for performing analyses of *All of Us* controlled tier data.

1. Create a Verily workbench workspace
2. Add the *All of Us* Controlled Tier v8 curated data repository to the workbench resources.
3. Add this repository to the workbench repositories with the name `igvf-cvfg-biobank` (this will automatically clone this repository to `~/repos/igvf-cvfg-biobank` in the computational environment; this path is used in scripts)
4. Create a write-accessible cloud storage bucket for saving analysis results. We use the resource ID `scalable-vus-validation` (this automatically mounts the bucket to `~/workspace/scalable-vus-validation` in the computational envoronment; this path is used in scripts)
5. Place the above-mentioned `integrated_variant_effect_dataset_analysis_excalibr_assays_oddspath_models.csv.gz` into a new folder `~/workspace/scalable-vus-validation/variant_tables/`
6. (Optionally) place `IGVFFI1328FVZD.csv.gz` into `~/workspace/scalable-vus-validation/variant_tables/`. The notebook `Annotate VEPs` contains code to attempt to download it if it is not present.

### Generation of the odds ratio estimate dataframe

The process to generate the odds ratio estimate dataframe requires running a series of notebooks in sequence.
We note which notebooks make use of [Hail](https://hail.is/) and should be run in a `JupyterLab Spark cluster for AoU` app environment with Hail installed (hereafter referred to as the "dataproc" environment). Unless indicated otherwise, a notebooks may be run in a `AoU Jupyter` app environment (hereafter referred to as the "standard" environment) to control costs or the dataproc environment. Notebooks should be run in the following order:

1. `Setup environment` Run once in each newly created environment. Installs necessary python packages in the standard environment. 
2. `Case pools with time to diagnosis` Collects phenotype-based case pools into tabular files for all conditions except those related to long QT and LDLR.
3. `Prepare control pools for loading` Collects phenotype-based control pools into tabular files for all conditions except those related to long QT and LDLR.
4. `Make long QT case and control pools` Collects case and control pools for long QT and and abnormal QT interval.
5. `Assemble gene cohorts` Asembles case-control cohorts for each of the genes we analyze (except for LDLR) from the phenotype-level case and control pools.
6. `LDLR get phenotypes` Produces a table of LDLR-related phenotypes used in [Tabet et al. (2026)](https://doi.org/10.1126/science.ady7186); this is an R notebook adapted from code provided by the authors of that paper.
7. `LDLR cohort` Creates a case-control assignment from the phenotypes defined by the previous notebook.
8. `Load common tables` Saves a table of participant demographic information and a table of the ancestry PCA components.
9. `Carrier definition Functional` **Run in dataproc** Creates carrier sets for variants grouped by functional class.
10. `OR calculation Functional` Computes ORs for variants grouped by functional class.
11. `Detailed variant statistics functional` **Run in dataproc** Counts variants in carrier cases/controls for variant groups defined by functional class.
12. `OR analysis ExCALIBR` **Run in dataproc** Performs OR analysis for variant groups defined by ExCALIBR-calibrated evidence point ranges.
13. `Detailed variant statistics ExCALIBR` **Run in dataproc** Counts variants in carrier cases/controls for variant groups defined by ExCALIBR-calibrated evidence point ranges.
14. `Per-gene VAT import` **Run in dataproc** Loads a subset of the AoU Variant Annotation Table (VAT) into files split up per gene.
15. `Annotate VEPs` Attaches molecular consequence from the VAT to the calibrated predictor table
16. `OR analysis VEP` **Run in dataproc** Performs OR analysis for variant groups defined by calibrated predictor evidence point ranges.
17. `Detailed variant statistics VEP` **Run in dataproc** Counts variants in carrier cases/controls for variant groups defined by calibrated predictor evidence point ranges.
18. `Filter for combined points` Filters the integrated variant effect table down to one-row-per-nucleotide-change.
19. `OR analysis total-points` **Run in dataproc** Performs OR analysis for variant groups defined by combined evidence points from ExCALIBR-calibrated assays and calibrated REVEL evidence point ranges.
20. `Detailed variant statistics total-points` **Run in dataproc** Counts variants in carrier cases/controls for variant groups defined by calibrated predictor evidence point ranges.
21. `OR postprocessing` Combines OR estimates and variant statistics into one table, performs multiple-comparison adjustment, generates censored table for export.