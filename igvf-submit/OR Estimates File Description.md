# OR Estimates File Description

This file describes the output format for our main output file, which contains statistical estimates of the ratio of odds of condition occurrence given that a person carries at least onne of the variants in a given class. See the methods section for more details.

The groups of variants come from a variaty of sources, including computational variant effect prediction, functional assays, and curated datasets. Each row in this table corresponds to a different variant group, and some of these variant groups overlap.

## Column definition

- **Dataset/Predictor** (*string*) The experimental assay dataset or the name of the variant effect predictor used to define the variant group
- **Gene** (*string*) Official HGNC symbol for the gene analyzed
- **Grouping/Calibration** (*string*) The method used to define the variant group; either functional class for functional categories or calibration method for calibrated evidence point ranges
- **Consequence** (*string*) The molecular consequence to which the variant group has been filtered
- **Group/Range** (*string*) The functional class or calibrated evidence point range defining the variant group
- **Odds Ratio** (*float*) The point estimate of the odds ratio computed for association between the group of variants and the gene-associated phenotype; missing if there are 20 or fewer carrier cases or carrier controls
- **OR_LI** (*float*) The lower bound of the (unadjusted) 95% confidence interval of the odds ratio estimate; missing if there are 20 or fewer carrier cases or carrier controls
- **OR_UI** (*float*) The upper bound of the (unadjusted) 95% confidence interval of the odds ratio estimate; missing if there are 20 or fewer carrier cases or carrier controls
- **LogOR** (*float*) The natural log of the point estimate of the odds ratio; missing if there are 20 or fewer carrier cases or carrier controls
- **LogOR_LI** (*float*) The lower bound of the (unadjusted) 95% confidence interval of the natural log of the odds ratio estimate; missing if there are 20 or fewer carrier cases or carrier controls
- **LogOR_UI** (*float*) The upper bound of the (unadjusted) 95% confidence interval of the natural log of the odds ratio estimate; missing if there are 20 or fewer carrier cases or carrier controls
- **Wald p-value** (*float*) The Wald test p-value for the OR estimate.
- **p-value (Holm adjusted)** (*float*) The p-value adjusted for multiple comparisons using the Holm-Bonferroni method
- **Cohort size** (*integer*) The total size of the cohort (cases and controls meeting phenotyping criteria)
- **Total cases** (*integer*) The number cases in the cohort
- **Total controls** (*integer*) The number of controls in the cohort
- **Carrier cases** (*string*) The number of cases with at least one variant from the variant group; displayed as "≤ 20" if the count is greater than 0 but no greater than 20
- **Carrier controls** (*string*) The number of controls with at least one variant from the variant group; displayed as "≤ 20" if the count is greater than 0 but no greater than 20
- **Variants in group** (*integer*) The number of distinct nucleotide variants in the group
- **Case only variant count** (*integer*) The number of distinct nucleotide variants from the group that appear only in cases
- **Control only variant count** (*integer*) The number of distinct nucleotide variants from the group that appear only in controls
- **Overlap variant count** (*integer*) The number of distinct nucleotide variants from the group that appear both in cases and in controls
- **Variants in cohort** (*integer*) The total number of distinct nucleotide variants from the group that appear in the cohort
- **ENSG** (*string*) The Ensembl Gene ID
- **Case inclusion phenotypes** (*list of strings*) The phenotypes the presence of which include a participant in the case cohort
- **Control exclusion phenotypes** (*list of strings*) The phenotypes the presence of which exclude a participant from the control cohort

## Methods

### Definition of carrier status for a variant class

We consider variant classes defined by (i) author-provided classifications and (ii) calibrated classifications of assays and of (iii) variant effect predictors. For each variant class, we consider the set of variants that are members of the class and use the All of Us short-read whole genome sequencing matrix to identify participants with variants in the class.

### Odds ratio estimation

We fit a logistic regression model to estimate the odds ratio of disease phenotype given carrier status. Specifically, we fit a model where the response variable is case/control status, and explanatory variables are carrier status for the variant class, sex assigned at birth, age at the time of AoU data release, and 15 genetic ancestry PCA components computed by the AoU research program. This yields a point estimate and standard deviation for the log odds ratio associated with the carrier status variable, and we report 95% confidence intervals constructed using the normal distribution approximation of the distribution of the log odds based on these values.


