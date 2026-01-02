# OR Estimates File Description

This file describes the output format for our main output file, which contains statistical estimates of the ratio of odds of condition occurrence given that a person carries at least onne of the variants in a given class. See the methods section for more details.

The classes of variants come from a variaty of sources, including computational variant effect prediction, functional assays, and curated datasets. Each row in this table corresponds to a different variant class, and some of these variant classes overlap.

## Column definition

- Dataset - string - The dataset on which the variant classification is based. This can be the name of a VEP (e.g. "AlphaMissense") or a functional study (e.g. "BRCA1_Findlay_2018").
- Gene Symbol - string - The gene symbol of the gene evaluated, (e.g. "BRCA1").
- ENSG - string - The Ensembl Gene ID
- Classifier - string - The method by which variant classes are determined, (e.g. "Author Reported" or "Calibrated (gene-specific)")
- Classification - string - The classification that defines variants in the evaluated class, e.g. "≥ +1"
- LogOR - float - The point estimate of the natural log of the odds ratio for the variant class (see Methods for details).
- p-value - float - The two tailed p-value under the null hypothesis that the log odds ratio is zero.
- LogOR_LI - float - The lowerbound of the two-sided 95% confidence interval of the log odds ratio estimate. 
- LogOR_UI - float - The upperbound of the two-sided 95% confidence interval of the log odds ratio estimate. 
- Cases with variants - int - The number of participants in the case cohort that carry a variant in the class, rounded up to a multiple of 20.
- Controls with variants - int - The number of participants in the control cohort that carry a variant in the class, rounded up to a multiple of 20.
- Cases without variants - int - The number of participants in the case cohort that do not carry a variant in the class, rounded up to a multiple of 20.
- Controls without variants - int - The number of participants in the control cohort that do not carry a variant in the class, rounded up to a multiple of 20.
- Variants per case - float - The mean number of variants in the class that a participant in the case cohort carries. 
- Variants per control - float - The mean number of variants in the class that a participant in the case cohort carries.
- Case variant min AF - float - The minimum allele frequency seen in a class variant in a participant in the case cohort.
- Control variant min AF - float - The minimum allele frequency seen in a class variant in a participant in the control cohort.
- Case variant max AF - float - The maximum allele frequency seen in a class variant in a participant in the case cohort.
- Control variant max AF - float - The maximum allele frequency seen in a class variant in a participant in the control cohort.
- Case only variant count - int - The number of variants from the variant class that appear only in the case cohort.
- Control only variant count - int - The number of variants from the variant class that appear only in the control cohort.
- Overlap variant count - int - The number of variants from the variant class that appear in both the case cohort and the control cohort.
- Variants in class - int - The number of variants that the class of interest contains.
- Variants in cohort - int - The number of variants from the variant class that were actually observed in the biobank cohort.
- Class ClinVar Pathogenic - int - The number of variants from the class that are annotated as Pathogenic in CinVar.
- Class ClinVar Likely pathogenic - int - The number of variants from the class that are annotated as Likely pathogenic in CinVar.
- Class ClinVar Likely benign - int - The number of variants from the class that are annotated as Likely benign in ClinVar.
- Class ClinVar Benign - int - The number of variants from the class that are annotated as Benign in CinVar.
- Class ClinVar Conflicting - int - The number of variants from the class that are annotated as having conflicting evidence in CinVar.
- Class ClinVar Uncertain significance - int - The number of variants from the class that are annotated as variants of uncertain significance in ClinVar.
- Class ClinVar Other / not in ClinVar - int - The number of variants from the class that don't have a ClinVar annotation or don't fit into one of the other ClinVar categories above.
- Cohort ClinVar Pathogenic - int - The number of variants from the the variant class that were actually observed in the biobank cohort that are annotated as Pathogenic in CinVar.
- Cohort ClinVar Likely pathogenic - int - The number of variants from the variant class that were actually observed in the biobank cohort that are annotated as Likely pathogenic in CinVar.
- Cohort ClinVar Likely benign - int - The number of variants from the variant class that were actually observed in the biobank cohort that are annotated as Likely benign in ClinVar.
- Class ClinVar Benign - int - The number of variants from the variant class that were actually observed in the biobank cohort that are annotated as Benign in CinVar.
- Class ClinVar Conflicting - int - The number of variants from the variant class that were actually observed in the biobank cohort that are annotated as having conflicting evidence in CinVar.
- Class ClinVar Uncertain significance - int - The number of variants from the variant class that were actually observed in the biobank cohort that are annotated as variants of uncertain significance in ClinVar.
- Class ClinVar Other / not in ClinVar - int - The number of variants from the variant class that were actually observed in the biobank cohort that don't have a ClinVar annotation or don't fit into one of the other ClinVar categories above.
- Case inclusion phenotypes - list of string - The conditions the presence of which include a participant in the case cohort.
- Control exclusion phenotypes - list of string - The conditions the presence of which exclude a participant from the control cohort.
- SpliceAI filter max - float - If present, participants carrying splice variants in the gene were excluded from the analysis. The value in this column is the SpliceAI score threshold used for determining this exclusion criterion.
- Data Version - string - Identifies relevant version of functional data, if applicable.
- AoU CDR - string - The version of the All of Us data release used.

## Methods

### Definition of carrier status for a variant class

We consider variant classes defined by (i) author-provided classifications and (ii) calibrated classifications of assays and of (iii) variant effect predictors. For each variant class, we consider the set of variants that are members of the class and use the All of Us short-read whole genome sequencing matrix to identify participants with variants in the class.

### Odds ratio estimation

We fit a logistic regression model to estimate the odds ratio of disease phenotype given carrier status. Specifically, we fit a model where the response variable is case/control status, and explanatory variables are carrier status for the variant class, sex assigned at birth, age at the time of AoU data release, and 15 genetic ancestry PCA components computed by the AoU research program. This yields a point estimate and standard deviation for the log odds ratio associated with the carrier status variable, and we report 95% confidence intervals constructed using the normal distribution approximation of the distribution of the log odds based on these values.


