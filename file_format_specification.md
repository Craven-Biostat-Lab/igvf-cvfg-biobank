# File format specification for IGVF

## Log odds ratio files

Columns:

1. Dataset
    - string
    - Unique identifier for the assay/predictor/database which provides the classification
2. Gene
    - string
    - HGNC Symbol of the gene studied
3. Classifier
    - string
    - Refers to how the classification was made. Possible values include "Curated", "Calibrated", "Author Label."
4. Classification
    - string
    - Describes the classification of the class of variants. E.g. "Likely Pathogenic" for ClinVar, or "Pathogenic Moderate" for a calibrated classifier.
5. All of Us CDR Version
    - string
6. LogOR
    - float
    - Log odds ratio estimate
7. LogOR_LI
    - float
    - Lower boundary of the 95% confidence interval of the log odds ratio
8. LogOR_UI
    - float
    - Upper boundary of the 95% confidence interval of the log odds ratio

## Prevalence curve files

Columns:

1. Dataset
    - string
    - Unique identifier for the assay/predictor/database which provides the classification
2. Gene
    - string
    - HGNC Symbol of the gene studied
3. Classifier
    - string
    - Refers to how the classification was made. Possible values include "Curated", "Calibrated", "Author Label."
4. Classification
    - string
    - Describes the classification of the class of variants. E.g. "Likely Pathogenic" for ClinVar, or "Pathogenic Moderate" for a calibrated classifier.
5. All of Us CDR Version
    - string
6. Age
    - float
    - Age value corresponding to prevalence estimate
7. Survival to onset
    - float
    - Prevalence estimate at age point
8. Survival_UI
    - float
    - Upper boundary of the 95% confidence envelope of the prevalence estimate
9. Survival_LI
    - float
    - Lower boundary of the 95% confidence envelope of the prevalence estimate
