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
9. cases_with_variants
    - nullable int
    - Number of cases with variants in the specified class, null when true count is <= 20
10. cases_without_variants
    - nullable int
    - Number of cases with variants in the specified class, null when true count is <= 20
11. controls_with_variants
    - nullable int
    - Number of cases with variants in the specified class, null when true count is <= 20
12. controls_without_variants
    - nullable int
    - Number of cases with variants in the specified class, null when true count is <= 20

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
    - Survival to onset estimate at age point
8. Survival_UI
    - float
    - Upper boundary of the 95% confidence envelope of the survival estimate
9. Survival_LI
    - float
    - Lower boundary of the 95% confidence envelope of the survival estimate
