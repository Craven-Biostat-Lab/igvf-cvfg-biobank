""" Data that is project-specific
Currently planning to move a lot of this to resource files.
"""

project_genes = [
    'ASPA',
    'BAP1',
    'BARD1',
    'BRCA1',
    'BRCA2',
    'CALM1',
    'CALM2',
    'CALM3',
    'CARD11',
    'CBS',
    'CHEK2',
    'CRX',
    'CTCF',
    'DDX3X',
    'F9',
    'FKRP',
    'G6PD',
    'GCK',
    'HMBS',
    'JAG1',
    'KCNE1',
    'KCNH2',
    'KCNQ4',
    'LARGE1',
    'MSH2',
    'NDUFAF6',
    'OTC',
    'PALB2',
    'PAX6',
    'PTEN',
    'RAD51C',
    'RAD51D',
    'RHO',
    'SCN5A',
    'SFPQ',
    'SGCB',
    'TARDBP',
    'TP53',
    'TPK1',
    'TSC2',
    'VHL',
    'XRCC2'
]

# Mapping of genes to phenotypes:
# gene: (case phenotypes, control phenotypes, inheritence)
gene_phenotypes = {
    'BAP1': (
        {
            'melanoma',
            'mesothelioma',
            'renal cancer',
            'basal cell carcinoma',
            'intracranial meningioma',
            'choloangiocarcinoma of the billiary tract'
        },
        {'cancer'},
        'autosomal dominant'
    ),
    'BARD1': (
        {'breast cancer'},
        {'cancer'},
        'autosomal dominant'
    ),
    'BRCA1': (
        {
            'breast cancer',
            'ovarian cancer',
            'prostate cancer',
            'pancreatic cancer'
        },
        {'cancer'},
        'autosomal dominant'
    ),
    'BRCA2': (
        {
            'breast cancer',
            'ovarian cancer',
            'prostate cancer',
            'pancreatic cancer'
        },
        {'cancer'},
        'autosomal dominant'
    ),
    'CALM1': (
        {'long QT syndrome'},
        {'abnormal QT interval'},
        'autosomal dominant'
    ),
    'CALM2': (
        {'long QT syndrome'},
        {'abnormal QT interval'},
        'autosomal dominant'
    ),
    'CALM3': (
        {'long QT syndrome'},
        {'abnormal QT interval'},
        'autosomal dominant'
    ),
    #'CARD11',
    #'CBS',
    'CHEK2': (
        {
            'breast cancer',
            'prostate cancer'
        },
        {'cancer'},
        'autosomal dominant'
    ),
    #'CRX',
    #'CTCF',
    #'DDX3X',
    #'F9',
    #'FKRP',
    'G6PD': (
        {'deficiency of glucose-6-phosphate dehydrogenase'},
        {'deficiency of glucose-6-phosphate dehydrogenase'},
        'X-linked dominant'
    ),
    'GCK': (
        {'MODY'},
        {'diabetes'},
        'autosomal dominant'
    ),
    #'HMBS',
    #'JAG1',
    'KCNE1': (
        {'long QT syndrome'},
        {'abnormal QT interval'},
        'autosomal dominant'
    ),
    'KCNH2': (
        {'long QT syndrome'},
        {'abnormal QT interval'},
        'autosomal dominant'
    ),
    'KCNQ4': (
        {'nonsyndromic genetic hearing loss'},
        {'hearing loss'},
        'autosomal dominant'
    ),
    #'LARGE1',
    'MSH2': (
        {
            'colorectal cancer',
            'ovarian cancer',
            'stomach cancer',
            'small bowel cancer',
            'urinary tract cancer',
            'prostate cancer',
            'brain cancer',
            'endometrial cancer'
        },
        {'cancer'},
        'autosomal dominant'
    ),
    #'NDUFAF6',
    'OTC': (
        {'ornithine carbamoyltransferase deficiency'},
        {'disorder of the urea cycle metabolism', 'hyperammonemia'},
        'X-linked' # Males severely affected from birth, females asymptomatic to severe
    ),
    'PALB2': (
        {
            'breast cancer',
            'ovarian cancer',
            'pancreatic cancer'
        },
        {'cancer'},
        'autosomal dominant'
    ),
    #'PAX6',
    'PTEN': (
        {
            'breast cancer',
            'endometrial cancer',
            'thyroid cancer',
            'colorectal cancer',
            'renal cancer',
            'melanoma',
        },
        {'cancer'},
        'autosomal dominant'
    ),
    'RAD51C': (
        {
            'breast cancer',
            'ovarian cancer',
        },
        {'cancer'},
        'autosomal dominant'
    ),
    'RAD51D': (
        {
            'breast cancer',
            'ovarian cancer',
        },
        {'cancer'},
        'autosomal dominant'
    ),
    #'RHO',
    'SCN5A':(
        {
            'long QT syndrome',
            'dilated cardiomyopathy'
        },
        {
            'abnormal QT interval',
            'cardiomyopathy'
        },
        'autosomal dominant'
    ),
    #'SFPQ',
    #'SGCB',
    'TARDBP': (
        {
            'amyotrophic lateral sclerosis',
            'frontotemporal dementia'
        },
        {
            'amyotrophic lateral sclerosis',
            'frontotemporal dementia'
        },
        'autosomal dominant'
    ),
    'TP53': (
        {'cancer'},
        {'cancer'},
        'autosomal dominant'
    ),
    #'TPK1',
    'TSC2': (
        {'tuberous sclerosis'},
        {'tuberous sclerosis'},
        'autosomal dominant'
    ),
    #'VHL',
    #'XRCC2'
}

# Bergquist et al. thresholds 10.1016/j.gim.2025.101402

# AlphaMissense 

# Benign Moderate+: <= 0.070
# Benign Moderate: <= 0.099
# Benign Supporting: <= 0.169
# Pathogenic Supporting: >= 0.792
# Pathogenic Moderate: >= 0.906
# Pathogenic Moderate+: >= 0.972
# Pathogenic Strong: >= 0.990

# MutPred2

# REVEL