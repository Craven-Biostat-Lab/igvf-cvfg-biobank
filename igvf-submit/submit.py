"""Script for submitting to the IGVF portal"""
# This script requires the igvf_utils package and assumes that the submitter's
# API keys are available in IGVF_API_KEY and IGVF_SECRET_KEY environmenr variables

if __name__ != '__main__':
    raise 'This is a script. Do not import.'

DEBUG = False
if DEBUG:
    print("Debug / dry-run mode")

import pandas as pd
from igvf_utils.connection import Connection
import igvf_utils as iu

common_properties = {
    "lab": "mark-craven",
    "award": "HG012039",
}

# Pull gene list from submission file
estimates_file_path = 'data/all-or-estimates_2026-02-13.csv.gz'
ensg_list = pd.read_csv(estimates_file_path, usecols=['ENSG'])['ENSG'].drop_duplicates().to_list()
ensg_list = [s.split('.')[0] for s in ensg_list]

payloads = [
    { # Checked 2026-01-02
        Connection.PROFILE_KEY: 'document',
        'aliases': ["mark-craven:cvfg-aou-or-estimates-documentation-v1"],
        'document_type': 'file format specification',
        'description': 'File format specification for Biobank validation OR estimates.',
        'attachment': {'path': 'OR Estimates File Description.pdf'}
    },
    { # Checked 2026-02-13
        Connection.PROFILE_KEY: 'prediction_set',
        'aliases': ["mark-craven:cvfg-aou-or-estimates-fileset-v1"],
        'description':
            'Statistical estimates of the ratio of odds of condition occurrence given that '
            'a person carries at least one of the variants in a given class.',
        'donors': ["igvf:virtual_human_donor"],
        'input_file_sets': [
            "IGVFDS4364ECDX", # Calibration sets
            "IGVFDS4011IBUI" # Assay variant sets
        ],
        'file_set_type': 'functional effect',
        'scope': 'genes',
        'small_scale_gene_list': ensg_list
    },
    { # Checked 2026-01-02
        Connection.PROFILE_KEY: 'analysis_step',
        'aliases': ["mark-craven:cvfg-aou-analysis-step-v1"],
        "title": 'Biobank validation of variant classifications',
        "analysis_step_types": ["logistic regression"],
        "input_content_types": ['variant effects'],
        "output_content_types": ['pathogenicity validation'],
        "step_label": "cvfg-aou-step"
    },
    { # Checked 2026-01-02
        Connection.PROFILE_KEY: 'software',
        'aliases': ["mark-craven:cvfg-aou-software-v1"],
        "title": "Biobank validation of variant classifications",
        "name": "igvf-cvfg-biobank",
        "source_url": 'https://github.com/Craven-Biostat-Lab/igvf-cvfg-biobank',
        "description":
            "Python package for performing biobank validation of variant classifications "
            "in the All of Us Workbench."
    },
    { # Checked 2026-02-13
        Connection.PROFILE_KEY: 'software_version',
        'aliases': ["mark-craven:cvfg-aou-software-version-2026-02-13"],
        'download_id': 'https://github.com/Craven-Biostat-Lab/igvf-cvfg-biobank/releases/tag/v1.1.0',
        'version': 'v1.1.0',
        'software': 'mark-craven:cvfg-aou-software-v1'
    },
    { # Checked 2026-02-13
        Connection.PROFILE_KEY: 'analysis_step_version',
        'aliases': ['mark-craven:cvfg-aou-avalysis-step-version-2026-02-13'],
        'analysis_step': 'mark-craven:cvfg-aou-analysis-step-v1',
        'software_versions': ['mark-craven:cvfg-aou-software-version-2026-02-13']
    },
    { # Checked 2025-12-08
        Connection.PROFILE_KEY: 'workflow',
        'aliases': ["mark-craven:cvfg-aou-workflow-v1"],
        'name': 'Biobank validation of variant classifications',
        'source_url': 'https://github.com/Craven-Biostat-Lab/igvf-cvfg-biobank',
        'analysis_step_versions': ['mark-craven:cvfg-aou-avalysis-step-version-v1']
    },
    { # Checked 2026-02-13
        Connection.PROFILE_KEY: 'tabular_file',
        "aliases": ["mark-craven:cvfg-aou-or-estimates-v3"],
        "description": "Statistical estimates of the ratio of odds of condition occurrence given that a person carries at least one of the variants in a given class.",
        "analysis_step_version": "mark-craven:cvfg-aou-avalysis-step-version-2026-02-13",
        "content_type": "pathogenicity validation",
        "controlled_access": False,
        "derived_from": [
            "IGVFFI1443TQDN", # Aggregated calibrations (has predictor and combined points)
            "IGVFFI2521UGYG", # Assay variant scores table (has author scores)
            "IGVFFI7610PCPU" # Assay variant calibrations
        ],
        "file_format": "csv",
        "file_format_specifications": ["mark-craven:cvfg-aou-or-estimates-documentation-v1"],
        "file_set": "mark-craven:cvfg-aou-or-estimates-fileset-v1",
        "submitted_file_name": estimates_file_path
    },
]

conn = Connection('prod', dry_run=DEBUG)

for payload in payloads:
    payload.update(common_properties)
    print(payload)
    print(conn.get_profile_from_payload(payload).properties)
    if DEBUG:
        try:
            conn.post(payload)
        except Exception as e:
            print(e)
    else:
        conn.post(payload)
