"""Script for submitting to the IGVF portal"""
# This script requires the igvf_utils package and assumes that the submitter's
# API keys are available in IGVF_API_KEY and IGVF_SECRET_KEY environmenr variables

if __name__ != '__main__':
    raise 'This is a script. Do not import.'

from igvf_utils.connection import Connection

common_properties = {
    "lab": "mark-craven",
    "award": "HG012039",
}

payloads = [
    {
        Connection.PROFILE_KEY: 'prediction_set'
    },
    {
        Connection.PROFILE_KEY: 'workflow',
        'aliases': ["mark-craven:cvfg-aou-workflow-v1"],
        'name': 'Biobank validation of variant classifications',
        'source_url': 'https://github.com/Craven-Biostat-Lab/igvf-cvfg-biobank'
    },
    {
        Connection.PROFILE_KEY: 'analysis_step',
        'aliases': ["mark-craven:cvfg-aou-analysis-step-v1"],
        "title": 'Biobank validation of variant classifications',
        "analysis_step_types": [], # Need to create new types?
        "input_content_types": [], # Not sure about matching types
        "output_content_types": [], # Not sure about matching types
        "step_label": "esm-1v-substitution-scoring-step",
        "workflow": "mark-craven:esm1v-wokflow-v1"
    },
    {
        Connection.PROFILE_KEY: 'software',
        'aliases': ["mark-craven:cvfg-aou-software-v1"],
        "title": "Biobank validation of variant classifications",
        "name": "igvf-cvfg-biobank",
        "source_url": 'https://github.com/Craven-Biostat-Lab/igvf-cvfg-biobank',
        "description": "Python package for performing biobank validation of variant classifications in the All of Us Workbench."
    },
    {
        Connection.PROFILE_KEY: 'software_version',
        'aliases': ["mark-craven:cvfg-aou-software-version-v1"],
        'download_id': None,
        'version': None,
        'software': 'mark-craven:cvfg-aou-software-v1'
    },
    {
        Connection.PROFILE_KEY, 'tabular_file'
    },
]

conn = Connection('dev', dry_run=True)

for payload in payloads:
    payload.update(common_properties)
    conn.post(payload)


