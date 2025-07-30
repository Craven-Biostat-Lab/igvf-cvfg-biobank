"""Utilities for odds ratio estimation"""

import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm

def estimate_logOR(exposure_series, cohort_df, alpha=0.05, variants_series=None):
    """Estimate log odds ratio given exposure and cohort

    exposure_series: A series of person IDs

    cohort_df: A dataframe indexed by person ID containing columns indicating case,
        sex_at_birth, age_at_cdr, and ancestry pca dimensions (lebeled pca_*)

    alpha: confidence level for confidence intervals

    variant_series (optional): A series of lists of dicts identifying variants
        for each person ID. If present, the output will include variant statistics.
    
    Returns a dictionary listing
        LogOR, LogOR_LI, LogOR_UI,
        and a breakdown of cases/controls with/without variants.
    """

    # Annotate cohort with exposure feature
    model_df = cohort_df.assign(exposure = cohort_df.index.isin(exposure_series))

    pca_cols = cohort_df.columns[cohort_df.columns.str.startswith('pca_')]

    lr_model = smf.glm(
        formula='case ~ exposure + sex_at_birth + age_at_cdr + ' + ' + '.join(pca_cols),
        data=model_df,
        family=sm.families.Binomial()
    ).fit()

    estimate = lr_model.params['exposure[T.True]']
    ci = lr_model.conf_int(alpha=alpha).loc['exposure[T.True]']

    result_dict = {
        'LogOR': estimate,
        'LogOR_LI': ci[0],
        'LogOR_UI': ci[1],
        'cases_with_variants': (model_df.exposure & model_df.case).sum(),
        'controls_with_variants': (model_df.exposure & ~model_df.case).sum(),
        'cases_without_variants': (~model_df.exposure & model_df.case).sum(),
        'controls_without_variants': (~model_df.exposure & ~model_df.case).sum()
    }

    return result_dict


def collect_variant_stats(exposure_series, cohort_df, variants_series, af_map, clinvar_class_map):
    """Collect variant statistics
    
    exposure_series: A series of person IDs

    cohort_df: A dataframe indexed by person ID containing columns indicating case,
        sex_at_birth, age_at_cdr, and ancestry pca dimensions (lebeled pca_*)

    variant_series:
        A series of lists of strings identifying variants for each person ID.
        Variant strings are formatted contig:pos:ref>alt.

    af_map: A mapping from variant string to Allele Frequency
        
    clinvar_class_map: A mapping from variant string (formatted contig:pos:ref>alt)
        to Clinvar significance. The set of variants in the keys of this map is
        assumed to be the set of all variants in the "exposure" variant class.
            
    Returns dictionary of variant statistics
    """

    result_dict = {}

    # Assemble variant statistics
    # N. variants in cases/controls/both
    case_vs = variants_series[exposure_series.isin(cohort_df[cohort_df['case'] == 1].index)]
    control_vs = variants_series[exposure_series.isin(cohort_df[cohort_df['case'] == 0].index)]
    result_dict['variants_per_case'] = case_vs.apply(len).mean()
    result_dict['variants_per_control'] = control_vs.apply(len).mean()

    # Make case, control, and shared variant sets.
    case_variants = {variant for v_list in case_vs for variant in v_list}
    control_variants = {variant for v_list in control_vs for variant in v_list}
    total_variants = case_variants | control_variants
    overlap_variants = case_variants & control_variants
    case_only_variants = case_variants - overlap_variants
    control_only_variants = control_variants - overlap_variants
    
    result_dict['case_variant_min_af'] = min(af_map[v] for v in case_variants) if case_variants else None
    result_dict['control_variant_min_af'] = min(af_map[v] for v in control_variants) if control_variants else None
    result_dict['case_variant_max_af'] = max(af_map[v] for v in case_variants) if case_variants else None
    result_dict['control_variant_max_af'] = max(af_map[v] for v in control_variants) if control_variants else None
    result_dict['case_only_variant_count'] = len(case_only_variants)
    result_dict['control_only_variant_count'] = len(control_only_variants)
    result_dict['overlap_variant_count'] = len(overlap_variants)

    # Variant count breakdowns, in class vs in AoU, and by ClinVar class
    clinvar_df = pd.Series(clinvar_class_map).to_frame(name='ClinVar')
    clinvar_df['Cohort'] = clinvar_df.index.isin(total_variants)

    result_dict['Variants in class'] = clinvar_df['ClinVar'].count()
    result_dict['Variants in cohort'] = clinvar_df['Cohort'].sum()

    for clinvar_class, count in clinvar_df['ClinVar'].value_counts().items():
        result_dict[f'Class ClinVar {clinvar_class}'] = count
    
    for clinvar_class, count in clinvar_df.loc[clinvar_df['Cohort'], 'ClinVar'].value_counts().items():
        result_dict[f'Class ClinVar {clinvar_class}'] = count
    
    return result_dict