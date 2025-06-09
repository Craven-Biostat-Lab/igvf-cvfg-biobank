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
    
    Returns a 1-row dataframe with columns LogOR, LogOR_LI, LogOR_UI
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

    if variants_series is not None:
        # Assemble variant statistics
        # N. variants in cases/controls/both, AF range
        case_vs = variants_series[exposure_series.isin(cohort_df[cohort_df['case'] == 1].index)]
        control_vs = variants_series[exposure_series.isin(cohort_df[cohort_df['case'] == 0].index)]
        result_dict['variants_per_case'] = case_vs.apply(len).mean()
        result_dict['variants_per_control'] = control_vs.apply(len).mean()
        # Make case, control, and shared variant sets.
        # We don't use sets for this because dicts aren't hashable
        case_variants = []
        control_variants = []
        total_variants = []
        for v_list in case_vs:
            for variant in v_list:
                if variant not in case_variants:
                    case_variants.append(variant)
                if variant not in total_variants:
                    total_variants.append(variant)
        for v_list in control_vs:
            for variant in v_list:
                if variant not in control_variants:
                    control_variants.append(variant)
                if variant not in total_variants:
                    total_variants.append(variant)
        case_only_variants = [v for v in case_variants if v not in control_variants]
        control_only_variants = [v for v in control_variants if v not in case_variants]
        overlap_variants = [
            v for v in total_variants
            if (v in case_variants) and (v in control_variants)]
        
        result_dict['case_variant_min_af'] = min(v['af'] for v in case_variants) if case_variants else None
        result_dict['control_variant_min_af'] = min(v['af'] for v in control_variants) if control_variants else None
        result_dict['case_variant_max_af'] = max(v['af'] for v in case_variants) if case_variants else None
        result_dict['control_variant_max_af'] = max(v['af'] for v in control_variants) if control_variants else None
        result_dict['case_only_variant_count'] = len(case_only_variants)
        result_dict['control_only_variant_count'] = len(control_only_variants)
        result_dict['overlap_variant_count'] = len(overlap_variants)

    return pd.DataFrame(
        result_dict,
        index=[0]
    )