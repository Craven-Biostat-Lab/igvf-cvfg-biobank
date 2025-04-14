"""Utilities for odds ratio estimation"""

import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm

def estimate_logOR(exposure_series, cohort_df, alpha=0.05):
    """Estimate log odds ratio given exposure and cohort

    exposure_series: A series of person IDs

    cohort_df: A dataframe indexed by person ID containing columns indicating case,
        sex_at_birth, age_at_cdr, and ancestry pca dimensions (lebeled pca_*)
    
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

    return pd.DataFrame(
        {
            'LogOR': estimate,
            'LogOR_LI': ci[0],
            'LogOR_UI': ci[1],
            'cases_with_variants': (model_df.exposure & model_df.case).sum(),
            'controls_with_variants': (model_df.exposure & ~model_df.case).sum(),
            'cases_without_variants': (~model_df.exposure & model_df.case).sum(),
            'controls_without_variants': (~model_df.exposure & ~model_df.case).sum()
        },
        index=[0]
    )