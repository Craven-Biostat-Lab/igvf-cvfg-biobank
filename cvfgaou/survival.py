from sksurv.nonparametric import kaplan_meier_estimator

def get_kaplan_meier_estimate(cohort_with_timepoints_df):

    time, prob_survival, conf_int = kaplan_meier_estimator(
        event=cohort_with_timepoints_df.case,
        time_exit=cohort_with_timepoints_df.timepoints,
        conf_type='log-log'
    )

    return pd.DataFrame(
        {
            'Age': time,
            'Survival to onset': prob_survival,
            'Survival_LI': conf_int[0],
            'Survival_UI': conf_int[1]
        }
    )