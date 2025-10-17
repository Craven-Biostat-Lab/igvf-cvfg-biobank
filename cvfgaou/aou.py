"""Utilities specific to how All of Us does things"""

from functools import cache
import pandas as pd

class CohortLoader():

    def __init__(
        self,
        gene_cohort_map,
        ancestry_df,
        demo_df,
        data_dir
    ):
        self.gene_cohort_map = gene_cohort_map
        self.ancestry_df = ancestry_df
        self.demo_df = demo_df
        self.data_dir = data_dir
        self.trichotomize_age = True
    

    @cache
    def gene_cohort(self, gene):

        case_dfs = [
            pd.read_table(
                f'{self.data_dir}/{case_cohort}',
                usecols=['person_id'],
                index_col='person_id'
            ).assign(case=1)
            for case_cohort in self.gene_cohort_map[gene][0]
        ]

        controls = None
        for control_cohort in self.gene_cohort_map[gene][1]:
            df = pd.read_table(f'{self.data_dir}/{control_cohort}',
                usecols=['person_id'],
                index_col='person_id'
            )
            if controls is None:
                controls = df
            else:
                controls = controls.join(df, how='inner')
        
        result = pd.concat(
            case_dfs + [controls.assign(case=0)]
        ).join(self.ancestry_df).join(self.demo_df)
        
        if self.trichotomize_age:
            result.sex_at_birth.where(
                result.sex_at_birth.isin(['Male', 'Female']),
                'Other',
                inplace=True
            )
        
        return result
