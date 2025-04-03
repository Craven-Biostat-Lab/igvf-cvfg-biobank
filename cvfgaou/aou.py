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
    

    @cache
    def gene_cohort(self, gene):
        return pd.concat(
            [
                pd.read_table(
                    f'{self.data_dir}/{cohort}',
                    usecols=['person_id'],
                    index_col='person_id'
                ).assign(case=case)
                for case, cohorts in zip([1,0], self.gene_cohort_map[gene])
                for cohort in cohorts
            ]
        ).join(self.ancestry_df).join(self.demo_df)
