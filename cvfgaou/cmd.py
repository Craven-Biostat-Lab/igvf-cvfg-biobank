"""CLI that don't require optional packages"""

from argparse import ArgumentParser
from cvfgaou.odds import estimate_logOR
import pandas as pd
import os
import ast

def collect_aou_covariates(args=None):
    if args is None:
        parser = ArgumentParser(
            description="Collect covariates of interest for all AoU participants."
        )
        parser.add_argument('out')

        args = parser.parse_args()

    # Move these to args?
    dataset = os.environ["WORKSPACE_CDR"]
    use_bq = 'BIGQUERY_STORAGE_API_ENABLED' in os.environ
    ancestry_table = 'gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv'
    

    # Demographics and event count table for all WGS participants

    count_events_sql = f"""
        SELECT *
        FROM (
            SELECT
                person_id,
                COUNT(DISTINCT condition_start_date) AS n_events,
                MIN(condition_start_date) AS first_event,
                MAX(condition_start_date) AS last_event
            FROM `{dataset}.condition_occurrence`
            WHERE
                person_id IN (
                    SELECT person_id FROM `{dataset}.cb_search_person`
                    WHERE has_whole_genome_variant = 1
                )
            GROUP BY person_id
        ) INNER JOIN (
            SELECT
                person_id,
                sex_at_birth,
                dob,
                age_at_cdr
            FROM `{dataset}.cb_search_person`
        ) USING (person_id)
    """

    count_events_df = pd.read_gbq(
        count_events_sql,
        index_col='person_id',
        dialect='standard',
        use_bqstorage_api=use_bq
    ).assign(
        first_event = lambda df: pd.to_datetime(df.first_event),
        last_event = lambda df: pd.to_datetime(df.last_event),
        dob = lambda df: pd.to_datetime(df.dob)
    )

    ancestry_df = pd.read_table(
        ancestry_table,
        storage_options={'requester_pays': True}
    )

    pca_df = pd.concat(
        [
            ancestry_df[['research_id']],
            ancestry_df.pca_features.apply(
                lambda x: pd.Series(ast.literal_eval(x))
            ).rename(
                columns=lambda c: f'pca_{c}'
            )
        ],
        axis='columns'
    )

    result_df = count_events_df.merge(
        pca_df,
        left_index=True,
        right_on='research_id' 
    )

    result_df.to_parquet(args.out)


def estimate_odds_ratios(args=None):

    if args is None:
        parser = ArgumentParser(
            description="Estimates odds ratios using LR"
        )
        parser.add_argument('--cohort', required=True)
        parser.add_argument('--covariates', required=True)
        parser.add_argument('--carriers', required=True)
        parser.add_argument('--odds-ratios', required=True)

        args = parser.parse_args()
    
    # Open input frames
    carriers_df = pd.read_parquet(args.carriers)

    # Grab cohort, join to covariates
    # Loop through carriers columns, join to dataframe
    # Compute ORs (each carrier column becomes an OR row)
    # Build ouptut table

    df = pd.read_parquet(args.cohort).merge(
        pd.read_parquet(args.covariates),
        left_index=True,
        right_index=True
    )

    pd.DataFrame.from_dict(
        {
            col:
                estimate_logOR(
                    carriers_df[carriers_df[col]].index.to_series,
                    df
                )
            for col in carriers_df.columns
        },
        orient='index'
    ).to_parquet(args.odds_ratios)