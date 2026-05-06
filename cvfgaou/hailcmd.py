"""Entrypoints that require Hail"""

from argparse import ArgumentParser
import pandas as pd
import hail as hl
from cvfgaou.hailtools import get_detailed_cols_with_variants

def identify_carrier_groups(args=None):

    if args is None:
        parser = ArgumentParser(
            prog="Identify Carrier Groups",
            description="Given a table of variants with scores and a WGS matrix, extracts the participants that have the annotated variants."
        )
        parser.add_argument('--wgs', required=True)
        parser.add_argument('--variants', required=True)
        parser.add_argument('--carriers', required=True)

        args = parser.parse_args()
    
    wgs_mt = hl.read_matrix_table(args.wgs)
    variants_df = pd.read_table(args.variants)

    carriers_df, af_map = get_detailed_cols_with_variants(
        variant_df=variants_df,
        wgs_mt=wgs_mt
    )

    carriers_df.to_csv(args.carriers)

