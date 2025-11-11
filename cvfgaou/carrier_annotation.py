"""Variant Carrier Annotations"""

import logging
from itertools import chain

import pandas as pd
import numpy as np
from cvfgaou import hailtools, gctools, notation

log = logging.getLogger(__name__)

class CarrierAnnotatorAlphaMissense:
    """Class for performing carrier status annotations for AlphaMissense"""

    def __init__(
        self,
        variant_level_calibrations_df,
        am_gene_thresholds_df,
        output_location,
        wgs_mt,
        clinvar_bins_df,
        vat_loader, # Function that takes a gene and returns a dataframe
        splice_ai_filter_max=np.inf,
        af_filter_max=np.inf,
        progress_tracker=None
    ):
        self.variant_level_calibrations_df = variant_level_calibrations_df
        self.output_location = output_location
        self.load_vat = vat_loader

        self.wgs_mt = wgs_mt
        self.clinvar_bins_df = clinvar_bins_df

        self.splice_ai_filter_max = splice_ai_filter_max
        self.af_filter_max = af_filter_max

        self.progress_tracker = (lambda x: x) if progress_tracker is None else progress_tracker

        # Load gene-specific Calibration table
        # This table is indexed by gene
        self.am_gene_thresholds_df = am_gene_thresholds_df

        # Mapping from our thresholds to thresholds used in the table and the respective comparison direction
        self.gs_threshold_map = {
            f'{inequality} {sign}{points}': (f'{criterion}_{strength.title()}', comparator)
            for inequality, sign, criterion, _, comparator, _ in notation.DOE_TABLE
            for points, strength in notation.SOE_TABLE
        }

    def build_exposure_package(self):

        for gene, per_gene_df in self.progress_tracker(self.variant_level_calibrations_df.groupby('gene_symbol')):
            
            exposures_file = f'{self.output_location}/exposures/alphamissense_{gene}.parquet'
            clinvar_file = f'{self.output_location}/clinvar_maps/alphamissense_{gene}.parquet'
            af_file = f'{self.output_location}/af_maps/alphamissense_{gene}.parquet'
            
            if gctools.blob_exists(exposures_file):
                log.info(f'Skipping {gene} because {exposures_file} exists.')
                continue
            
            # Filter spliceAI score
            gene_vat = self.load_vat(gene)
            gene_vat = gene_vat[
                (gene_vat.gnomad_all_af <= self.af_filter_max) &
                ~(
                    (gene_vat.splice_ai_acceptor_gain_score > self.splice_ai_filter_max) |
                    (gene_vat.splice_ai_acceptor_loss_score > self.splice_ai_filter_max) |
                    (gene_vat.splice_ai_donor_gain_score > self.splice_ai_filter_max) |
                    (gene_vat.splice_ai_donor_loss_score > self.splice_ai_filter_max)
                )
            ]
            per_gene_df = per_gene_df.merge(
                gene_vat[['contig', 'position', 'ref_allele', 'alt_allele']].drop_duplicates(),
                how='inner',
                left_on=['#CHROM', 'POS', 'REF', 'ALT'],
                right_on=['contig', 'position', 'ref_allele', 'alt_allele']
            )
            
            gene_result_dfs = []
            clinvar_classes_dfs = []
            joint_af_map = {}
            
            variant_classes = (
                (
                    method_label,
                    f'{notation.GEQ_CHAR if points > 0 else notation.LEQ_CHAR} {points:+d}',
                    per_gene_df[
                        per_gene_df[colname].apply(
                            notation.evidence_meets_strength,
                            args=(points,)
                        )
                    ]
                )
                for points in range(-4, 4) if points != 0
                for method_label, colname in (
                    ('Calibrated (genome-wide)', 'old_call'),
                    ('Calibrated (cluster-based)', 'cluster_call')
                )
            )

            # Gene-specific thresholds
            if gene in self.am_gene_thresholds_df.index:
                variant_classes = chain(
                    variant_classes,
                    (
                        ('Calibrated (gene-specific)', c, per_gene_df[selection])
                        for c, selection in (
                            (classification, compare(per_gene_df['am_pathogenicity'], threshold))
                            for classification, (label, compare) in self.gs_threshold_map.items()
                            for threshold in (self.am_gene_thresholds_df.loc[gene, label], )
                            if not np.isnan(threshold)
                        )
                    )
                )
            
            for classifier, classification, variant_df in self.progress_tracker(variant_classes):

                if variant_df.empty: continue
                    
                try:
                    exposure_df, af_map, clinvar_df = hailtools.get_exposure_package(
                        variant_df,
                        self.wgs_mt,
                        self.clinvar_bins_df,
                        contig_col='#CHROM',
                        pos_col='POS',
                        ref_col='REF',
                        alt_col='ALT',
                        metadata_dict = {
                            'Dataset': 'AlphaMissense',
                            'Gene': gene,
                            'Classifier': classifier,
                            'Classification': classification,
                            'SpliceAI filter max': self.splice_ai_filter_max,
                            'AF filter max': self.af_filter_max
                        }
                    )

                    clinvar_classes_dfs.append(clinvar_df)
                    joint_af_map.update(af_map)
                    gene_result_dfs.append(exposure_df)
                    
                except:
                    log.error(f'Failed on {gene}, {classifier}, {classification}:')
                    log.error(variant_df)
                    raise
                
            if clinvar_classes_dfs:
                pd.concat(clinvar_classes_dfs, ignore_index=True).to_parquet(clinvar_file, index=False)
            if joint_af_map:
                pd.Series(joint_af_map).to_frame(name='AF').to_parquet(af_file)
            if gene_result_dfs:
                pd.concat(gene_result_dfs, ignore_index=True).to_parquet(exposures_file, index=False)
