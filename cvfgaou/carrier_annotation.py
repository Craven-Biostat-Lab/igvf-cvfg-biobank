"""Variant Carrier Annotations"""

import logging
from itertools import chain

import pandas as pd
import numpy as np
from cvfgaou import hailtools, gctools, notation

log = logging.getLogger(__name__)

class CarrierAnnotatorVEP:
    """Class for performing carrier status annotations for AlphaMissense"""

    def __init__(
        self,
        variant_level_calibrations_df,
        gene_thresholds_df,
        output_location,
        wgs_mt,
        clinvar_bins_df,
        vat_loader, # Function that takes a gene and returns a dataframe
        variant_table_col_map,
        predictor_name,
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

        self.predictor_name = predictor_name
        self.outfile_prefix = predictor_name.lower()

        self.cols = variant_table_col_map

        # Load gene-specific Calibration table
        # This table is indexed by gene
        self.gene_thresholds_df = gene_thresholds_df

        # Mapping from our thresholds to thresholds used in the table and the respective comparison direction
        self.gs_threshold_map = {
            f'{inequality} {sign}{points}': (f'{criterion}_{strength.title()}', comparator)
            for inequality, sign, criterion, _, comparator, _ in notation.DOE_TABLE
            for points, strength in notation.SOE_TABLE
        }

    def build_exposure_package(self):

        for gene, per_gene_df in self.progress_tracker(self.variant_level_calibrations_df.groupby('gene_symbol')):
            
            exposures_file = f'{self.output_location}/exposures/{self.outfile_prefix}_{gene}.parquet'
            clinvar_file = f'{self.output_location}/clinvar_maps/{self.outfile_prefix}_{gene}.parquet'
            af_file = f'{self.output_location}/af_maps/{self.outfile_prefix}_{gene}.parquet'
            
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
                left_on=[
                    self.cols['contig'],
                    self.cols['position'],
                    self.cols['ref_allele'],
                    self.cols['alt_allele']
                ],
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
            if gene in self.gene_thresholds_df.index:
                variant_classes = chain(
                    variant_classes,
                    (
                        ('Calibrated (gene-specific)', c, per_gene_df[selection])
                        for c, selection in (
                            (classification, compare(per_gene_df[self.cols['score']], threshold))
                            for classification, (label, compare) in self.gs_threshold_map.items()
                            for threshold in (self.gene_thresholds_df.loc[gene, label], )
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
                        contig_col=self.cols['contig'],
                        pos_col=self.cols['position'],
                        ref_col=self.cols['ref_allele'],
                        alt_col=self.cols['alt_allele'],
                        metadata_dict = {
                            'Dataset': self.predictor_name,
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


class CarrierAnnotatorAlphaMissense(CarrierAnnotatorVEP):
    """Class for performing carrier status annotations for AlphaMissense"""

    def __init__(
        self,
        variant_level_calibrations_df,
        gene_thresholds_df,
        output_location,
        wgs_mt,
        clinvar_bins_df,
        vat_loader,
        splice_ai_filter_max=np.inf,
        af_filter_max=np.inf,
        progress_tracker=None
    ):
        super().__init__(
            variant_level_calibrations_df,
            gene_thresholds_df,
            output_location,
            wgs_mt,
            clinvar_bins_df,
            vat_loader,
            {
                'gene': 'gene_symbol',
                'contig': '#CHROM',
                'position': 'POS',
                'ref_allele': 'REF',
                'alt_allele': 'ALT',
                'score': 'am_pathogenicity'
            },
            'AlphaMissense',
            splice_ai_filter_max,
            af_filter_max,
            progress_tracker
        )


class CarrierAnnotatorREVEL(CarrierAnnotatorVEP):
    """Class for performing carrier status annotations for REVEL"""

    def __init__(
        self,
        variant_level_calibrations_df,
        gene_thresholds_df,
        output_location,
        wgs_mt,
        clinvar_bins_df,
        vat_loader,
        splice_ai_filter_max=np.inf,
        af_filter_max=np.inf,
        progress_tracker=None
    ):
        super().__init__(
            variant_level_calibrations_df,
            gene_thresholds_df,
            output_location,
            wgs_mt,
            clinvar_bins_df,
            vat_loader,
            {
                'gene': 'gene_symbol',
                'contig': '#CHROM',
                'position': 'POS',
                'ref_allele': 'REF',
                'alt_allele': 'ALT',
                'score': 'am_pathogenicity'
            },
            'AlphaMissense',
            splice_ai_filter_max,
            af_filter_max,
            progress_tracker
        )