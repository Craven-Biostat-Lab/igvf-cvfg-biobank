"""Variant Carrier Annotations"""

import logging
from itertools import chain

import pandas as pd
import numpy as np
from cvfgaou import hailtools, gctools, notation

log = logging.getLogger(__name__)


class VariantGrouper:
    """Interface for generating variant groups for analysis"""

    def __init__(self):
        self.classifier_name = None
        self.bundles = set()
        self.cols = {
            'contig': 'contig',
            'position': 'position',
            'ref_allele': 'ref_allele',
            'alt_allele': 'alt_allele'
        }
    
    def metadata(self, dataset, gene):
        """Metadata dictionary for this bundle"""
        return {}
    
    def get_groups(self, dataset, gene):
        return []
    
    @property
    def variant_cols(self):
        return [
            self.cols['contig'],
            self.cols['position'],
            self.cols['ref_allele'],
            self.cols['alt_allele']
        ]


class TableVariantGrouper(VariantGrouper):

    def __init__(
        self,
        variants_df,
        gene_col,
        classifier_name,
        dataset = None,
        dataset_col = None,
        col_mapping = None,
        metadata_cols = None,
        fixed_metadata = None
    ):
        super().__init__()

        self.df = variants_df
        self.classifier_name = classifier_name
        self.metadata_cols = metadata_cols
        self.fixed_metadata = fixed_metadata

        self.cols = {
            'contig': 'contig',
            'position': 'position',
            'ref_allele': 'ref_allele',
            'alt_allele': 'alt_allele'
        }

        if col_mapping is not None:
            self.cols.update(col_mapping)
        
        self.cols['gene'] = gene_col
        self.cols['dataset'] = dataset_col

        self.dataset = dataset
        if self.dataset is not None:
            if dataset_col is not None:
                raise ValueError("Provide exactly one of dataset or dataset_col")
                        
            self.bundles = {
                (self.dataset, gene) for gene in self.df[self.cols['gene']].drop_duplicates()
            }
        else:
            if dataset_col is None:
                raise ValueError("Provide exactly one of dataset or dataset_col")
            
            self.bundles = {
                (dataset, gene)
                for dataset, gene
                in self.df[[self.cols['dataset'], self.cols['gene']]].drop_duplicates().itertuples(index=False)
            }
    
    def metadata(self, dataset, gene):
        metadata = {}
        if self.fixed_metadata is not None:
            metadata.update(self.fixed_metadata)
        if self.metadata_cols is not None:
            subframe = self.get_subframe(dataset, gene)
            if subframe is not None and not subframe.empty:
                for col in self.metadata_cols:
                    if col in subframe.columns:
                        info_val = subframe[col].iloc[0]
                        if not (subframe[col] == info_val).all():
                            log.warning(f'Multiple values for {col} in {dataset}, {gene}; Using the first.')
                        metadata[col] = info_val
                    else:
                        log.warning(f'Metadata column {col} not in dataframe')
        return metadata

    def get_subframe(self, dataset, gene):

        if self.dataset is not None:
            if dataset != self.dataset:
                return None
        
        subframe = self.df[
            (self.df[self.cols['gene']] == gene)
        ]
        if self.cols['dataset'] is not None:
            subframe = subframe[subframe[self.cols['dataset']] == dataset]
        
        return subframe

class TableClassVariantGrouper(TableVariantGrouper):

    def __init__(
        self,
        variants_df,
        class_col,
        gene_col,
        classifier_name=None,
        dataset=None,
        dataset_col=None,
        col_mapping=None,
        metadata_cols=None,
        fixed_metadata=None
    ):
        
        super().__init__(
            variants_df,
            gene_col,
            class_col if classifier_name is None else classifier_name,
            dataset,
            dataset_col,
            col_mapping,
            metadata_cols,
            fixed_metadata
        )

        self.cols['class'] = class_col

    def get_groups(self, dataset, gene):

        subframe = self.get_subframe(dataset, gene)
        if subframe is None or subframe.empty:
            return

        return subframe.groupby(self.df[self.cols['class']])


class TablePointsVariantGrouper(TableVariantGrouper):

    def __init__(
        self,
        variants_df,
        points_col,
        gene_col,
        classifier_name=None,
        dataset=None,
        dataset_col=None,
        points_limits=(None,None),
        col_mapping=None,
        metadata_cols=None,
        fixed_metadata=None
    ):
        
        super().__init__(
            variants_df,
            gene_col,
            points_col if classifier_name is None else classifier_name,
            dataset,
            dataset_col,
            col_mapping,
            metadata_cols,
            fixed_metadata
        )

        self.cols['points'] = points_col
        self.points_min, self.points_max = points_limits

    def get_groups(self, dataset, gene):
        
        subframe = self.get_subframe(dataset, gene)
        if subframe is None or subframe.empty:
            return

        # Establish upper and lower bounds
        points_min = int(subframe[self.cols['points']].min())
        points_max = int(subframe[self.cols['points']].max())
        if self.points_min is not None and points_min > self.points_min:
            points_min = self.points_min
        if self.points_max is not None and points_max < self.points_max:
            points_min = self.points_max

        for points in range(points_min, points_max+1):
            if points < 0:
                yield f'{notation.LEQ_CHAR} {points:+d}', subframe.loc[
                    subframe[self.cols['points']] <= points,
                    self.variant_cols
                ]
            elif points > 0:
                yield f'{notation.GEQ_CHAR} {points:+d}', subframe.loc[
                    subframe[self.cols['points']] >= points,
                    self.variant_cols
                ]
            else: # points == 0 or NaN
                yield '0', subframe.loc[
                    subframe[self.cols['points']] == 0,
                    self.variant_cols
                ]


class TableScoresVariantGrouper(TableVariantGrouper):

    def __init__(
        self,
        variants_df,
        score_col,
        gene_col,
        classifier_name=None,
        dataset=None,
        dataset_col=None,
        fixed_thresholds=None,
        gene_thresholds=None,
        dataset_thresholds=None,
        col_mapping=None,
        metadata_cols=None,
        fixed_metadata=None
    ):
        
        super().__init__(
            variants_df,
            gene_col,
            score_col if classifier_name is None else classifier_name,
            dataset,
            dataset_col,
            col_mapping,
            metadata_cols,
            fixed_metadata
        )

        self.cols['score'] = score_col

        # Save thresholds
        if fixed_thresholds is not None and gene_thresholds is None and dataset_thresholds is None:
            self.threshold_type = 'fixed'
            self.thresholds = fixed_thresholds
        elif fixed_thresholds is None and gene_thresholds is not None and dataset_thresholds is None:
            self.threshold_type = 'gene'
            self.thresholds = gene_thresholds
        elif fixed_thresholds is None and gene_thresholds is None and dataset_thresholds is not None:
            self.threshold_type = 'dataset'
            self.thresholds = dataset_thresholds
        else:
            raise ValueError("Exactly one of fixed_ gene_ or dataset_ thresholds may be provided.")

    def get_groups(self, dataset, gene):

        subframe = self.get_subframe(dataset, gene)
        if subframe is None or subframe.empty:
            return

        if self.threshold_type == 'fixed':
            thresholds = self.thresholds
        if self.threshold_type == 'gene':
            thresholds = self.thresholds.get('gene')
        if self.threshold_type == 'dataset':
            thresholds = self.thresholds.get('dataset')
        
        if thresholds is not None:
            for label, compare, threshold in thresholds:
                variants = subframe.loc[compare(subframe[self.cols['score']], threshold), self.variant_cols]
                if not variants.empty:
                    yield label, variants


class CarrierAnnotator:
    """Class for performing carrier status annotations"""

    def __init__(
        self,
        output_location,
        wgs_mt,
        clinvar_bins_df,
        vat_loader, # Function that takes a gene and returns a dataframe
        variant_groupers=None,
        restrict_to_genes=None,
        progress_tracker=None
    ):
        self.output_location = output_location
        self.load_vat = vat_loader

        self.wgs_mt = wgs_mt
        self.clinvar_bins_df = clinvar_bins_df

        self.progress_tracker = (lambda x: x) if progress_tracker is None else progress_tracker

        self.variant_groupers = [] if variant_groupers is None else variant_groupers

        self.restrict_to_genes = restrict_to_genes
    

    def data_bundles(self):
        bundle_set = set()
        for grouper in self.variant_groupers:
            bundle_set |= grouper.bundles
        
        if self.restrict_to_genes is not None:
            bundle_set = {
                (dataset, gene)
                for dataset, gene in bundle_set
                if gene in self.restrict_to_genes
            }

        return bundle_set


    def build_exposure_package(self):

        for (dataset, gene) in self.progress_tracker(self.data_bundles()):

            bundle_name = f'{dataset}_{gene}'
            
            exposures_file = f'{self.output_location}/exposures/{bundle_name}.parquet'
            clinvar_file = f'{self.output_location}/clinvar_maps/{bundle_name}.parquet'
            af_file = f'{self.output_location}/af_maps/{bundle_name}.parquet'
            
            if gctools.blob_exists(exposures_file):
                log.info(f'Skipping {gene} because {exposures_file} exists.')
                continue
            
            gene_result_dfs = []
            clinvar_classes_dfs = []
            joint_af_map = {}

            for grouper in self.variant_groupers:
                classifier = grouper.classifier_name

                for classification, variant_df in self.progress_tracker(grouper.get_groups(dataset, gene)):      

                    if variant_df.empty: continue
                    
                    try:
                        exposure_df, af_map, clinvar_df = hailtools.get_exposure_package(
                            variant_df,
                            self.wgs_mt,
                            self.clinvar_bins_df,
                            contig_col=grouper.cols['contig'],
                            pos_col=grouper.cols['position'],
                            ref_col=grouper.cols['ref_allele'],
                            alt_col=grouper.cols['alt_allele'],
                            metadata_dict = {
                                'Dataset': dataset,
                                'Gene': gene,
                                'Classifier': classifier,
                                'Classification': classification
                            } | grouper.metadata(dataset, gene)
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
