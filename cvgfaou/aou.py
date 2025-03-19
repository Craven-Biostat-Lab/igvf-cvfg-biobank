# AoU Analysis Code
#
# Functions used for analysis on the All of Us researcher workbench.

import os
from datetime import date
from pathlib import Path
from functools import cache, cached_property
import json
import gzip

import pandas as pd
import hail as hl
import statsmodels.formula.api as smf
import statsmodels.api as sm


# General AoU helpers


def load_ancestry():
    """Load AoU ancestry table"""
    return pd.read_table(
        'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv',
        storage_options={'requester_pays': True}
    )


def save_json_gz(obj, path):
    with gzip.open(path, 'wt') as ostream:
        json.dump(obj, ostream)


# Our utility functions


def load_json_ht(path, key='hgvs_p'):
    return hl.Table.from_pandas(pd.read_json(path), key=key)


def load_intervals(path, gene_set):
    """Load interval definition table.
    
    path: path to file where a Biomart table with gene IDs and gene coordinates is stored.
    gene_set: the set of gene symbols to which to filter the table down
    """

    intervals_df = pd.read_table(
        path,
        usecols=[
            'Gene stable ID version',
            'Chromosome/scaffold name',
            'Gene start (bp)',
            'Gene end (bp)',
            'Gene name'
        ],
        dtype=str
    ).drop_duplicates()

    # Filter to genes of interest

    intervals_df = intervals_df[intervals_df['Gene name'].isin(gene_set)]

    # Cleanup

    intervals_df = intervals_df[intervals_df['Chromosome/scaffold name'] != 'HG2334_PATCH']

    # Build interval specs

    intervals_df['interval_spec'] = intervals_df.apply(
        lambda row: f"chr{row['Chromosome/scaffold name']}:{row['Gene start (bp)']}-{row['Gene end (bp)']}",
        axis='columns'
    )
    intervals_df


def get_mave_variant_counts(missense_mt, cohort_ht, mave_ht):
    
    # Filter only to mave variants
    filtered_mt = missense_mt.filter_rows(
        hl.is_defined(mave_ht.index(missense_mt.hgvs_p))
    )
    
    # For each variant (row), annotate if it has calls in cases or controls
    row_annotated_mt = filtered_mt.annotate_rows(
        case_variant=hl.agg.any(cohort_ht[filtered_mt.s].case & filtered_mt.GT.is_non_ref()),
        control_variant=hl.agg.any(cohort_ht[filtered_mt.s].control & filtered_mt.GT.is_non_ref()),
    )
    
    # For each participant, annotate the if it has calls in a variant
    col_annotated_mt = filtered_mt.annotate_cols(
        has_variant=hl.agg.any(filtered_mt.GT.is_non_ref())
    )
    
    return {
        # Number of mave variants present in cases and controls
        'case_variants': row_annotated_mt.aggregate_rows(hl.agg.sum(row_annotated_mt.case_variant)),
        'control_variants': row_annotated_mt.aggregate_rows(hl.agg.sum(row_annotated_mt.control_variant)),
        
        # Count number of cases and controls that have mave variants 
        'variant_cases': col_annotated_mt.aggregate_cols(hl.agg.sum(
            col_annotated_mt.has_variant & cohort_ht[col_annotated_mt.s].case
        )),
        'variant_controls': col_annotated_mt.aggregate_cols(hl.agg.sum(
            col_annotated_mt.has_variant & cohort_ht[col_annotated_mt.s].control
        ))
    }


# We use a class to share variables across analysis steps.
class AoU_OR_Analysis:
    """Main class used for OR analysis"""

    def __init__(
        self,
        genes_to_cohorts,
        mave_genes,
        bucket = os.getenv("WORKSPACE_BUCKET"),
        version = date.today().isoformat(),
        intervals_df = None,
        wgs_mt = None,
    ) -> None:
        
        self.genes_to_cohorts = genes_to_cohorts
        self.mave_genes = mave_genes
        self.bucket = bucket
        self.version = version

        # Options that override defaults:

        if intervals_df is not None:
            self.intervals_df = intervals_df
        
        if wgs_mt is not None:
            self.wgs_mt = wgs_mt

    
    @cached_property
    def intervals_df(self):
        return load_intervals(
            f'{self.bucket}/aux_data/gene_metadata.txt',
            self.genes_to_cohorts.keys()
        )

    
    @cached_property
    def wgs_mt(self):
        return hl.read_matrix_table(os.getenv('WGS_EXOME_SPLIT_HAIL_PATH'))


    @cache
    def filter_wgs_to_gene(self, gene_name):
        return hl.filter_intervals(
            self.wgs_mt,
            [
                hl.parse_locus_interval(interval, 'GRCh38')
                for interval
                in self.intervals_df[self.intervals_df['Gene name'] == gene_name]['interval_spec']
            ]
        )


    @cache
    def load_gene_vat(self, gene_name):
        vat_table = hl.import_table(
            f'{self.bucket}/aux_data/{gene_name.lower()}_vat.tsv',
            types={'position': hl.tint32}
        )

        vat_table = vat_table.annotate(
            locus=hl.locus(vat_table.contig, vat_table.position, 'GRCh38'),
            alleles=hl.array([vat_table.ref_allele, vat_table.alt_allele]),
            missense=vat_table.consequence.contains('missense')
        )

        # Use only missense
        vat_table = vat_table.filter(vat_table.missense, keep=True).key_by('locus', 'alleles')
        
        return vat_table


    @cached_property
    def demo_df(self):
        return pd.read_table(
            f'{self.bucket}/cohorts/wgs_demo.tsv.gz',
            parse_dates = ['first_event', 'last_event', 'dob']
        )


    @cache
    def get_pool(self, pool_name):
        return pd.read_table(f'{self.bucket}/cohorts/{pool_name}.tsv.gz')


    @cached_property
    def pca_df(self):
        return pd.concat(
            [
                self.ancestry_df[['research_id']].astype(str),
                self.ancestry_df.pca_features.apply(
                    lambda x: pd.Series(json.loads(x))
                ).rename(columns=lambda c: f'pca_{c}')
            ],
            axis='columns'
        )


    def build_cohort_mt(self, case_pool_df, control_pool_df, gene):
        
        cohort_pool_df = self.demo_df.assign(
            case = self.demo_df.person_id.isin(case_pool_df.person_id),
            control = self.demo_df.person_id.isin(control_pool_df.person_id),
        )
        
        # Remove short length of record from control pool, and remove non-cohort rows
        cohort_pool_df = cohort_pool_df[
            cohort_pool_df.case |
            (
                cohort_pool_df.control &
                ((cohort_pool_df.last_event.dt.year - cohort_pool_df.first_event.dt.year) >= 5)
            )
        ]
        
        # Trichotomize sex at birth
        cohort_pool_df.sex_at_birth.where(
            cohort_pool_df.sex_at_birth.isin(['Male', 'Female']),
            'Other',
            inplace=True
        )
        
        # Filter WGS columns to cohort
        cohort_ht = hl.Table.from_pandas(
            cohort_pool_df.drop(
                columns=['first_event', 'last_event', 'dob']
            ).astype({'person_id': str}),
            key='person_id'
        )
        
        # Load WGS and VAT, constrain to cohort
        gene_vat = self.load_gene_vat(gene)
        missense_mt = self.filter_wgs_to_gene(gene).semi_join_rows(gene_vat).semi_join_cols(cohort_ht)
        
        # Annotate WGS with clinvar classification
        clinvar_lists = gene_vat[missense_mt.locus, missense_mt.alleles].clinvar_classification.split(r',\s*')

        missense_mt = missense_mt.annotate_rows(
            clinvar_pathogenic = clinvar_lists.contains('pathogenic'),
            clinvar_likely_pathogenic = clinvar_lists.contains('pathogenic') | clinvar_lists.contains('likely pathogenic'),
            hgvs_p = gene_vat[missense_mt.locus, missense_mt.alleles].aa_change.replace(r'^.+p.\(([A-Za-z0-9]+)\)', 'p.$1'),
            revel = hl.parse_float(gene_vat[missense_mt.locus, missense_mt.alleles].revel)
        )
        
        return missense_mt, cohort_ht, gene_vat


    def estimate_OR_for_set(self, missense_mt, cohort_ht, selection, alpha=0.05):

        # Annotate case control and exposure data
        exposure_mt = missense_mt.annotate_cols(
            exposure = hl.agg.any(
                missense_mt.GT.is_non_ref() &
                selection
            )
        )

        exposure_df = cohort_ht.annotate(exposure = exposure_mt.cols()[cohort_ht.person_id].exposure).to_pandas()
        
        # Fix types
        exposure_df = exposure_df.astype({
            'age_at_cdr': int,
            'sex_at_birth': str,
            'case': int,
            'exposure': bool
        })

        if exposure_df.exposure.any():
        
            # Attach PCA
            exposure_df = exposure_df.merge(
                self.pca_df,
                how='left',
                left_on='person_id',
                right_on='research_id'
            )

            lr_model = smf.glm(
                formula='case ~ exposure + sex_at_birth + age_at_cdr + ' + ' + '.join(self.pca_df.columns[1:]),
                data=exposure_df,
                family=sm.families.Binomial()
            ).fit()

            estimate = lr_model.params['exposure[T.True]']
            ci = lr_model.conf_int(alpha=alpha).loc['exposure[T.True]']
        
            return estimate, ci[0], ci[1]
        
        else:
            return None, None, None


    def estimate_OR_for_collection(self, missense_mt, cohort_ht, classes, alpha=0.05):

        result = []

        for classifier, label, selection in classes:

            estimate, lower, upper = self.estimate_OR_for_set(missense_mt, cohort_ht, selection, alpha)

            if estimate is not None:
                result.append({
                    'Classifier': classifier,
                    'Classification': label,
                    'LogOR': estimate,
                    'LogOR_LI': lower,
                    'LogOR_UI': upper
                })


    def get_clinvar_logors(self, missense_mt, cohort_ht, alpha=0.05):
        
        classes = [
            ('ClinVar', 'pathogenic', missense_mt.clinvar_pathogenic),
            ('ClinVar', 'likely pathogenic', missense_mt.clinvar_likely_pathogenic)
        ]
        
        return self.estimate_OR_for_collection(missense_mt, cohort_ht, classes, alpha)


    # REVEL thresholds based on Pejaver et al. https://pubmed.ncbi.nlm.nih.gov/36413997/
    def get_revel_logors(self, missense_mt, cohort_ht, alpha=0.05):
        
        classes = [
            ('Calibrated', 'Benign Very Strong', missense_mt.revel <= 0.003),
            ('Calibrated', 'Benign Strong', missense_mt.revel <= 0.016),
            ('Calibrated', 'Benign Moderate', missense_mt.revel <= 0.183),
            ('Calibrated', 'Benign Supporting', missense_mt.revel <= 0.290),
            ('Calibrated', 'Pathogenic Supporting', missense_mt.revel >= 0.644),
            ('Calibrated', 'Pathogenic Moderate', missense_mt.revel >= 0.773),
            ('Calibrated', 'Pathogenic Strong', missense_mt.revel >= 0.932)
        ]
        
        return self.estimate_OR_for_collection(missense_mt, cohort_ht, classes, alpha)


    def get_mave_logors(self, missense_mt, cohort_ht, mave_ht, alpha=0.05):
        
        author_labels = mave_ht[missense_mt.hgvs_p].author_labels
        oob_pred = mave_ht[missense_mt.hgvs_p].oob_pred
        
        classes = [
            ('Author Label', 'Intermediate', hl.set({'Functionally_Abnormal', 'Intermediate'}).contains(author_labels)),
            ('Author Label', 'Functionally Abnormal', author_labels == 'Functionally_Abnormal'),
            ('Calibrated', 'Pathogenic Supporting', oob_pred >= 1),
            ('Calibrated', 'Pathogenic Moderate', oob_pred >= 2),
            ('Calibrated', 'Pathogenic Strong', oob_pred >= 4),
            ('Calibrated', 'Pathogenic Very Strong', oob_pred >= 8)
        ]
        
        return self.estimate_OR_for_collection(missense_mt, cohort_ht, classes, alpha)

    @cache
    def get_mcv(self, gene):
        """Missense matrixtable, cohort table, VAT cache"""
        return self.build_cohort_mt(
            self.get_pool[self.genes_to_cohorts[gene][0]],
            self.get_pool[self.genes_to_cohorts[gene][1]],
            gene
        )


    @cache
    def get_mave_ht(mave):
        """Mave cache"""
        return load_json_ht(Path('figs_10_09_24', mave, 'oob_df.json'), key='hgvs_p')


    def compute_and_run_all(self, output_dir):

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        for gene in self.genes_to_cohorts:
            
            save_target = output_path / f'clinvar-logors_{gene}_{self.version}.json.gz'
            if not save_target.exists():
            
                print(f"Computing ClinVar ORs for {gene}...")
                missense_mt, cohort_ht, _ = self.get_mcv(gene)
                clors = self.get_clinvar_logors(missense_mt, cohort_ht)

                print('Saving results...')
                save_json_gz(clors, save_target)

        # Generate revel estimates
        for gene in self.genes_to_cohorts:
            
            save_target = output_path / f'revel-logors_{gene}_{self.version}.json.gz'
            if not save_target.exists():
            
                print(f"Computing REVEL ORs for {gene}...")
                missense_mt, cohort_ht, _ = self.get_mcv(gene)
                rlors = self.get_revel_logors(missense_mt, cohort_ht)

                print('Saving results...')
                save_json_gz(rlors, save_target)

        # Generate MAVE estimates
        for mave, gene in self.mave_genes.items():
            if gene in self.genes_to_cohorts:
                
                save_target = output_path / f'mave-logors_{gene}_{mave}_{self.version}.json.gz'
                if not save_target.exists():
            
                    print(f"Computing MAVE logORs for {mave} ({gene})")
                    missense_mt, cohort_ht, _ = self.get_mcv(gene)
                    mave_ht = self.get_mave_ht(mave)
                    mlors = self.get_mave_logors(missense_mt, cohort_ht, mave_ht)
                    
                    print('Saving results...')
                    save_json_gz(mlors, save_target)
                        
                counts_target = output_path / f'mave-variant-counts_{gene}_{mave}_{self.version}.json.gz'
                if not Path(counts_target).exists():
                    
                    print(f'Counting variants for {mave} ({gene})')
                    missense_mt, cohort_ht, _ = self.get_mcv(gene)
                    mave_ht = self.get_mave_ht(mave)
                    variant_counts = get_mave_variant_counts(missense_mt, cohort_ht, mave_ht)
                    
                    print('Saving results...')
                    save_json_gz(variant_counts, counts_target)