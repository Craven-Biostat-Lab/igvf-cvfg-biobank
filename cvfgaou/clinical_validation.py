import os
import numpy as np
import pandas as pd
import hail as hl
import json
import matplotlib.pyplot as plt
import seaborn as sns
from cvfgaou.hailtools import get_cols_with_variants

def generate_fig1_variant_counts(srwes_split_df, functional_dataset, phenotype_df, test_intervals=None, output_bucket=None):
    """
    Generate variant counts for Figure 1 visualization based on Findlay dataset and AoU WGS data.
    
    Parameters:
    -----------
    srwes_split_path : str
        Path to the WGS exome split Hail matrix table
    findlay_csv_path : str
        Path to the Findlay CSV file with variant annotations
    phenotype_df : pandas.DataFrame
        DataFrame containing phenotype data with person_id/s column
    test_intervals : list, optional
        List of genomic intervals to filter by (default: ['chr17:43044294-43125482'] for BRCA1)
    output_bucket : str, optional
        Output bucket path (uses env var WORKSPACE_BUCKET if None)
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with variant count statistics grouped by functional class
    """

    hl.init(default_reference='GRCh38', idempotent=True)
    
    # Set default values
    if test_intervals is None:
        test_intervals = ['chr17:43044294-43125482']  # this needs to change
        
    if output_bucket is None:
        output_bucket = os.getenv("WORKSPACE_BUCKET")
    
    # Filter WES data down to specified intervals
    mt_Findlay = hl.filter_intervals(
        srwes_split_df,
        [hl.parse_locus_interval(x) for x in test_intervals]
    )
    
    # 2. Reading functional CSV for further filtering - this is written based on Findlay dataset but should be adaptable to other datasets with similar structure
    # to-do: make a copy of the input df and process that instead of modifying the original
    functional_dataset['Chromosome'] = 'chr' + functional_dataset['Chrom'].astype(str)
    functional_dataset['locus'] = 'chr' + functional_dataset['Chrom'].astype(str) + ":" + functional_dataset['hg38_start'].astype(str)
    functional_dataset['alleles'] = functional_dataset.apply(lambda x: [x['ref_allele'], x['alt_allele']], axis=1)
    functional_dataset_subset = functional_dataset[[
        'Chromosome', 'hg38_start', 'ref_allele', 'alt_allele', 'locus',
        'alleles', 'auth_reported_score', 'auth_reported_func_class', 'Ensembl_transcript_ID'
    ]]
    hail_findlay = hl.Table.from_pandas(functional_dataset_subset)
    hail_findlay = hail_findlay.annotate(locus=hl.parse_locus(hail_findlay.locus))
    hail_findlay = hail_findlay.key_by('locus', 'alleles')
    
    filtered_mt = mt_Findlay.filter_rows(hl.is_defined(hail_findlay[mt_Findlay.row_key]))
    
    mt_findlay_annotated = filtered_mt.annotate_rows(
        findlay_metadata=hail_findlay[filtered_mt.locus, filtered_mt.alleles]
    )
    
    # 3. Process phenotype DataFrame
    df_persons_conditions = phenotype_df.copy()
    
    # Ensure the DataFrame has the expected structure
    if 'person_id' in df_persons_conditions.columns and 's' not in df_persons_conditions.columns:
        df_persons_conditions = df_persons_conditions.rename(columns={"person_id": "s"})

    ht_pheno = hl.Table.from_pandas(df_persons_conditions)
    ht_pheno = ht_pheno.annotate(s=hl.str(ht_pheno.s)).key_by('s')
    
    mt_filtered_ehr = mt_findlay_annotated.filter_cols(hl.is_defined(ht_pheno[mt_findlay_annotated.s]))
    
    # Attach phenotype data to column keys "s"
    mt_Findlay_with_pheno = mt_filtered_ehr.annotate_cols(phenotype_metadata=ht_pheno[mt_filtered_ehr.s])
    
    # Select only columns we need and flatten the MatrixTable into Table
    mt_select = mt_Findlay_with_pheno.select_rows('findlay_metadata')
    mt_flattened = mt_select.entries()
    
    # Export to workspace bucket in tsv.bgz format
    out_path = f"{output_bucket}/mt_fig1.tsv.bgz"
    mt_flattened.export(out_path)
    
    ## make this a separate function that takes in the path to the exported file and processes it for plotting
    # Read it back in using pandas
    df_plot = pd.read_csv(out_path, sep="\t", compression="gzip") 
    
    # Process for plotting
    df_plot_filtered = df_plot[~df_plot['GT'].isin(['0/0', '0|0'])]  # filter out non-carriers
    
    df_plot_filtered['findlay_metadata'] = df_plot_filtered['findlay_metadata'].apply(json.loads)
    df_plot_filtered['auth_reported_func_class'] = df_plot_filtered['findlay_metadata'].apply(
        lambda x: x['auth_reported_func_class']
    )
    
    # Extract alleles and create variant identifier
    df_plot_filtered[['allele1', 'allele2']] = df_plot_filtered['alleles'].str.extract(r'\["([^"]+)"\s*,\s*"([^"]+)"\]')
    df_plot_filtered['variant'] = df_plot_filtered['locus'] + ':' + df_plot_filtered['allele1'] + '_' + df_plot_filtered['allele2']
    
    # Create final dataframe for analysis
    df_plot_fig1 = df_plot_filtered[['variant', 's', 'auth_reported_func_class', 'GT']]
    
    # Calculate per-variant counts by functional class
    person_count_per_variant_func = df_plot_fig1.groupby(["variant", "auth_reported_func_class"])["s"].nunique().reset_index(name='count')
    person_count_per_variant_func['count_bin'] = np.where(
        person_count_per_variant_func["count"] >= 20, 
        ">20", 
        person_count_per_variant_func["count"]
    )
    
    # Create summary statistics for plotting
    stacked_plot_fig1 = person_count_per_variant_func[['auth_reported_func_class', 'count_bin', 'variant']]
    stacked_plot_fig1 = stacked_plot_fig1.groupby(['auth_reported_func_class', 'count_bin']).nunique('variant').reset_index()
    stacked_plot_fig1['percent_each_class'] = stacked_plot_fig1.groupby('auth_reported_func_class')['variant'].transform(lambda x: x/x.sum()*100)
    
    return stacked_plot_fig1


def plot_fig1_variant_counts(stacked_plot_df, plot_style="barplot", save_path=None):
    """
    Create a plot of variant counts by functional class.
    
    Parameters:
    -----------
    stacked_plot_df : pandas.DataFrame
        DataFrame with variant count statistics from generate_fig1_variant_counts
    plot_style : str, optional
        Plot style: "barplot" or "stacked" (default: "barplot")
    save_path : str, optional
        If provided, save the figure to this path
        
    Returns:
    --------
    matplotlib.axes.Axes
        The plot axes object
    """
    # Sort count_bin values numerically
    fig1_order = sorted([x for x in stacked_plot_df['count_bin'].unique() if x != ">20"], 
                        key=lambda x: float(x)) + ['>20']
    
    if plot_style == "barplot":
        fig, ax = plt.subplots(figsize=(12, 8))
        sns.barplot(
            data=stacked_plot_df,
            x='count_bin',
            y='percent_each_class',
            hue='auth_reported_func_class',
            order=fig1_order,
            ax=ax
        )
        plt.legend(
            title='Functional Class',
            bbox_to_anchor=(0.5, -0.15),
            loc="upper right",
            ncol=3
        )
    else:
        fig1_plot = stacked_plot_df.pivot(
            index='auth_reported_func_class',
            columns='count_bin',
            values='percent_each_class'
        ).fillna(0)
        
        # Ensure columns are in the right order
        fig1_plot = fig1_plot.reindex(columns=fig1_order)
        
        fig, ax = plt.subplots(figsize=(12, 8))
        ax = fig1_plot.plot(kind="bar", stacked=True, ax=ax)
        ax.legend(title='Count Bin', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.title("Variant Counts by Functional Class")
    plt.xlabel("Count Bin")
    plt.ylabel("Percent Each Class")
    
    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
    
    return ax


def generate_fig2_variant_counts(srwes_split_df, functional_dataset, phenotype_df, test_intervals=None, output_bucket=None):
    """
    Generate variant counts for Figure 1 visualization based on Findlay dataset and AoU WGS data.
    
    Parameters:
    -----------
    srwes_split_path : str
        Path to the WGS exome split Hail matrix table
    findlay_csv_path : str
        Path to the Findlay CSV file with variant annotations
    phenotype_df : pandas.DataFrame
        DataFrame containing phenotype data with person_id/s column
    test_intervals : list, optional
        List of genomic intervals to filter by (default: ['chr17:43044294-43125482'] for BRCA1)
    output_bucket : str, optional
        Output bucket path (uses env var WORKSPACE_BUCKET if None)
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with variant count statistics grouped by functional class
    """

    hl.init(default_reference='GRCh38', idempotent=True)
    
    # Set default values
    if test_intervals is None:
        test_intervals = ['chr17:43044294-43125482']  # this needs to change
        
    if output_bucket is None:
        output_bucket = os.getenv("WORKSPACE_BUCKET")

    # Filter WES data down to specified intervals
    mt_Findlay = hl.filter_intervals(
        srwes_split_df,
        [hl.parse_locus_interval(x) for x in test_intervals]
    )
    
    # 2. Reading functional CSV for further filtering - this is written based on Findlay dataset but should be adaptable to other datasets with similar structure
    # to-do: make a copy of the input df and process that instead of modifying the original
    functional_dataset['Chromosome'] = 'chr' + functional_dataset['Chrom'].astype(str)
    functional_dataset['locus'] = 'chr' + functional_dataset['Chrom'].astype(str) + ":" + functional_dataset['hg38_start'].astype(str)
    functional_dataset['alleles'] = functional_dataset.apply(lambda x: [x['ref_allele'], x['alt_allele']], axis=1)
    functional_dataset_subset = functional_dataset[[
        'Chromosome', 'hg38_start', 'ref_allele', 'alt_allele', 'locus',
        'alleles', 'auth_reported_score', 'auth_reported_func_class', 'Ensembl_transcript_ID'
    ]]

    # only people with findlay score variants
    persons_with_findlay_score = get_cols_with_variants(functional_dataset_subset, mt_Findlay, contig_col = 'Chromosome', pos_col = "hg38_start", ref_col = "ref_allele", alt_col = "alt_allele")

    ## phenotype processing
    dataset = %env WORKSPACE_CDR # ???
    death_data = pd.read_gbq(f'''SELECT distinct person_id, death_date, death_type_concept_id, cause_concept_id
                                , primary_death_record, src_id
                                FROM `{dataset}.aou_death` ''')
    