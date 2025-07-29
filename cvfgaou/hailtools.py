"""Helper utilities for working with Hail"""

from cvfgaou.notation import var_str, var_col

import hail as hl
import logging

log = logging.getLogger(__name__)

def get_filtered_mt(
    variant_df, wgs_mt,
    contig_col,
    pos_col,
    ref_col,
    alt_col,
    reference_genome
):
    """Helper used in get_cols_with_variants and get_variant_sets_mt_row"""

    # Performance is much improved when the WGS table is first filtered to
    # intervals.

    intervals = [
        hl.locus_interval(
            contig=contig,
            start=int(df[pos_col].min()),
            end=int(df[pos_col].max()),
            includes_start=True,
            includes_end=True,
            reference_genome=reference_genome
        )
        for contig, df in variant_df.groupby(contig_col)
    ]

    # Make hail table for variants
    variant_ht = hl.Table.from_pandas(variant_df)
    variant_ht = variant_ht.annotate(
        locus = hl.locus(
            contig=variant_ht[contig_col],
            pos=variant_ht[pos_col],
            reference_genome=reference_genome
        ),
        alleles = [variant_ht[ref_col], variant_ht[alt_col]]
    )
    variant_ht = variant_ht.key_by(variant_ht.locus, variant_ht.alleles)

    # Filter down wgs_mt
    return hl.filter_intervals(wgs_mt, intervals).semi_join_rows(variant_ht)


def get_cols_with_variants(
    variant_df, wgs_mt,
    contig_col='Chromosome',
    pos_col='Genomic_coordinates',
    ref_col='Ref_allele',
    alt_col='Alt_allele',
    reference_genome='GRCh38'
):
    """Get the set of MatrixTable columns that have calls
    
    Given a dataframe describing variants and a MatrixTable with variant rows
    and sample columns (i.e. like the AoU WGS MT), output the columns that have
    the specified variants as a dataframe.
    """

    # Get person table
    filtered_mt = get_filtered_mt(
        variant_df=variant_df, wgs_mt=wgs_mt,
        contig_col=contig_col,
        pos_col=pos_col,
        ref_col=ref_col,
        alt_col=alt_col,
        reference_genome=reference_genome
    )
    
    return filtered_mt.filter_cols(
        hl.agg.any(filtered_mt.GT.is_non_ref())
    ).cols().to_pandas()


def get_detailed_cols_with_variants(
    variant_df, wgs_mt,
    contig_col='Chromosome',
    pos_col='Genomic_coordinates',
    ref_col='Ref_allele',
    alt_col='Alt_allele',
    reference_genome='GRCh38'
):
    """Get annotated set of MatrixTable columns that have calls
    
    Given a dataframe describing variants and a MatrixTable with variant rows
    and sample columns (i.e. like the AoU WGS MT), output the columns that have
    the specified variants as a dataframe.
    For each column, also annotate the list of variants from the class that the
    participant has, and whether any variant hits are homozygous.
    Additionally, return a dictionary specifying the allele frequency of each
    variant.

    Returns result_df, af_map
    """

    # Get person table
    filtered_mt = get_filtered_mt(
        variant_df=variant_df, wgs_mt=wgs_mt,
        contig_col=contig_col,
        pos_col=pos_col,
        ref_col=ref_col,
        alt_col=alt_col,
        reference_genome=reference_genome
    )
    
    mt_annotated = filtered_mt.annotate_cols(
        variants=hl.agg.filter(
            filtered_mt.GT.is_non_ref(),
            hl.agg.collect((
                filtered_mt.locus.contig,
                filtered_mt.locus.position,
                filtered_mt.alleles[0],
                filtered_mt.alleles[1],
                filtered_mt.info.AF[0]
            ))
        ),
        homozygous_variants=hl.agg.any(filtered_mt.GT.is_hom_var())
    )

    result_df = mt_annotated.filter_cols(
        hl.agg.any(mt_annotated.GT.is_non_ref())
    ).cols().to_pandas().set_index('s')

    # Convert to string for downstream compatibility
    af_map = {
        var_str(contig, pos, ref, alt): af
        for variants in result_df['variants']
        for contig, pos, ref, alt, af in variants
    }

    result_df['variants'] = result_df['variants'].apply(
        lambda l: [
            var_str(contig, pos, ref, alt)
            for contig, pos, ref, alt, _ in l
        ]
    )

    return result_df, af_map


def get_variant_set_mt_row(
    label_dict, variant_df, wgs_mt,
    contig_col,
    pos_col,
    ref_col,
    alt_col,
    reference_genome='GRCh38'
):
    
    # Filter down wgs_mt
    filtered_mt = get_filtered_mt(
        variant_df=variant_df, wgs_mt=wgs_mt,
        contig_col=contig_col,
        pos_col=pos_col,
        ref_col=ref_col,
        alt_col=alt_col,
        reference_genome=reference_genome
    ).annotate_rows(select=True)

    return (
        filtered_mt.
        group_rows_by('select').
        aggregate_entries(exposure=hl.agg.any(filtered_mt.GT.is_non_ref())).
        result().
        annotate_rows(**label_dict).
        key_rows_by(*label_dict.keys())
    )
 

def get_exposure_package(
    variant_df,
    wgs_mt,
    clinvar_bins_df,
    contig_col='Chromosome',
    pos_col='Genomic_coordinates',
    ref_col='Ref_allele',
    alt_col='Alt_allele',
    reference_genome='GRCh38',
    metadata_dict=None
):

    if variant_df.empty:
        log.warning("Empty variant set given")
        return None, {}, None
    
    clinvar_df = var_col(
        variant_df,
        contig_col,
        pos_col,
        ref_col,
        alt_col
    ).to_frame(name='Variant').merge(
        clinvar_bins_df,
        on='Variant',
        how='left'
    )
    clinvar_df['Clinvar significance'].fillna('Other / not in ClinVar', inplace=True)

    exposure_df, af_map = get_detailed_cols_with_variants(
        variant_df, wgs_mt, contig_col, pos_col, ref_col, alt_col, reference_genome
    )

    exposure_df.reset_index(names='person_id', inplace=True)

    # Attach metadata to frames
    if metadata_dict is not None:
        exposure_df=exposure_df.assign(**metadata_dict)
        clinvar_df=clinvar_df.assign(**metadata_dict)

    return exposure_df, af_map, clinvar_df