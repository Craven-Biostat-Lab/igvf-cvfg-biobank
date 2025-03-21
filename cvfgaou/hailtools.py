"""Helper utilities for working with Hail"""

import hail as hl

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
    filtered_mt = hl.filter_intervals(wgs_mt, intervals).semi_join_rows(variant_ht)

    # Get person table
    return filtered_mt.filter_cols(
        hl.agg.any(filtered_mt.GT.is_non_ref())
    ).cols().to_pandas()