# setup
import warnings
import os
import subprocess
import numpy as np
import pandas as pd
import hail as hl
import bokeh.io
from bokeh.io import output_notebook, show
from bokeh.layouts import gridplot
from bokeh.resources import INLINE
import re
from lifelines import KaplanMeierFitter
from cvfgaou.hailtools import get_cols_with_variants

## just for AoU purposes, can be removed
hl.init(default_reference="GRCh38", idempotent=True)
bucket = os.getenv("WORKSPACE_BUCKET")
genomic_location = os.getenv("CDR_STORAGE_PATH")
bokeh.io.output_notebook(INLINE)


## to-do: explore the options to write unittests
def generate_variant_counts(
    input_functional_dataset, phenotype_df=None, test_intervals=None, output_bucket=None
):
    """
    Generate variant counts for visualization based on Findlay dataset and AoU WGS data.

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

    ## to-do: assert columns exist in functional dataset
    # clarify how we should ask users to generate auth_reported_func_class
    expected_columns_from_functional_dataset = [
        "chromosome",
        "hg38_start",
        "ref_allele",
        "alt_allele",
        "auth_reported_func_class",
    ]
    # to-do: assert datatypes (e.g. chromosome should be '11')
    # assert input_functional_dataset.columns == expected_columns_from_functional_dataset, 'Input columns do not match expected naming'

    # infer columns
    input_functional_dataset["locus"] = (
        "chr"
        + input_functional_dataset["chromosome"].astype(str)
        + ":"
        + input_functional_dataset["hg38_start"].astype(str)
    )
    input_functional_dataset["alleles"] = input_functional_dataset.apply(
        lambda x: [x["ref_allele"], x["alt_allele"]], axis=1
    )

    hl.init(default_reference="GRCh38", idempotent=True)

    # Set default values
    if test_intervals is None:
        warnings.warn("No interval specified, deriving using hg38 column")
        min_intervals = input_functional_dataset["hg38_start"].min()
        max_intervals = input_functional_dataset["hg38_start"].max()
        chromosome = "chr" + str(input_functional_dataset["chromosome"].unique()[0])
        test_intervals = [f"{chromosome}:{min_intervals}-{max_intervals}"]

    if output_bucket is None:
        output_bucket = os.getenv("WORKSPACE_BUCKET")

    # Get WES split file from directory
    srwes_split_path = os.getenv("WGS_EXOME_SPLIT_HAIL_PATH")
    split_mt = hl.read_matrix_table(srwes_split_path)

    # Filter WES data down to specified intervals
    mt_filtered = hl.filter_intervals(
        split_mt, [hl.parse_locus_interval(x) for x in test_intervals]
    )

    # 2. Reading functional CSV for further filtering - this is written based on Findlay dataset but should be adaptable to other datasets with similar structure
    hail_functional_df = hl.Table.from_pandas(input_functional_dataset)
    hail_functional_df = hail_functional_df.annotate(
        locus=hl.parse_locus(hail_functional_df.locus)
    )
    hail_functional_df = hail_functional_df.key_by("locus", "alleles")

    mt_filtered = mt_filtered.filter_rows(
        hl.is_defined(hail_functional_df[mt_filtered.row_key])
    )

    mt_functional_annotated = mt_filtered.annotate_rows(
        functional_metadata=hail_functional_df[mt_filtered.locus, mt_filtered.alleles]
    )

    if phenotype_df is not None:
        # selects all participants with variants present in func dataset
        input_functional_dataset["chromosome"] = (
            "chr" + input_functional_dataset["chromosome"]
        )
        persons_with_findlay_score = get_cols_with_variants(
            input_functional_dataset,
            split_mt,
            contig_col="chromosome",
            pos_col="hg38_start",
            ref_col="ref_allele",
            alt_col="alt_allele",
        )

        # 3. Process phenotype DataFrame
        df_persons_conditions = phenotype_df.copy()

        phenotype_df = phenotype_df[
            phenotype_df["s"].isin(persons_with_findlay_score["s"])
        ]

        # Ensure the DataFrame has the expected structure
        if (
            "person_id" in df_persons_conditions.columns
            and "s" not in df_persons_conditions.columns
        ):
            df_persons_conditions = df_persons_conditions.rename(
                columns={"person_id": "s"}
            )

        ht_pheno = hl.Table.from_pandas(df_persons_conditions)
        ht_pheno = ht_pheno.annotate(s=hl.str(ht_pheno.s)).key_by("s")

        mt_filtered_ehr = mt_functional_annotated.filter_cols(
            hl.is_defined(ht_pheno[mt_functional_annotated.s])
        )

        # Attach phenotype data to column keys "s"
        mt_functional_with_pheno = mt_filtered_ehr.annotate_cols(
            phenotype_metadata=ht_pheno[mt_filtered_ehr.s]
        )

        # Select only columns we need and flatten the MatrixTable into Table
        mt_select = mt_functional_with_pheno.select_rows("functional_metadata")

    mt_select = mt_functional_annotated.select_rows("functional_metadata")
    mt_select.describe(widget=True)
    mt_flattened = mt_select.entries()

    # Export to workspace bucket in tsv.bgz format
    out_path = f"{output_bucket}/mt_fig1.tsv.bgz"
    mt_flattened.export(out_path)


def plot_fig_1(file_name):
    bucket = os.getenv("WORKSPACE_BUCKET")
    out_path = f"{bucket}/{file_name}"

    df_plot = pd.read_csv(out_path, sep="\t", compression="gzip")

    df_plot_filtered = df_plot[
        ~df_plot["GT"].isin(["0/0", "0|0"])
    ]  # filter out non-carriers

    df_plot_filtered["functional_metadata"] = df_plot_filtered[
        "functional_metadata"
    ].apply(json.loads)
    df_plot_filtered["auth_reported_func_class"] = df_plot_filtered[
        "functional_metadata"
    ].apply(lambda x: x["auth_reported_func_class"])
    df_plot_filtered[["allele1", "allele2"]] = df_plot_filtered["alleles"].str.extract(
        r'\["([^"]+)"\s*,\s*"([^"]+)"\]'
    )
    df_plot_filtered["variant"] = (
        df_plot_filtered["locus"]
        + ":"
        + df_plot_filtered["allele1"]
        + "_"
        + df_plot_filtered["allele2"]
    )
    df_plot_fig1 = df_plot_filtered[["variant", "s", "auth_reported_func_class", "GT"]]

    return df_plot_fig1


def plot_fig_2(file_name):
    bucket = os.getenv("WORKSPACE_BUCKET")
    out_path = f"{bucket}/{file_name}"
    df_plot = pd.read_csv(
        out_path, sep="\t", compression="gzip"
    )  # 31876051 rows, unfiltered
    df_plot[["allele1", "allele2"]] = df_plot["alleles"].str.extract(
        r'\["([^"]+)"\s*,\s*"([^"]+)"\]'
    )
    df_plot["variant"] = (
        df_plot["locus"] + ":" + df_plot["allele1"] + "_" + df_plot["allele2"]
    )

    df_fig2_filtered = df_plot[~df_plot["GT"].isin(["0/0", "0|0"])]  # 78601 rows
    df_fig2_filtered["phenotype_metadata"] = df_fig2_filtered[
        "phenotype_metadata"
    ].apply(json.loads)
    df_fig2_filtered["age_at_diagnosis"] = df_fig2_filtered["phenotype_metadata"].apply(
        lambda x: x["age_at_diagnosis"]
    )
    df_fig2_filtered["breast_or_ovarian_diagnosis"] = df_fig2_filtered[
        "phenotype_metadata"
    ].apply(lambda x: x["breast_or_ovarian_diagnosis"])
    df_fig2_filtered["sex_at_birth"] = df_fig2_filtered["phenotype_metadata"].apply(
        lambda x: x["sex_at_birth"]
    )

    df_fig2_filtered["findlay_metadata"] = df_fig2_filtered["findlay_metadata"].apply(
        load_strings
    )
    df_fig2_filtered["auth_reported_func_class"] = df_fig2_filtered[
        "findlay_metadata"
    ].apply(lambda x: x["auth_reported_func_class"] if x else None)
    df_fig2_filtered["auth_reported_score"] = df_fig2_filtered[
        "findlay_metadata"
    ].apply(lambda x: x["auth_reported_score"] if x else None)
    df_fig2_filtered["clinvar_class"] = df_fig2_filtered["findlay_metadata"].apply(
        lambda x: x["clinvar_class"] if x else None
    )

    df_fig2_plot = df_fig2_filtered[
        [
            "s",
            "age_at_diagnosis",
            "sex_at_birth",
            "variant",
            "auth_reported_func_class",
            "auth_reported_score",
            "clinvar_class",
            "breast_or_ovarian_diagnosis",
        ]
    ]

    filtered_df_plot_fig2 = (
        df_fig2_plot.groupby("s").apply(select_variant).reset_index(drop=True)
    )

    df_plot_females = filtered_df_plot_fig2[
        filtered_df_plot_fig2["sex_at_birth"] == "Female"
    ]
    df_plot_breast_cancer = df_plot_females[
        [
            "s",
            "age_at_diagnosis",
            "breast_or_ovarian_diagnosis",
            "auth_reported_func_class",
        ]
    ].drop_duplicates()
    df_plot_breast_cancer["breast_or_ovarian_diagnosis"] = df_plot_breast_cancer[
        "breast_or_ovarian_diagnosis"
    ].astype(bool)

    # LOF vs. unscored
    km_res_func = km_risk_analysis(
        df_plot_breast_cancer, "FUNC", "breast_or_ovarian_diagnosis"
    )
    km_res_func
    km_res_lof = km_risk_analysis(
        df_plot_breast_cancer, "LOF", "breast_or_ovarian_diagnosis"
    )
    km_res_lof


def select_variant(val):
    # priority 1 - LOF
    lof_row = val[val["auth_reported_func_class"] == "LOF"]
    if len(lof_row):
        return lof_row.iloc[[0]]

    # priority 2 - INT
    int_row = val[val["auth_reported_func_class"] == "INT"]
    if len(int_row):
        return int_row.iloc[[0]]

    # priority 3 - findlay score, but FUNC
    int_row = val[val["auth_reported_func_class"] == "FUNC"]
    if len(int_row):
        return int_row.iloc[[0]]

    # if none of above, take first row
    return val.iloc[[0]]


def load_strings(x):
    import json

    if isinstance(x, str):
        return json.loads(x)
    return None


def fld(x, key):
    return x.get(key) if isinstance(x, dict) else None


def export_to_tsv(hail_df, file_name):
    bucket = os.getenv("WORKSPACE_BUCKET")
    out_path = f"{bucket}/{file_name}.tsv.bgz"
    hail_df.export(out_path)


def set_diagnosis_flag(df, column_name, keywords_list, new_column):
    pattern = "|".join(re.escape(keyword) for keyword in keywords_list)
    df[new_column] = np.where(
        df[column_name].str.contains(pattern, case=False, na=False), 1, 0
    )

    return df


def generate_phenotype_df(df_persons_conditions_control, df_persons_conditions_case):
    ## some logic to massage the cohort data
    df_pheno = pd.concat([df_persons_conditions_control, df_persons_conditions_case])
    ## To get the kaplan-meier curve, we'll need age at diagnosis
    df_pheno["condition_start_datetime"] = df_pheno[
        "condition_start_datetime"
    ].dt.tz_localize(None)
    df_pheno["date_of_birth"] = df_pheno["date_of_birth"].dt.tz_localize(None)
    df_pheno["age_at_diagnosis"] = (
        (
            (df_pheno["condition_start_datetime"] - df_pheno["date_of_birth"]).dt.days
            / 365
        )
        .round()
        .astype("Int64")
    )
    df_pheno["s"] = df_pheno["s"].astype(str)
    df_pheno = df_pheno[
        [
            "s",
            "standard_concept_name",
            "standard_concept_code",
            "race",
            "ethnicity",
            "sex_at_birth",
            "age_at_diagnosis",
        ]
    ]
    df_pheno.head(5)


def kmf_test(df, group, event_column):
    from lifelines import KaplanMeierFitter

    g1 = df[df["auth_reported_func_class"] == group]
    g1["age_at_diagnosis"] = g1["age_at_diagnosis"].astype(int)
    print(g1)
    f = KaplanMeierFitter()
    f.fit(durations=g1["age_at_diagnosis"], event_observed=g1[event_column])

    return f


# total pop at risk=total LOF variants / number of cases of cancers
def km_risk_analysis(df, func, event_column):
    from lifelines import KaplanMeierFitter
    from lifelines.statistics import logrank_test
    import matplotlib.pyplot as plt

    g1 = df.query(f"auth_reported_func_class == '{func}'")

    kmf = KaplanMeierFitter()
    plt.figure(figsize=(8, 6))

    f = kmf.fit(
        durations=g1["age_at_diagnosis"], event_observed=g1[event_column], label=func
    )
    cumulative_risk = 1 - kmf.survival_function_
    plt.step(cumulative_risk.index, cumulative_risk[func], where="post", label=func)
    return f
