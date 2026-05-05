version 1.0

task assemble_case_control_cohort {
    # Should the input here be SQL?
    input {
        Array[File] case_pools
        Array[File] control_pools
    }
    output {
        case_control_cohort = "cohort.parquet"
    }
}

task collect_participant_covariates {
    # Maybe should be direct query to AoU tables?
    input {
        File ancestry_pca
        File demographics_df
    }
    output {
        File covariates
    }
}

task identify_carrier_groups {
    input {
        File variant_groups
    }
    output {
        File carrier_groups
    }
}

task estimate_odds_ratios {
    input {
        File case_control_cohort
        File covariates
        File carrier_groups
    }

    command <<<
        python estimate_odds_ratios
    >>>

    output {
        File odds_ratio_df = "odds_ratio_df.parquet"
    }
}

workflow single_dataset_analysis {
    input {
        
    }

    call assemble_case_control_cohort {}

    call collect_participant_covariates {}

    call identify_carrier_groups {}

    call estimate_odds_ratios {
        input:
            case_control_cohort = assemble_case_control_cohort.case_control_cohort,
            covariates = collect_participant_covariates.covariates,
            carrier_groups = identify_carrier_groups.carrier_groups
    }

    output {
        File odd_ratio_estimates = estimate_odds_ratios.odds_ratio_df
    }
}