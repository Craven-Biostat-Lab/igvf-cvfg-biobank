version 1.0

task collect_participant_covariates {
    # Maybe should be direct query to AoU tables?
    #input {
    #    File ancestry_pca
    #    File demographics_df
    #}

    command <<<
        collect-aou-covariates covariates.parquet
    >>>

    output {
        File covariates = "covariates.parquet"
    }
}

task identify_carrier_groups {
    input {
        File variant_groups
    }

    command <<<
        pip install git+https://github.com/Craven-Biostat-Lab/igvf-cvfg-biobank.git@cromwell[hail]
        identify-carrier-groups --variants ~{variant_groups} --carriers carrier_groups.parquet
    >>>

    output {
        File carrier_groups = "carrier_groups.parquet"
    }

    runtime {
        docker: "spark:python3"
    }
}

task estimate_odds_ratios {
    input {
        File case_control_cohort
        File covariates
        File carrier_groups
    }

    command <<<
        python estimate-odds-ratios \
            --cohort ~{case_control_cohort} \
            --covariates ~{covariates} \
            --carriers ~{carrier_groups} \
            --odds-ratios odds_ratio_df.parquet
    >>>

    output {
        File odds_ratio_df = "odds_ratio_df.parquet"
    }
}

workflow single_dataset_analysis {
    input {
        File case_control_cohort
        File variant_groups
    }

    call collect_participant_covariates {}

    call identify_carrier_groups {
        input:
            variant_groups = variant_groups
    }

    call estimate_odds_ratios {
        input:
            case_control_cohort = case_control_cohort,
            covariates = collect_participant_covariates.covariates,
            carrier_groups = identify_carrier_groups.carrier_groups
    }

    output {
        File odds_ratio_estimates = estimate_odds_ratios.odds_ratio_df
    }
}