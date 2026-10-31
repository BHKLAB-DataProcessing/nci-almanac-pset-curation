rule curate:
    input:
        fitted="data/procdata/fitted_response.rds",
        source_response="data/procdata/response.rds",
        annotations="data/procdata/annotations.rds",
        molecular="data/procdata/molecular_profiles.rds",
        receipts=expand("data/rawdata/downloads/{source}/receipt.json", source=list(SOURCES)),
        naming_sources=extracted(["combocompoundnames_small","combocompoundset"]),
        naming_policy="config/naming-policy.json",
        naming_helpers=["workflow/scripts/almanac_source_names.R","workflow/scripts/named_pset_helpers.R","workflow/scripts/compatibility_v2.R","workflow/scripts/parse_compound_names.R"],
        script="workflow/scripts/curate.R",
        aggregate="workflow/scripts/aggregate_response.R",
        observed_index="workflow/scripts/construct_observed_tre.R",
        molecular_annotations="workflow/scripts/molecular_annotations.R",
        portable_annotations="workflow/scripts/portable_annotations.R",
        environment_guard="workflow/scripts/check_environment.R",
        common="workflow/scripts/common.R"
    output:
        pset=config["outputs"]["pharmacoset"],
        samples=config["outputs"]["tables"] + "/samples.tsv",
        treatments=config["outputs"]["tables"] + "/treatments.tsv",
        sensitivity=config["outputs"]["tables"] + "/sensitivity.tsv",
        synergy=config["outputs"]["tables"] + "/synergy.tsv",
        mono=config["outputs"]["tables"] + "/monotherapy_profiles.tsv",
        combo=config["outputs"]["tables"] + "/combination_profiles.tsv",
        molecular_summary=config["outputs"]["tables"] + "/molecular_summary.tsv",
        validation="data/procdata/metadata/object_validation.tsv"
    params:
        dataset=config["dataset"]
    shell:
        "Rscript --vanilla workflow/scripts/curate.R"
