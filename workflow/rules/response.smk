rule process_response:
    input:
        packages="data/procdata/metadata/github_packages.tsv",
        sources=extracted(RESPONSE_SOURCES),
        script="workflow/scripts/process_response.R",
        name_parser="workflow/scripts/parse_compound_names.R",
        environment_guard="workflow/scripts/check_environment.R",
        common="workflow/scripts/common.R"
    output:
        "data/procdata/response.rds"
    shell:
        "Rscript --vanilla workflow/scripts/process_response.R"

rule annotate:
    input:
        packages="data/procdata/metadata/github_packages.tsv",
        response="data/procdata/response.rds",
        sources=extracted(["cellosaurus"] + (["cmp_models"] if DATASET.startswith("GDSC") else [])),
        script="workflow/scripts/annotate.R",
        annotation_cache="workflow/scripts/annotation_cache.R",
        replay_checksums="config/annotationdb-checksums.json",
        environment_guard="workflow/scripts/check_environment.R",
        common="workflow/scripts/common.R"
    output:
        "data/procdata/annotations.rds"
    params:
        api=config["annotationdb"]
    shell:
        "Rscript --vanilla workflow/scripts/annotate.R"

rule fit:
    input:
        response="data/procdata/response.rds",
        script="workflow/scripts/fit.R",
        environment_guard="workflow/scripts/check_environment.R",
        common="workflow/scripts/common.R",
        packages="data/procdata/metadata/github_packages.tsv"
    output:
        "data/procdata/curve_response.rds"
    params:
        fitting=config["parameters"]["fitting"]
    shell:
        "Rscript --vanilla workflow/scripts/fit.R"

rule synergy:
    input:
        curves="data/procdata/curve_response.rds",
        script="workflow/scripts/synergy.R",
        environment_guard="workflow/scripts/check_environment.R",
        common="workflow/scripts/common.R",
        packages="data/procdata/metadata/github_packages.tsv"
    output:
        "data/procdata/fitted_response.rds"
    shell:
        "Rscript --vanilla workflow/scripts/synergy.R"
