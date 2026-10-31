rule molecular:
    input:
        response="data/procdata/response.rds",
        sources=extracted(MOLECULAR_SOURCES),
        script="workflow/scripts/molecular.R",
        genes="workflow/scripts/gene_helpers.R",
        packages="data/procdata/metadata/github_packages.tsv",
        environment_guard="workflow/scripts/check_environment.R",
        common="workflow/scripts/common.R"
    output:
        "data/procdata/molecular_profiles.rds"
    shell:
        "Rscript --vanilla workflow/scripts/molecular.R"
