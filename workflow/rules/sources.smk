rule setup_packages:
    input:
        "workflow/scripts/setup.R",
        "renv.lock",
        "pixi.lock"
    output:
        "data/procdata/metadata/github_packages.tsv"
    params:
        packages=config["packages"]
    shell:
        "Rscript --vanilla workflow/scripts/setup.R"

rule download:
    input:
        script="workflow/scripts/download.py"
    output:
        payload="data/rawdata/downloads/{source}/payload",
        receipt="data/rawdata/downloads/{source}/receipt.json"
    params:
        specification=lambda w: SOURCES[w.source]
    shell:
        "python workflow/scripts/download.py {wildcards.source:q}"

rule extract:
    input:
        payload="data/rawdata/downloads/{source}/payload",
        receipt="data/rawdata/downloads/{source}/receipt.json",
        script="workflow/scripts/extract.py"
    output:
        "data/procdata/extracted/{source}/manifest.tsv"
    shell:
        "python workflow/scripts/extract.py {wildcards.source:q}"
