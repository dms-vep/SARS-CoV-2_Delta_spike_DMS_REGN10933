"""Top-level ``snakemake`` file that runs analysis."""


import os


configfile: "config.yaml"


# include `dms-vep-pipeline` pipeline Snakemake file
include: os.path.join(config["pipeline_path"], "pipeline.smk")


rule all:
    input:
        variant_count_files,
        rules.check_adequate_variant_counts.output.passed,
        antibody_escape_files,
        (
            [config["muteffects_observed"], config["muteffects_latent"]]
            if len(func_selections)
            else []
        ),
        config["docs"],


# Arbitrary other rules should be added here
rule site_numbering_map:
    """Map sequential numbering of protein in experiments to standard reference."""
    input:
        prot=config["gene_sequence_protein"],
    output:
        reference="results/site_numbering/numbering_reference.fa",
        alignment="results/site_numbering/alignment.fa",
        to_align="results/site_numbering/to_align.fa",
        site_numbering_map=config["site_numbering_map"],
    params:
        numbering_reference_accession=config["numbering_reference_accession"],
    log:
        os.path.join(config["logdir"], "site_numbering_map.txt"),
    conda:
        "dms-vep-pipeline/environment.yml"
    script:
        "scripts/site_numbering_map.py"

rule validation_IC50s:
    """Get ``polyclonal`` predicted IC50s for validated mutations."""
    input:
        config["validation_ic50s"],
        [
            os.path.join(config["escape_dir"], f"{antibody}.pickle")
            for antibody in pd.read_csv(config["validation_ic50s"])["antibody"].unique()
        ],
        nb="notebooks/validation_IC50s.ipynb",
    output:
        nb="results/notebooks/validation_IC50s.ipynb",
    log:
        os.path.join(config["logdir"], "validation_IC50s.txt"),
    conda:
        "dms-vep-pipeline/environment.yml"
    shell:
        "papermill {input.nb} {output.nb} &> {log}"


# Add any extra data/results files for docs with name: file
extra_data_files = {
    "sequential to reference site numbering": config["site_numbering_map"],
}


# include `dms-vep-pipeline` docs building Snakemake file
include: os.path.join(config["pipeline_path"], "docs.smk")
