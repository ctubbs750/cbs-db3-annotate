"""Ruleset description"""

from pathlib import Path
from snakemake.utils import min_version

# ----------- #
#  Settings   #
# ----------- #

min_version("7.32.4")


# ------------- #
# Config        #
# ------------- #

# IO DIRS
PROCESS_DIR = Path(config["FILEPATHS"]["PROCESS_DIRS"]["TABLES"])
BENCHMARK_DIR = Path(config["FILEPATHS"]["BENCHMARK_DIR"])
LOG_DIR = Path(config["FILEPATHS"]["LOG_DIR"])

# IO
GENCODE_GTF = Path(config["RESOURCES"]["GENCODE"]["V40"]["GTF"])

# ------------- #
#     Rules     #
# ------------- #


rule all:
    """Target rule"""
    input:
         PROCESS_DIR / "gencode.v40.basic.annotation.parsed.pc.longest_tx.tsv.table.gz"


# ------------- #
#    Gencode    #
# ------------- #


rule gencode_gene_table:
    """
    Parses Gencode GTF into delimited table with 1:1 gene to transcript mapping. 
    """
    input:
        GENCODE_GTF
    output:
        PROCESS_DIR / "gencode.v40.basic.annotation.parsed.pc.longest_tx.tsv.table.gz"
    benchmark:
        BENCHMARK_DIR / "gencode_gene_table.benchmark"
    log:
        LOG_DIR / "gencode_gene_table.log"
    # envmodules:
    #     # List of modules here
    # params:
    #     # Parameters here
    resources:
        mem_mb=3000,
        runtime=5,
    threads: 2
    conda:
        "install"
    script:
        "../scripts/tables/gencode/gencode.py"