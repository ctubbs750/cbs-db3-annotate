"""Ruleset description"""

from pathlib import Path
from snakemake.utils import min_version

# ----------- #
#  Settings   #
# ----------- #

# min_version("7.32.4")

# ----------- #
#  Config     #
# ----------- #

# Parse config information here
PARAM = Path(config["PARAM"])

# ------------- #
#     Rules     #
# ------------- #


rule all:
    """Target rule"""
    input:
        # Path to target file here


# ------------- #
#     ----      #
# ------------- #


rule rule_X:
    """
    Notes on rule_X here
    """
    input:
        # Input file here,
    output:
        # Output file here
    benchmark:
        BENCHMARKS_DIR / "rule_X.benchmark"
    logs:
        LOGS_DIR / "rule_X.log"
    envmodules:
        # List of modules here
    params:
        # Parameters here
    resources:
        # Resources here
    conda:
        # Conda environment here
    shell:
        "# Shell command here, for example"