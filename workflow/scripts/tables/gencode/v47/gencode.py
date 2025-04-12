"""Creates Gencode gene database with 1:1 mapping from gene to transcript"""

import pandas as pd
from gtfparse import read_gtf

# Snakemake
# GENCODE_IP = snakemake.input[0]  # type: ignore
GENCODE_IP = "/panfs/accrepfs.vampire/data/ruderferlab/projects/ctcf/resources/data/gencode/v47/gencode.v47.basic.annotation.gtf"
# GENCODE_OP = snakemake.output[0]  # type: ignore
GENCODE_OP = "/data/ruderferlab/projects/ctcf/projects/cbs-db3/cbs-db3-annotate/results/tables/gencode/gencode.v47.basic.annotation.parsed.pc.longest_tx.tsv.gz"


# ------------- #
# Functions     #
# ------------- #


def print_process():
    """Print process information."""
    print("Creating GENCODE gene datatable with 1:1 mapping from gene to transcript")
    print("-------------------------------------------------------------")
    print("Input GTF file: ", GENCODE_IP)
    print("Output TSV file: ", GENCODE_OP)
    print("-------------------------------------------------------------")


def parse_gencode() -> pd.DataFrame:
    """Parse Gencode GTF file to extract gene and transcript information."""
    # Read GTF file
    return read_gtf(GENCODE_IP).to_pandas()


def main():
    # Print process information
    print_process()

    # Parse Gencode GTF file
    gencode = parse_gencode()

    # Filter for protein-coding transcripts
    fields = [
        "transcript_id",
        "gene_id",
        "gene_name",
        "gene_type",
        "seqname",
        "start",
        "end",
    ]
    gencode_tx = gencode[
        (gencode["feature"] == "transcript")
        & (gencode["gene_type"] == "protein_coding")
    ][fields]
    gencode_tx["length"] = gencode_tx["end"] - gencode_tx["start"]

    # Select longest transcript per gene
    gencode_tx.sort_values(
        by=["gene_id", "length"], ascending=[True, False], inplace=True
    )
    gencode_tx.drop_duplicates(subset=["gene_id"], keep="first", inplace=True)
   
    # Trim trailing "." from gene_id and tx_id
    gencode_tx.loc[:, "transcript_id"] = gencode_tx["transcript_id"].str.split(".").str[0]
    gencode_tx.loc[:, "gene_id"] = gencode_tx["gene_id"].str.split(".").str[0]

    # Sort to tidy and subset
    gencode_tx.sort_values(by=["seqname", "start", "end"], inplace=True)
    gencode_tx[["gene_id", "gene_name", "transcript_id"]].to_csv(
        GENCODE_OP, sep="\t", index=False, compression="gzip"
    )


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
