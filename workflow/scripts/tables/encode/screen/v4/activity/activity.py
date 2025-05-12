"""rdhs-activity.py"""

import numpy as np
import pandas as pd

# Snakemake
CCRES = snakemake.input[0]  # type: ignore
SIGNAL_SUMMARY = snakemake.input[1]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore


# ------------- #
# Functions     #
# ------------- #


def read_rdhslist(filepath: str) -> pd.DataFrame:
    """Returns CCRE bed file as DataFrame"""
    return pd.read_csv(
        filepath,
        sep="\t",
        header=None,
        engine="c",
        usecols=[3],
        names=["rDHS"],
        dtype={"rDHS": str},
    )["rDHS"].to_list()


def read_signal_summary(filepath: str) -> pd.DataFrame:
    dtype = {
        "rDHS": str,
        "sum_signal": float,
        "num_signal": int,
    }
    return pd.read_csv(
        filepath,
        sep="\t",
        dtype=dtype,
        engine="c",
    )


def main():
    """Main"""
    # Read inputs
    rdhs_list = read_rdhslist(CCRES)
    signal_summary = read_signal_summary(SIGNAL_SUMMARY)

    # Subset to rDHS to known ccre
    signal_summary = signal_summary[signal_summary["rDHS"].isin(rdhs_list)]

    # Calculate activity
    signal_summary["activity"] = signal_summary["sum_signal"] / np.sqrt(
        signal_summary["num_signal"]
    )

    # Add quantile across all rDHS
    signal_summary["quantile_rdhswide"] = pd.qcut(
        signal_summary["activity"], 100, labels=[i for i in range(1, 101)]
    )

    # Tidy up
    signal_summary = signal_summary[["rDHS", "activity", "quantile_rdhswide"]].round(4)
    dtypes = {"rDHS": str, "activity": float, "quantile_rdhswide": int}
    signal_summary = signal_summary.astype(dtypes)

    # Write output
    signal_summary.to_csv(OUTPUT, sep="\t", index=False)
