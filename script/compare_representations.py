"""
compare_representations.py

3. Decide which (bin_size, statistic) combination
to use as the PRIMARY representation for Phase 1, using chr1-only
output from extract_signal_multiscale.py.

Criterion: a good primary representation should (a) have low
sensitivity to arbitrary choices (its binary-presence calls shouldn't
flip wildly if you nudge the bin size slightly) and (b) preserve
enough spatial resolution to be biologically meaningful. This script
just reports the numbers; the choice of primary representation is
yours to make and record before Phase 1 starts.
"""

from pathlib import Path
from itertools import combinations

import pandas as pd
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SIGNAL_DIR = PROJECT_ROOT / "results" / "multiscale_signal"
OUTPUT_FILE = PROJECT_ROOT / "results" / "representation_comparison.csv"

STATS = ["mean", "median", "max", "frac_above_bg"]


def load_all():
    files = sorted(SIGNAL_DIR.glob("*_signal.parquet"))
    records = []
    for f in files:
        # filename pattern: {accession}_{binsize}bp_signal.parquet
        parts = f.stem.split("_")
        bin_size = int(parts[-2].replace("bp", ""))
        accession = "_".join(parts[:-2])
        df = pd.read_parquet(f)
        records.append((accession, bin_size, df))
    return records


def main():
    records = load_all()
    accessions = sorted(set(r[0] for r in records))

    rows = []
    for accession in accessions:
        by_binsize = {r[1]: r[2] for r in records if r[0] == accession}
        bin_sizes = sorted(by_binsize.keys())

        for stat in STATS:
            # Correlate the same statistic across adjacent bin-size
            # choices, after aligning on genomic position via
            # nearest-bin membership at the coarser resolution.
            for a, b in combinations(bin_sizes, 2):
                df_a = by_binsize[a][["chromosome", "start", "end", stat]].copy()
                df_b = by_binsize[b][["chromosome", "start", "end", stat]].copy()

                # Aggregate the finer bin size up to the coarser one
                # (mean of the statistic) so they're comparable.
                fine, coarse = (df_a, df_b) if a < b else (df_b, df_a)
                fine_size, coarse_size = min(a, b), max(a, b)

                fine["coarse_bin"] = (fine["start"] // coarse_size)
                agg = fine.groupby(["chromosome", "coarse_bin"])[stat].mean().reset_index()
                coarse["coarse_bin"] = (coarse["start"] // coarse_size)

                merged = agg.merge(
                    coarse[["chromosome", "coarse_bin", stat]],
                    on=["chromosome", "coarse_bin"],
                    suffixes=("_fine_agg", "_coarse"),
                )
                if len(merged) > 1:
                    corr = merged[f"{stat}_fine_agg"].corr(merged[f"{stat}_coarse"])
                else:
                    corr = np.nan

                rows.append({
                    "accession": accession,
                    "statistic": stat,
                    "bin_size_a": fine_size,
                    "bin_size_b": coarse_size,
                    "cross_scale_correlation": corr,
                })

    result = pd.DataFrame(rows)
    result.to_csv(OUTPUT_FILE, index=False)

    print("Mean cross-scale correlation by statistic (higher = more scale-stable):")
    print(result.groupby("statistic")["cross_scale_correlation"].mean().sort_values(ascending=False))
    print(f"\nFull table saved to {OUTPUT_FILE}")
    print("\nPick the (statistic, bin_size) pair that is both stable across "
          "scales AND retains resolution you actually need — record the "
          "choice explicitly before starting Phase 1.")


if __name__ == "__main__":
    main()