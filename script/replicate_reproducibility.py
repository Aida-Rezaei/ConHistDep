"""

biological question: comparing the siganl tracks
from biological replictes of the same cell type and 
histone marker, using Pearson correlation across genomic bins.
"""

from pathlib import Path
from itertools import combinations

import pandas as pd
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parent.parent
REPLICATE_METADATA_FILE = PROJECT_ROOT / "results" / "replicate_signal_metadata.csv"
SIGNAL_DIR = PROJECT_ROOT / "results" / "multiscale_signal"
RESULTS_DIR = PROJECT_ROOT / "results" / "replicate_reproducibility"
FIGURE_DIR = PROJECT_ROOT / "figures" / "replicate_reproducibility"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

PRIMARY_BIN_SIZE = 10000 #10kb size: dividing the genome into bins instead of comparing individual bases
PRIMARY_STAT = "mean" #genomic bins are represented by the mean signal

#write the file-loading function
def load_signal(accession):
    file = SIGNAL_DIR / f"{accession}_{PRIMARY_BIN_SIZE}bp_signal.parquet"
    if not file.exists():
        return None
    return pd.read_parquet(file)[["chromosome", "start", "end", PRIMARY_STAT]]

#load the metadata
metadata = pd.read_csv(REPLICATE_METADATA_FILE)
#initialize an empty list to collect one dict per replicate pair
rows = []
grouped_metadata = metadata.groupby(by = ["cell_type", "histone_mark"])
for (cell, mark), group in grouped_metadata:
    accessions = group["file_accession"].unique().tolist()
    if len(accessions) < 2:
        print(f"{cell}|{mark}: only N per-replicate track(s) - not checkable: skipping!")
        continue
    signals = {}
    for acc in accessions:
        sig = load_signal(acc)
        if sig is None:
            print("missing extracted signal for {acc}. run ")
            continue
        signals[acc] = sig
        for acc_a, acc_b in combinations(signals.keys(), 2):
            merged = pd.merge(signals[acc_a], 
                              signals[acc_b], 
                              how="inner", 
                              on=["chromosome", "start", "end"],
                              suffixes=("_a", "_b"))
            if merged.empty:
                print("empty merged df for {acc_a} vs {acc_b}, skipping.")
                continue
            corraltion = merged[f"{PRIMARY_STAT}_a"].corr(merged[f"{PRIMARY_STAT}_b"])
            rows.append({
                "cell_type": cell, "histone_mark": mark,
                "replicate_a": acc_a, "replicate_b": acc_b,
                "pearson_r": corraltion, "n_bins": len(merged)
            })
            plt.figure(figsize=(5, 5))
            plt.scatter(x= merged[f"{PRIMARY_STAT}_a"],
                        y = merged[f"{PRIMARY_STAT}_b"],
                        s=1, alpha=0.3)
            plt.title(f"{cell} {mark}: r = {corraltion:.3f}")
            plt.xlabel(f"{acc_a} signal")
            plt.ylabel(f"{acc_b} signal")
            plt.tight_layout()
            plt.savefig(FIGURE_DIR / f"{cell}_{mark}_{acc_a}_vs_{acc_b}.png")
            plt.close()
report = pd.DataFrame(rows)
report.to_csv(path_or_buf=RESULTS_DIR / "replicate_correlations.csv", index=False)
if report.empty:
    print("nothing checkable")
summary = report.groupby(["cell_type","histone_mark"])["pearson_r"].mean().sort_values()
print("Summary", summary)
low_quality = report[report["pearson_r"] < 0.7]
if len(low_quality) != 0:
    print(f"WARNING: {len(low_quality)} replicate pairs have r < 0.7:")
    print(low_quality[["cell_type", "histone_mark", "pearson_r"]].to_string(index=False))