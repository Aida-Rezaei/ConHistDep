"""
extract_signal_multiscale.py

Phase 0.3 - Multi-scale signal representation.

Extends your original extract_signal.py: instead of committing to
200 bp / mean-signal up front, compute several bin sizes and several
summary statistics per bin, so Phase 0 analysis can decide which
representation is stable before Phase 1 modeling locks one in.

Run this on chr1 only first (CHROM_SUBSET below) to compare
representations quickly, before scaling to the full genome.
"""

from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import pybigtools

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = [PROJECT_ROOT / "data" / "raw",
            PROJECT_ROOT / "data" / "raw" / "per_replicate"]
OUTPUT_DIR = PROJECT_ROOT / "results" / "multiscale_signal"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

BIN_SIZES = [200, 1000, 5000, 10000]

CHROM_SUBSET = ["chr1"]# for fast iteration. later to None for full genome but only after picking a representation

def bin_statistics(signal_values_np, background_threshold=0.0):
    """
    takes in a numpy array of signal values of a signal bin
    returns a dictionary with 4 keys
    """
    if len(signal_values_np) == 0: #if an empty array
        return {
            "mean": 0, "median": 0, "max": 0,
            "fraction_of_values_above_background": 0.0
        }
    return {
        "mean": float(signal_values_np.mean()),
        "median": float(np.median(signal_values_np)),
        "max": float(signal_values_np.max()),
        "fraction_of_values_above_background": float((signal_values_np > background_threshold).mean())
    }
def extract_multiscale(bigwig_file, bin_size, output_dir, chrom_subset=None):
    """
    one file at one bin size
    """
    output_file = OUTPUT_DIR / f"{bigwig_file.stem}_{bin_size}bp_signal.parquet"
    if output_file.exists():
        print(f"skipping {output_file} (exists)")
        return
    print(f"opening {bigwig_file.name} at {bin_size} bp resolution.")
    bw = pybigtools.open(str(bigwig_file))
    chromos = chrom_subset if chrom_subset else (
        [f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]]
    )
    rows = []
    bin_id = 0
    for chrom in chromos:
        if chrom not in bw.chroms():
            continue
        chrom_length = bw.chroms()[chrom]
        for start in range(0, chrom_length, bin_size):
            end = min(start+ bin_size, chrom_length)
            values = bw.values(chrom=chrom, start=start, end=end, fillna=0)
            stats = bin_statistics(values)
            rows.append({
                "bin_id": bin_id, "chromosome": chrom,
                "start": start, "end": end, **stats,
            })
            bin_id += 1
    df = pd.DataFrame(rows)
    df.to_parquet(output_file, index=False)
    print(f"saved {len(df)} bins to {output_file}")
    bw.close()

if __name__ == "__main__":
    bigwigs = []
    for data_dir in DATA_DIR:
        bigwigs.extend(sorted(data_dir.glob("*.bigWig")))
    print(f"Found {len(bigwigs)} bigWig files across {len(DATA_DIR)} directories.")
 
    jobs = [
        (bw_file, bin_size, OUTPUT_DIR, CHROM_SUBSET)
        for bw_file in bigwigs
        for bin_size in BIN_SIZES
    ]
    print(f"total (file, bin_size) jobs: {len(jobs)}")
    NUM_WORKERS = 3
    with ProcessPoolExecutor(max_workers=NUM_WORKERS) as executor:
        list(executor.map(
            extract_multiscale,
            [j[0] for j in jobs],
            [j[1] for j in jobs],
            [j[2] for j in jobs],
            [j[3] for j in jobs],
        ))
    print("DONE. Next step: run compare_representations.py")