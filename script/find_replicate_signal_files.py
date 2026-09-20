"""
find_replicate_signal_files.py

Fixes the root cause behind replicate_reproducibility.py finding
"only 1 file" everywhere: your existing metadata kept ENCODE's
POOLED signal track per experiment (biological_replicates: [1,2]),
not the separate per-replicate tracks. You can't measure replicate
reproducibility from a file that already has the replicates merged.

This queries each experiment you already have and pulls out any
files where biological_replicates has exactly ONE element — i.e.
genuinely replicate-specific signal tracks — if ENCODE released them
separately. Not every experiment will have these; some pipelines
only ever publish the pooled track, in which case that mark/cell
type genuinely can't be reproducibility-checked this way, and that
itself is worth noting in the paper rather than silently skipped.
"""

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import requests

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROVENANCE_FILE = PROJECT_ROOT / "results" / "provenance_table.csv"
OUTPUT_METADATA = PROJECT_ROOT / "results" / "replicate_signal_metadata.csv"
DOWNLOAD_DIR = PROJECT_ROOT / "data" / "raw" / "per_replicate"
DOWNLOAD_DIR.mkdir(parents=True, exist_ok=True)

BASE_URL = "https://www.encodeproject.org"
HEADERS = {"accept": "application/json"}


OUTPUT_TYPES = ("signal p-value", "fold change over control", "read-depth normalized signal")


def find_replicate_files(experiment_accession):
    url = f"{BASE_URL}/experiments/{experiment_accession}/?format=json"
    response = requests.get(url, headers=HEADERS, timeout=30)
    response.raise_for_status()
    exp = response.json()

    candidates = []
    for f in exp.get("files", []):
        if not isinstance(f, dict):
            continue
        if f.get("file_format") != "bigWig":
            continue
        if f.get("output_type") not in OUTPUT_TYPES:
            continue
        if f.get("status") != "released":
            continue
        bio_reps = f.get("biological_replicates", [])
        if len(bio_reps) == 1:  # single-replicate track, not pooled
            candidates.append({
                "experiment_accession": experiment_accession,
                "file_accession": f.get("accession"),
                "biological_replicate": bio_reps[0],
                "output_type": f.get("output_type"),
                "date_created": f.get("date_created", ""),
                "download_url": BASE_URL + f.get("href", ""),
            })

    # Multiple released files can exist for the same replicate when an
    # experiment was reprocessed under a newer pipeline version. Keep
    # the most recently created file per (replicate, output_type) and
    # report what got dropped, instead of either grabbing every
    # version (the original bug) or requiring preferred_default (which
    # excluded almost everything, since that flag mostly lives on
    # pooled files, not per-replicate ones).
    best_by_key = {}
    for c in candidates:
        key = (c["biological_replicate"], c["output_type"])
        if key not in best_by_key or c["date_created"] > best_by_key[key]["date_created"]:
            best_by_key[key] = c

    dropped = len(candidates) - len(best_by_key)
    if dropped > 0:
        print(f"  Kept {len(best_by_key)} file(s), dropped {dropped} older "
              f"duplicate(s) for {experiment_accession} (kept most recent "
              f"per replicate/output-type).")

    return list(best_by_key.values())


# Priority order when a replicate has more than one output type
# available: prefer signal p-value (most directly comparable across
# replicates), then fold-change, then raw read-depth normalized signal.
OUTPUT_TYPE_PRIORITY = ["signal p-value", "fold change over control", "read-depth normalized signal"]


def deduplicate_and_align(hits):
    """Given all per-replicate hits for one experiment (possibly mixing
    output types across replicates), pick ONE output type and use it
    for every replicate that has it — so replicate correlations always
    compare like with like. Returns the aligned subset of hits."""
    if not hits:
        return []

    replicates = sorted(set(h["biological_replicate"] for h in hits))
    by_type = {}
    for h in hits:
        by_type.setdefault(h["output_type"], {})[h["biological_replicate"]] = h

    # Pick the highest-priority output type that covers the most
    # replicates; ties broken by priority order.
    best_type, best_coverage = None, -1
    for output_type in OUTPUT_TYPE_PRIORITY:
        reps_covered = len(by_type.get(output_type, {}))
        if reps_covered > best_coverage:
            best_type, best_coverage = output_type, reps_covered

    if best_type is None or best_coverage < 2:
        return []  # no single output type covers >=2 replicates

    aligned = list(by_type[best_type].values())
    dropped_types = set(by_type.keys()) - {best_type}
    if dropped_types:
        print(f"  Using '{best_type}' for {best_coverage}/{len(replicates)} "
              f"replicates; ignoring other output type(s) found: {dropped_types}")

    return aligned


def main():
    provenance = pd.read_csv(PROVENANCE_FILE)
    exp_col = "Experiment_Accession" if "Experiment_Accession" in provenance.columns else "experiment_accession"
    cell_col = "Cell_Type" if "Cell_Type" in provenance.columns else "cell_type"
    mark_col = "Histone_Mark" if "Histone_Mark" in provenance.columns else "histone_mark"

    lookup = provenance.set_index(exp_col)[[cell_col, mark_col]].to_dict("index")

    rows = []
    experiments = provenance[exp_col].unique()
    print(f"Fetching metadata for {len(experiments)} experiments (parallel)...")

    def process_experiment(exp_acc):
        try:
            hits = find_replicate_files(exp_acc)
        except Exception as e:
            return exp_acc, None, str(e)
        hits = deduplicate_and_align(hits)
        return exp_acc, hits, None

    done_count = 0
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {executor.submit(process_experiment, acc): acc for acc in experiments}
        for future in as_completed(futures):
            exp_acc, hits, error = future.result()
            done_count += 1
            print(f"[{done_count}/{len(experiments)}] {exp_acc}", end="  ")

            if error:
                print(f"FAILED: {error}")
                continue

            if len(hits) < 2:
                print(f"only {len(hits)} usable track(s) "
                      f"({lookup[exp_acc][cell_col]} {lookup[exp_acc][mark_col]}) — skipping.")
            else:
                print(f"OK — {len(hits)} tracks "
                      f"({lookup[exp_acc][cell_col]} {lookup[exp_acc][mark_col]})")

            for hit in hits:
                hit["cell_type"] = lookup[exp_acc][cell_col]
                hit["histone_mark"] = lookup[exp_acc][mark_col]
                rows.append(hit)

    df = pd.DataFrame(rows)
    df.to_csv(OUTPUT_METADATA, index=False)
    print(f"\nSaved {len(df)} per-replicate file records to {OUTPUT_METADATA}")

    checkable = df.groupby(["cell_type", "histone_mark"]).size()
    checkable = checkable[checkable >= 2]
    print(f"\n{len(checkable)} / {provenance[[cell_col, mark_col]].drop_duplicates().shape[0]} "
          f"mark/cell-type combinations have >=2 per-replicate tracks and can be checked:")
    print(checkable.to_string())

    # Sanity check: after deduplicate_and_align, each (experiment,
    # replicate) should map to exactly ONE file. If this still fires,
    # something in ENCODE's metadata doesn't fit the assumptions here
    # (e.g. two files with the same output_type and same date_created)
    # — inspect manually rather than assuming it's safe to merge.
    dup_check = df.groupby(["experiment_accession", "biological_replicate"]).size()
    duplicates = dup_check[dup_check > 1]
    if len(duplicates):
        print(f"\nWARNING: {len(duplicates)} (experiment, replicate) pairs still "
              f"map to multiple files — inspect these manually before "
              f"downloading, do not assume they're safe to merge:")
        print(duplicates.to_string())
        return

    # Download them — in parallel, with visible per-file progress, so
    # this doesn't look stalled on a large bigWig for minutes at a time.
    to_download = [row for _, row in df.iterrows()
                   if not (DOWNLOAD_DIR / f"{row['file_accession']}.bigWig").exists()]
    already_have = len(df) - len(to_download)
    if already_have:
        print(f"\n{already_have} file(s) already downloaded, skipping.")
    print(f"Downloading {len(to_download)} file(s), 4 at a time...")

    def download_one(row):
        output_file = DOWNLOAD_DIR / f"{row['file_accession']}.bigWig"
        tmp_file = output_file.with_suffix(".bigWig.part")
        try:
            r = requests.get(row["download_url"], stream=True, timeout=120)
            r.raise_for_status()
            total = int(r.headers.get("content-length", 0))
            downloaded = 0
            with open(tmp_file, "wb") as fh:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        fh.write(chunk)
                        downloaded += len(chunk)
            tmp_file.rename(output_file)
            size_mb = downloaded / (1024 * 1024)
            return row["file_accession"], True, f"{size_mb:.0f} MB"
        except Exception as e:
            if tmp_file.exists():
                tmp_file.unlink()
            return row["file_accession"], False, str(e)

    done = 0
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = {executor.submit(download_one, row): row for row in to_download}
        for future in as_completed(futures):
            accession, success, info = future.result()
            done += 1
            status = "OK" if success else "FAILED"
            print(f"[{done}/{len(to_download)}] {accession}: {status} ({info})")


if __name__ == "__main__":
    main()