"""
1. Build a full experiment-level provenance table.
For every bigWig file already listed in results/encode_metadata.csv,
query the ENCODE REST API for the file and its parent experiment to
record everything needed to separate biological variation from
experiment/batch variation later.
Run this BEFORE any signal extraction. Every downstream file should
be traceable back to a row in this table.
"""
import sys
import time
from pathlib import Path
import pandas as pd
import requests

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
METADATA_FILE = PROJECT_ROOT / "results" / "encode_metadata.csv"
OUTPUT_FILE = PROJECT_ROOT / "results" / "provenance_table.csv"
BASE_URL = "https://www.encodeproject.org"
HEADERS = {
    "accept": "application/json",
}
def fetch_json(path):
    """
    GET a single ENCODE object and return its JSON.
    Includes basic retry handling.
    """
    if not path:
        return {}
    if path.startswith("http://") or path.startswith("https://"):
        url = path
    else:
        url = f"{BASE_URL}{path}"
    separator = "&" if "?" in url else "?"
    url = f"{url}{separator}format=json"
    last_error = None
    for attempt in range(3):
        try:
            response = requests.get(
                url,
                headers=HEADERS,
                timeout=30,
            )
            if response.status_code == 200:
                return response.json()
            last_error = (
                f"HTTP {response.status_code}: "
                f"{response.text[:300]}"
            )
        except requests.RequestException as exc:
            last_error = str(exc)
        if attempt < 2:
            time.sleep(2)
    raise RuntimeError(
        f"Failed to fetch ENCODE object {path}: {last_error}"
    )
def resolve_encode_object(value):
    """
    Normalize an ENCODE field that may be:
        - None
        - a string accession
        - an ENCODE relative URL
        - a full URL
        - an already-expanded dictionary
    Returns a dictionary when possible.
    """
    if value is None:
        return {}
    if isinstance(value, dict):
        return value
    if isinstance(value, str):
        if value.startswith("/") or value.startswith("http"):
            try:
                return fetch_json(value)
            except Exception:
                return {}
    return {}
def extract_identifier(value):
    """
    Convert an ENCODE object/reference into a useful string identifier.
    Preference: accession -> @id -> uuid -> name/title
    """
    if value is None:
        return None
    if isinstance(value, dict):
        return (
            value.get("accession")
            or value.get("@id")
            or value.get("uuid")
            or value.get("name")
            or value.get("title")
        )
    return str(value)
def extract_title(value):
    """
    Extract a human-readable title/name from an ENCODE object.
    """
    if value is None:
        return None
    if isinstance(value, dict):
        return (
            value.get("title")
            or value.get("name")
            or value.get("label")
        )
    return str(value)
def extract_donors(experiment_obj):
    """
    Extract donor identifiers from experiment replicates.
    ENCODE may represent donor as:
        - a dictionary
        - a URL/reference
        - an accession string
    """
    donors = set()
    replicates = experiment_obj.get("replicates", [])
    if not isinstance(replicates, list):
        return None
    for replicate in replicates:
        if not isinstance(replicate, dict):
            continue
        library = replicate.get("library", {})
        library = resolve_encode_object(library)
        if not library:
            continue
        biosample = library.get("biosample", {})
        biosample = resolve_encode_object(biosample)
        if not biosample:
            continue
        donor = biosample.get("donor")
        donor_obj = resolve_encode_object(donor)
        if donor_obj:
            donor_id = extract_identifier(donor_obj)
        else:
            donor_id = extract_identifier(donor)
        if donor_id:
            donors.add(str(donor_id))
    if not donors:
        return None
    return ";".join(sorted(donors))
def extract_cell_types(experiment_obj):
    """
    Extract cell type information from experiment replicates.
    """
    cell_types = set()
    replicates = experiment_obj.get("replicates", [])
    if not isinstance(replicates, list):
        return None
    for replicate in replicates:
        if not isinstance(replicate, dict):
            continue
        library = resolve_encode_object(
            replicate.get("library")
        )
        if not library:
            continue
        biosample = resolve_encode_object(
            library.get("biosample")
        )
        if not biosample:
            continue
        cell_type = biosample.get("biosample_ontology")
        cell_type_obj = resolve_encode_object(cell_type)
        if cell_type_obj:
            value = (
                cell_type_obj.get("term_name")
                or cell_type_obj.get("preferred_term")
                or cell_type_obj.get("name")
                or cell_type_obj.get("label")
            )
        else:
            value = cell_type
        if value:
            cell_types.add(str(value))
    if not cell_types:
        return None
    return ";".join(sorted(cell_types))
def extract_control_accessions(experiment_obj):
    """
    Extract control/input experiment or file references where present.
    """
    controls = set()
    # Possible experiment-level control fields.
    possible_fields = [
        "control",
        "controls",
        "control_experiments",
        "control_files",
    ]
    for field in possible_fields:
        value = experiment_obj.get(field)
        if value is None:
            continue
        if not isinstance(value, list):
            value = [value]
        for item in value:
            identifier = extract_identifier(item)
            if identifier:
                controls.add(str(identifier))
    if not controls:
        return None
    return ";".join(sorted(controls))
def extract_qc_metrics(file_obj):
    """
    Extract NRF, PBC1 and PBC2 from ENCODE quality metrics.
    Quality metrics may be represented as:
        - dictionaries embedded in the file object
        - references to quality metric objects
        - lists of either
    """
    nrf = None
    pbc1 = None
    pbc2 = None
    qc_metrics = file_obj.get("quality_metrics", [])
    if not isinstance(qc_metrics, list):
        qc_metrics = [qc_metrics]
    for qc in qc_metrics:
        qc_obj = resolve_encode_object(qc)
        if not qc_obj:
            continue
        nrf = qc_obj.get("NRF", nrf)
        pbc1 = qc_obj.get("PBC1", pbc1)
        pbc2 = qc_obj.get("PBC2", pbc2)
        nrf = qc_obj.get("nrf", nrf)
        pbc1 = qc_obj.get("pbc1", pbc1)
        pbc2 = qc_obj.get("pbc2", pbc2)
    return nrf, pbc1, pbc2
def extract_pipeline_version(file_obj):
    """
    Extract the most useful pipeline/analysis identifier available.
    ENCODE analysis schemas vary, so inspect several possible fields.
    """
    analyses = file_obj.get("analyses", [])
    if not isinstance(analyses, list):
        return None
    values = []
    for analysis in analyses:
        analysis_obj = resolve_encode_object(analysis)
        if not analysis_obj:
            continue
        candidates = [
            analysis_obj.get("pipeline_version"),
            analysis_obj.get("pipeline"),
            analysis_obj.get("pipeline_award_rfas"),
            analysis_obj.get("software_versions"),
            analysis_obj.get("title"),
        ]
        for value in candidates:
            if value is None:
                continue
            if isinstance(value, dict):
                identifier = extract_identifier(value)
            else:
                identifier = str(value)
            if identifier:
                values.append(identifier)
    if not values:
        return None
    unique_values = list(dict.fromkeys(values))
    return ";".join(unique_values)
def extract_file_provenance(file_accession):
    """
    Pull provenance fields for one ENCODE file accession.
    """
    file_obj = fetch_json(
        f"/files/{file_accession}/"
    )
    experiment_path = file_obj.get("dataset")
    experiment_obj = resolve_encode_object(
        experiment_path
    )
    file_status = file_obj.get("status")
    file_format = file_obj.get("file_format")
    assembly = (
        file_obj.get("assembly")
        or file_obj.get("genome_assembly")
    )
    output_type = file_obj.get("output_type")
    md5sum = file_obj.get("md5sum")
    experiment_accession = (
        experiment_obj.get("accession")
        or extract_identifier(experiment_path)
    )
    assay = (
        experiment_obj.get("assay_term_name")
        or experiment_obj.get("assay_title")
    )
    target = experiment_obj.get("target")
    target_obj = resolve_encode_object(target)
    if target_obj:
        target_label = (
            target_obj.get("label")
            or target_obj.get("name")
            or target_obj.get("title")
        )
    else:
        target_label = extract_title(target)
    lab = experiment_obj.get("lab")
    lab_obj = resolve_encode_object(lab)
    if lab_obj:
        lab_name = (
            lab_obj.get("title")
            or lab_obj.get("name")
        )
    else:
        lab_name = extract_title(lab)
    biological_replicates = file_obj.get(
        "biological_replicates",
        []
    )
    if isinstance(biological_replicates, list):
        biological_replicates_value = ";".join(
            str(rep)
            for rep in biological_replicates
        )
    elif biological_replicates:
        biological_replicates_value = str(
            biological_replicates
        )
    else:
        biological_replicates_value = None
    technical_replicates = file_obj.get(
        "technical_replicates",
        []
    )
    if isinstance(technical_replicates, list):
        technical_replicate_value = ";".join(
            str(rep)
            for rep in technical_replicates
        )
        technical_replicate_count = len(
            technical_replicates
        )
    elif technical_replicates:
        technical_replicate_value = str(
            technical_replicates
        )
        technical_replicate_count = 1
    else:
        technical_replicate_value = None
        technical_replicate_count = 0
    donor = extract_donors(experiment_obj)
    cell_type = extract_cell_types(
        experiment_obj
    )
    control_accession = extract_control_accessions(
        experiment_obj
    )
    nrf, pbc1, pbc2 = extract_qc_metrics(
        file_obj
    )
    pipeline_version = extract_pipeline_version(
        file_obj
    )
    return {
        "file_accession": file_accession,
        "experiment_accession": experiment_accession,
        "cell_type_api": cell_type,
        "assay": assay,
        "target_label": target_label,
        "lab": lab_name,
        "donor": donor,
        "biological_replicates": (
            biological_replicates_value
        ),
        "technical_replicates": (
            technical_replicate_value
        ),
        "technical_replicate_count": (
            technical_replicate_count
        ),
        "control_input_accession": (
            control_accession
        ),
        "output_type": output_type,
        "signal_type": output_type,
        "file_format": file_format,
        "assembly": assembly,
        "pipeline_version": pipeline_version,
        "md5sum": md5sum,
        "status": file_status,
        "nrf": nrf,
        "pbc1": pbc1,
        "pbc2": pbc2,
    }
def main():
    if not METADATA_FILE.exists():
        raise FileNotFoundError(
            f"Metadata file not found:\n{METADATA_FILE}"
        )
    df = pd.read_csv(
        METADATA_FILE
    )
    if "File_Accession" in df.columns:
        accession_col = "File_Accession"
    elif "file_accession" in df.columns:
        accession_col = "file_accession"
    else:
        raise ValueError(
            "Could not find a file accession column. "
            "Expected either 'File_Accession' or "
            "'file_accession'.\n"
            f"Available columns: {list(df.columns)}"
        )
    accessions = (
        df[accession_col]
        .dropna()
        .astype(str)
        .str.strip()
    )
    accessions = accessions[
        accessions != ""
    ].unique()
    total = len(accessions)
    print(
        f"Found {total} unique file accessions."
    )
    rows = []
    for i, accession in enumerate(
        accessions,
        start=1,
    ):
        print(
            f"[{i}/{total}] {accession}"
        )
        try:
            provenance = extract_file_provenance(
                accession
            )
            rows.append(
                provenance
            )
            print(
                "  OK"
            )
        except Exception as exc:
            error_message = (
                f"{type(exc).__name__}: {exc}"
            )
            print(
                f"  FAILED: {error_message}"
            )
            # IMPORTANT:
            # Always preserve the expected columns even
            # when an individual API request fails.
            rows.append(
                {
                    "file_accession": accession,
                    "status": f"FAILED: {error_message}",
                    "experiment_accession": None,
                    "cell_type_api": None,
                    "assay": None,
                    "target_label": None,
                    "lab": None,
                    "donor": None,
                    "biological_replicates": None,
                    "technical_replicates": None,
                    "technical_replicate_count": None,
                    "control_input_accession": None,
                    "output_type": None,
                    "signal_type": None,
                    "file_format": None,
                    "assembly": None,
                    "pipeline_version": None,
                    "md5sum": None,
                    "nrf": None,
                    "pbc1": None,
                    "pbc2": None,
                }
            )

        # Be polite to the ENCODE API.
        time.sleep(0.5)
    provenance_df = pd.DataFrame(
        rows
    )
    merged = df.merge(
        provenance_df,
        left_on=accession_col,
        right_on="file_accession",
        how="left",
        suffixes=("", "_api"),
    )
    OUTPUT_FILE.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    merged.to_csv(
        OUTPUT_FILE,
        index=False,
    )
    print()
    print(
        f"Saved provenance table with "
        f"{len(merged)} rows to:"
    )
    print(
        OUTPUT_FILE
    )
    required_columns = [
        "file_accession",
        "experiment_accession",
        "donor",
        "biological_replicates",
        "technical_replicates",
    ]
    missing_columns = [
        col
        for col in required_columns
        if col not in merged.columns
    ]
    if missing_columns:
        print()
        print(
            "WARNING: Expected columns are missing:"
        )
        for col in missing_columns:
            print(
                f"  - {col}"
            )
    else:
        missing = merged[
            merged["donor"].isna()
            | merged["biological_replicates"].isna()
        ]
        if len(missing):
            print()
            print(
                f"WARNING: {len(missing)} files are "
                "missing donor and/or biological "
                "replicate information:"
            )
            print(
                missing[
                    [accession_col]
                ].to_string(
                    index=False
                )
            )

        else:

            print()
            print(
                "Sanity check: all files have "
                "donor and biological replicate "
                "information."
            )
    if "status" in merged.columns:
        failed_mask = (
            merged["status"]
            .astype(str)
            .str.startswith("FAILED:")
        )
        failed = merged[
            failed_mask
        ]
        if len(failed):
            print()
            print(
                f"WARNING: {len(failed)} "
                "files failed ENCODE provenance "
                "retrieval."
            )
            print(
                failed[
                    [accession_col, "status"]
                ].to_string(
                    index=False
                )
            )
        else:
            print()
            print(
                "API check: all file requests "
                "completed successfully."
            )
if __name__ == "__main__":
    main()