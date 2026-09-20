import sys
import time
from pathlib import Path
 
import pandas as pd
import requests
 
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
from config.config import CELL_TYPES

BASE_URL = "https://www.encodeproject.org/search/"
FILES_BASE = "https://www.encodeproject.org"
HEADERS = {"accept": "application/json"}
ACCESSIBILITY_ASSAYS = ["DNase-seq", "ATAC-seq"]
OUTPUT_METADATA = PROJECT_ROOT / "results" / "accessibility_metadata.csv"
DOWNLOAD_DIR = PROJECT_ROOT / "data" / "raw" / "accessibility"
DOWNLOAD_DIR.mkdir(parents=True, exist_ok=True)

def find_accessibility_experiments()
    results = []
    for cell in CELL_TYPES:
        for assay in ACCESSIBILITY_ASSAYS:
            print(f"searching {cell}|{assay}")
            params = {
                "type": "experiment",
                "assay_title": assay,
                "biosample_ontology.term_name": cell,
                "status": "released",
                "assembly": "GRCh38",
                "format": "json",
                "limit": "all"
            }
            response = requests.get(
                BASE_URL,
                params=params,
                headers=HEADERS
            )
            response.raise_for_status()
            data = response.json()
            print(f"found {len(data['@graph'])} experiments")