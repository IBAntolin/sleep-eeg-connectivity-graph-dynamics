"""
ANPHY-Sleep Download Utilities
Compatible with current OSF file naming (EPCTLxx.zip or EPCTLxx-2025.zip)
"""

import re
import zipfile
from pathlib import Path
from typing import List, Optional

import requests

OSF_PROJECT_ID = "r26fh"
OSF_API_BASE = "https://api.osf.io/v2"

SKIP_SUBJECT_NUMBERS = {8}

# Flexible pattern to match both EPCTL29.zip and EPCTL29-2025.zip
SUBJECT_PATTERN = re.compile(r"EPCTL(\d{2})(?:-\d{4})?\.zip", re.IGNORECASE)


def list_osf_files(project_id: str = OSF_PROJECT_ID) -> list[dict]:
    """List all files in the OSF project (handles pagination)."""
    files = []
    url = f"{OSF_API_BASE}/nodes/{project_id}/files/osfstorage/?page[size]=100"
    while url:
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        payload = resp.json()
        for entry in payload.get("data", []):
            attrs = entry["attributes"]
            files.append({
                "name": attrs["name"],
                "size_bytes": attrs.get("size"),
                "download_url": entry["links"]["download"],
            })
        url = (payload.get("links") or {}).get("next")
    return files


def _stream_download(url: str, dest_path: Path):
    """Download file with progress."""
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=90) as r:
        r.raise_for_status()
        total = int(r.headers.get("content-length", 0))
        downloaded = 0
        with open(dest_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024*1024):
                f.write(chunk)
                downloaded += len(chunk)
                if total:
                    print(f"\r  {dest_path.name}: {100 * downloaded / total:5.1f}%", end="")
        print()


def download_shared_files(dest_dir: str | Path):
    """Download artifact matrix, .pos file, and subject details (once)."""
    dest_dir = Path(dest_dir)
    files = list_osf_files()

    wanted = {
        "artifact matrix.zip": "artifact_matrix.zip",
        "details information for healthy subjects.csv": "details_information_for_healthy_subjects.csv",
        "co-registered average positions.pos": "co-registered_average_positions.pos",
    }

    for f in files:
        key = f["name"].strip().lower()
        dest_name = wanted.get(key)
        if dest_name is None:
            continue

        dest_path = dest_dir / dest_name
        if dest_path.exists():
            print(f"  {dest_name} already exists.")
            continue

        print(f"Downloading shared file: {f['name']}")
        _stream_download(f["download_url"], dest_path)

    # Extract artifact matrix
    zip_path = dest_dir / "artifact_matrix.zip"
    extract_dir = dest_dir / "artifact_matrix"
    if zip_path.exists() and not extract_dir.exists():
        print("Extracting artifact matrices...")
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(extract_dir)


def download_subject(subject_num: int, dest_dir: str | Path, delete_zip: bool = True) -> Optional[Path]:
    """Download and extract one subject. Returns path to subject folder."""
    if subject_num in SKIP_SUBJECT_NUMBERS:
        print(f"Subject {subject_num:02d} is skipped (not available).")
        return None

    dest_dir = Path(dest_dir)
    files = list_osf_files()

    subject_file = None
    for f in files:
        if SUBJECT_PATTERN.search(f["name"]) and f"epctl{subject_num:02d}" in f["name"].lower():
            subject_file = f
            break

    if not subject_file:
        print(f"Could not find EPCTL{subject_num:02d} on OSF.")
        return None

    subject_id = f"EPCTL{subject_num:02d}"
    subject_dir = dest_dir / subject_id

    # Skip if already extracted
    if subject_dir.exists() and any(subject_dir.glob("**/*.edf")):
        print(f"{subject_id} already downloaded and extracted.")
        return subject_dir

    zip_path = dest_dir / subject_file["name"]
    print(f"Downloading {subject_file['name']} ({subject_file['size_bytes']/1e9:.2f} GB)...")
    _stream_download(subject_file["download_url"], zip_path)

    print(f"Extracting {subject_id}...")
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(subject_dir)

    if delete_zip:
        zip_path.unlink()

    return subject_dir


def get_subject_edf_path(subject_num: int, dest_dir: str | Path) -> Optional[Path]:
    """Return the main .edf file path for a subject."""
    dest_dir = Path(dest_dir)
    subject_dir = dest_dir / f"EPCTL{subject_num:02d}"
    
    if not subject_dir.exists():
        return None
    
    edf_files = list(subject_dir.glob("**/*.edf"))
    return edf_files[0] if edf_files else None


def get_subject_scoring_path(subject_num: int, dest_dir: str | Path) -> Optional[Path]:
    """
    Return the sleep scoring .txt file path for a subject. Filename inside
    each subject's zip isn't confirmed to follow one fixed pattern, so this
    takes any .txt file in the subject folder that isn't obviously
    something else -- if a subject folder ever contains more than one .txt
    file, this will need tightening (run discover_subject_files from
    preprocessing.py on that folder to see what's actually in there).
    """
    dest_dir = Path(dest_dir)
    subject_dir = dest_dir / f"EPCTL{subject_num:02d}"

    if not subject_dir.exists():
        return None

    txt_files = list(subject_dir.glob("**/*.txt"))
    return txt_files[0] if txt_files else None


def get_subject_artifact_path(subject_num: int, dest_dir: str | Path) -> Optional[Path]:
    """
    Return the *_artndxn.mat artifact matrix path for a subject, from the
    shared artifact_matrix/ folder extracted by download_shared_files()
    (one .mat per subject, all subjects in one zip -- not per-subject
    folders like the EDF/scoring files).
    """
    dest_dir = Path(dest_dir)
    artifact_dir = dest_dir / "artifact_matrix"

    if not artifact_dir.exists():
        return None

    matches = list(artifact_dir.glob(f"**/EPCTL{subject_num:02d}*artndxn*.mat"))
    return matches[0] if matches else None


def get_available_subjects() -> List[int]:
    """Return list of all available subject numbers."""
    return [i for i in range(1, 30) if i not in SKIP_SUBJECT_NUMBERS]