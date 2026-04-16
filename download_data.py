"""Download MetaCyc data files from Google Drive."""

from __future__ import annotations

import os
from pathlib import Path

import gdown

DATA_DIR = Path(__file__).parent / "data"


def download_file(file_id: str, output_path: str) -> None:
    """Download a file from Google Drive by its file ID."""
    gdown.download(id=file_id, output=output_path, quiet=False)


FILES_INFO = {
    "good_chems.txt":             "165xAHsCx1bD5PgSCFYbQon67YXlJWVQ9",
    "good_reactions.txt":         "1O7KMm8pEmyu6n8jkWAnqFzzp5H9YZf88",
    "minimal_metabolites.txt":    "17rmvSCeBm0ZRdfFmnOGwhLluF7rIz6Dz",
    "ubiquitous_metabolites.txt": "1n2YZcBSJI9hemkwZNrOO4WiJC4IbXB_f",
}


def main() -> None:
    DATA_DIR.mkdir(exist_ok=True)
    for filename, file_id in FILES_INFO.items():
        output_path = str(DATA_DIR / filename)
        if not os.path.exists(output_path):
            print(f"Downloading {filename}...")
            download_file(file_id, output_path)
        else:
            print(f"Already exists: {filename}")

    print("\nAvailable files:")
    for filename in sorted(FILES_INFO):
        print(f"  {DATA_DIR / filename}")


if __name__ == "__main__":
    main()
