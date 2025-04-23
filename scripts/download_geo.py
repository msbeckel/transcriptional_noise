#!/usr/bin/env python3
"""
download_geo.py

Downloads a GEO series and its supplementary files to a specified folder.
"""

import os
import argparse
import GEOparse  # GEOparse.get_GEO fetches GEO series metadata and files :contentReference[oaicite:0]{index=0}

def download_geo_data(geo_id: str, out_dir: str) -> str:
    """
    Download GEO series data.

    Args:
        geo_id: GEO accession (e.g., 'GSE141044').
        out_dir: Directory to store downloaded files.

    Returns:
        Path to the folder containing the raw matrix files.
    """
    os.makedirs(out_dir, exist_ok=True)
    print(f"Downloading GEO series {geo_id} into '{out_dir}'…")
    gse = GEOparse.get_GEO(geo=geo_id, destdir=out_dir, include_data=False)
    gse.download_supplementary_files()
    print("Download complete.")

    # Locate the Matrix Market files
    for root, dirs, files in os.walk(os.path.join(out_dir, geo_id)):
        if "matrix.mtx" in files:
            print(f"Found matrix directory: {root}")
            return root

    raise FileNotFoundError("Could not find 'matrix.mtx' under downloaded GEO files.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download GEO dataset")
    parser.add_argument("--geo_id",  required=True, help="GEO accession ID (e.g., GSE141044)")
    parser.add_argument("--out_dir", required=True, help="Output directory")
    args = parser.parse_args()

    matrix_dir = download_geo_data(args.geo_id, args.out_dir)
    print(f"Matrix files are in: {matrix_dir}")
