"""Installs the ENCODE SCREEN cCRE v4 resource directory"""
import os
import sys
import gzip
import shutil
import zipfile
from pathlib import Path
import yaml
import requests

# Globals
CONFIG = Path("config/config.yaml")

# ------------- #
# Functions     #
# ------------- #


def parse_config(config_path: str) -> dict:
    "Parses provided YAML configuration file as a dictionary."
    with open(config_path, "r", encoding="utf-8") as file:
        config = yaml.safe_load(file)
    return config


def download_file(url: str, destination: str) -> None:
    "Downloads a file from the specified URL to the destination path."
    filename = url.split("/")[-1]
    print(f"Downloading {filename}...")
    response = requests.get(url, stream=True, timeout=10)

    # Check for successful response
    if response.status_code == 200:
        with open(destination, "wb") as f:
            f.write(response.content)
        print(f"Downloaded {filename} to {destination}")
    else:
        print(f"Failed to download {filename}. Status code: {response.status_code}")


def compress_file(filepath: str) -> None:
    "Compresses the specified file using gzip."
    with open(filepath, "rb") as f_in:
        with gzip.open(filepath.with_suffix(".gz"), "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def main():
    """Main program"""
    # Parse config
    try:
        config_file = parse_config(CONFIG.resolve())
    except FileNotFoundError:
        print("Error: config not found.")

    # Register resource output directory and create if it doesn't exist
    op_dir = Path(config_file["FILEPATHS"]["RESOURCE_DIR"]).resolve()
    if op_dir.exists():
        print(f"Resource directory already exists: {op_dir}. Exiting.")
        sys.exit(1)
    else:
        print(f"Setting up resource directory: {op_dir}")
        op_dir.mkdir(parents=True, exist_ok=True)

        # Register encode screen urls
        data_urls = list(config_file["URLS"]["ENCODE"]["SCREEN"]["V4"].values())

        # Check if the directory is setup
        for url in data_urls:
            filename = url.split("/")[-1]
            filepath = op_dir / filename
            # Handle archive first
            if filename.endswith(".zip"):
                # Download as normal
                download_file(url, filepath)
                # Unzip the archive
                with zipfile.ZipFile(filepath, "r") as zip_ref:
                    zip_ref.extractall(op_dir)
                    # Remove the archive
                    os.remove(filepath)
                # Compress the extracted files
                for extracted_file in os.listdir(op_dir):
                    extracted_filepath = op_dir / extracted_file
                    if not extracted_filepath.name.endswith(
                        ".gz"
                    ) and not extracted_filepath.name.endswith(".md"):
                        compress_file(extracted_filepath)
                        # Remove the uncompressed file
                        os.remove(extracted_filepath)
            # Handle single files
            else:
                # Download the file
                download_file(url, filepath)
                # Compress the file
                if not filepath.name.endswith(".gz"):
                    compress_file(filepath)
                    # Remove the uncompressed file
                    os.remove(filepath)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
