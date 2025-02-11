import argparse
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


def process_file(file, output_directory):
    """Runs the Docker command for a single file."""
    smiles = str(file)
    short_inchikey = re.sub(r"^\d+_", "", file.stem)

    print(f"Processing file: {file}")
    command = [
        "docker", "run", "--rm=true",
        "-v", f"{os.getcwd()}:/cfmid/public/", "-i", "wishartlab/cfmid:latest",
        "sh", "-c", f"cd /cfmid/public/; cfm-predict '{smiles}' 0.001 "
                   f"/trained_models_cfmid4.0/[M+H]+/param_output.log "
                   f"/trained_models_cfmid4.0/[M+H]+/param_config.txt "
                   f"1 {output_directory}/{short_inchikey}.txt 0 0"
    ]
    
    subprocess.run(command)
    print(f"Processed {file} and output written to {output_directory}/{short_inchikey}.txt")


def run_docker_for_all_files(input_directory, output_directory, max_workers=4):
    """Runs Docker processing in parallel."""
    sorted_files = sorted(
        (file for file in Path(input_directory).rglob("*") if file.is_file()),
        key=lambda f: int(f.stem.split('_')[0])
    )

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_file, file, output_directory) for file in sorted_files]

        # Optionally wait for all to complete
        for future in futures:
            future.result()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process SMILES files using Docker in parallel.")
    parser.add_argument("--input-directory", type=str, default="in", help="Directory containing SMILES files")
    parser.add_argument("--output-directory", type=str, default="out", help="Directory for storing output files")
    parser.add_argument("--workers", type=int, default=4, help="Number of parallel workers")
    args = parser.parse_args()

    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)

    run_docker_for_all_files(args.input_directory, args.output_directory, args.workers)
