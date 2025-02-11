import aiofiles
import argparse
import asyncio
import polars
import os
import time
import shutil
from rdkit import Chem
from rdkit.Chem import Descriptors
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor


def smiles_to_mol_data(smiles):
    """Process SMILES to get both inchikey and molecular weight in one RDKit call."""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        inchikey = Chem.MolToInchiKey(mol).split("-")[0]
        weight = Descriptors.MolWt(mol)
        return inchikey, weight
    return None, float("inf")


async def cache_fragmented_files(search_directory: Path, skip_unfragmented: bool = True) -> set:
    """Pre-cache all fragmented files."""
    fragmented = set()
    if not search_directory.exists():
        return fragmented

    for file in search_directory.rglob("*"):
        if not file.is_file():
            continue
        
        if skip_unfragmented:
            fragmented.add(file.stem)
        else:
            try:
                async with aiofiles.open(file, "r", encoding="utf-8") as f:
                    lines = await f.readlines()
                    if len(lines) >= 3:
                        third_line = lines[2].strip()
                        if "#ID=" in third_line:
                            inchikey = third_line.split("#ID=")[1].strip()
                            fragmented.add(inchikey)
                    else:
                        print(f"File {file} has less than 3 lines.")
            except UnicodeDecodeError as e:
                print(f"UnicodeDecodeError reading file {file}: {e}")
            except OSError as e:
                print(f"OSError reading file {file}: {e}")
            except Exception as e:
                print(f"Error reading file {file}: {e}")

    return fragmented


def process_smiles_file(input_file, output_directory, search_directory, skip_unfragmented=True):
    start_time = time.time()

    print("Loading and processing SMILES file...")
    df_initial = polars.read_csv(input_file, separator="\t")
    total_initial = df_initial.height
    print(f"Loaded {total_initial} initial entries")
    
    # Remove entries with dots in SMILES
    df_no_dot = df_initial.filter(~df_initial["smiles"].str.contains("\\."))
    total_no_dot = df_no_dot.height
    print(f"Removed {total_initial - total_no_dot} entries with dots in SMILES")

    # Process SMILES in parallel
    print("Processing SMILES structures...")
    with ThreadPoolExecutor() as executor:
        results = list(executor.map(smiles_to_mol_data, df_no_dot["smiles"]))
        inchikeys, weights = zip(*results)
    
    # Add computed columns
    df_processed = df_no_dot.with_columns([
        polars.Series("short_inchikey", inchikeys),
        polars.Series("molecular_weight", weights)
    ])

    # Remove invalid entries
    df_valid = df_processed.filter(polars.col("short_inchikey").is_not_null() & (polars.col("short_inchikey") != ""))
    total_no_invalid = df_valid.height
    print(f"Removed {total_no_dot - total_no_invalid} invalid structures")

    # Sort by molecular weight and remove duplicates
    df_sorted = df_valid.sort("molecular_weight")
    df_unique = df_sorted.unique(subset=["short_inchikey"], keep="first", maintain_order=True)
    total_no_dup = df_unique.height
    print(f"Removed {total_no_invalid - total_no_dup} duplicate structures")

    # Cache fragmented files
    print("Checking fragmented files...")
    fragmented_files = asyncio.run(cache_fragmented_files(search_directory=Path(search_directory)))
    print(f"Found {len(fragmented_files)} fragmented files")

    # Split into non-fragmented and fragmented entries
    df_non_fragmented = df_unique.filter(~polars.col("short_inchikey").is_in(fragmented_files))
    df_fragmented = df_unique.filter(polars.col("short_inchikey").is_in(fragmented_files))
    
    total_no_fragmented = df_non_fragmented.height
    total_fragmented = df_fragmented.height
    print(f"Skipped {total_no_dup - total_no_fragmented} fragmented entries")

    # Prepare output directory
    output_dir = Path(output_directory)
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)

    # Write output files
    print("Writing output files...")
    single_output_file = "smiles4cfm.txt"
    fragmented_smiles_file = "smiles_fragmented.txt"

    with open(single_output_file, "w") as all_f, open(fragmented_smiles_file, "w") as frag_f:
        for idx, (short_inchikey, smiles) in enumerate(zip(df_non_fragmented["short_inchikey"], df_non_fragmented["smiles"]), 1):
            output_file = output_dir / f"{str(idx).zfill(8)}_{short_inchikey}.txt"
            
            # Write individual file
            with open(output_file, "w") as f:
                f.write(f"{short_inchikey} {smiles}\n")
            
            # Append to the single non-fragmented file
            all_f.write(f"{short_inchikey} {smiles}\n")

            if idx % 1000 == 0:
                print(f"Progress: {idx}/{total_no_fragmented}")
        
        # Write fragmented SMILES
        for short_inchikey, smiles in zip(df_fragmented["short_inchikey"], df_fragmented["smiles"]):
            frag_f.write(f"{short_inchikey} {smiles}\n")

    # Create unique (item, short_inchikey) file
    df_unique_pairs = df_valid.select(["item", "short_inchikey"]).unique()
    df_unique_pairs.write_csv("short_inchikeys.tsv", separator="\t")

    end_time = time.time()
    total_processing_time = (end_time - start_time) / 60

    print(f"\nSummary:")
    print(f"Initial entries: {total_initial}")
    print(f"After removing dots: {total_no_dot}")
    print(f"After removing invalid: {total_no_invalid}")
    print(f"After removing duplicates: {total_no_dup}")
    print(f"After skipping fragmented: {total_no_fragmented}")
    print(f"Fragmented entries saved: {total_fragmented}")
    print(f"Output written to {output_directory}")
    print(f"\nTotal processing time: {total_processing_time:.2f} minutes")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a SMILES file and match against already fragmented ones.")
    parser.add_argument("input_file", type=str, help="Path to the input SMILES file")
    parser.add_argument("--output-directory", type=str, default="in", help="Directory to store entries to fragment (default: in)")
    parser.add_argument("--search-directory", type=str, default="out", help="Directory to search for already fragmented entries (default: out)")
    parser.add_argument("--skip-unfragmented", action="store_true", help="Skip checking file contents and only check filenames")
    
    args = parser.parse_args()

    process_smiles_file(
        input_file=args.input_file,
        output_directory=args.output_directory,
        search_directory=args.search_directory,
        skip_unfragmented=args.skip_unfragmented
    )
