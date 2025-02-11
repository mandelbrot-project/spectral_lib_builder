import os
import shutil
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm


def check_and_move(entry, destination_dir, progress_bar):
    """Checks if 'energy2' is in the file. If not, move it to the destination."""
    try:
        with open(entry.path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if "energy2" in line:
                    progress_bar.update(1)  # Update progress and exit early
                    return  
        # Move file if 'energy2' was not found
        shutil.move(entry.path, os.path.join(destination_dir, entry.name))
    except Exception as e:
        print(f"Error processing {entry.name}: {e}")
    finally:
        progress_bar.update(1)  # Ensure progress updates even on errors


def process_files(source_dir, destination_dir, max_workers=8):
    """Scans files in source_dir and moves those without 'energy2' to destination_dir."""
    os.makedirs(destination_dir, exist_ok=True)

    # Collect all files first to know the total count
    files = [entry for entry in os.scandir(source_dir) if entry.is_file()]
    total_files = len(files)

    if total_files == 0:
        print("No files found in the source directory.")
        return

    # Use tqdm for a progress bar
    with tqdm(total=total_files, desc="Processing files", unit="file") as progress_bar:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(check_and_move, entry, destination_dir, progress_bar) for entry in files]

            # Ensure progress is updated even with multithreading
            for future in as_completed(futures):
                future.result()  # Handle any exceptions raised in threads


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Move files not containing 'energy2' to another folder.")
    parser.add_argument("--source", default="out", help="Source directory (default: out)")
    parser.add_argument("--destination", default="quarantine", help="Destination directory (default: quarantine)")
    parser.add_argument("--workers", type=int, default=8, help="Number of worker threads (default: 8)")
    
    args = parser.parse_args()
    process_files(args.source, args.destination, args.workers)
