import zipfile
import argparse
import numpy as np
import polars as pl
from collections import defaultdict
from matchms import Spectrum
from matchms.exporting import save_as_json, save_as_mgf, save_as_msp
from rdkit import Chem

def determine_adduct_and_charge(polarity):
    """Set adduct and charge based on polarity."""
    if polarity.lower() == "positive":
        return "[M+H]+", "1+"
    elif polarity.lower() == "negative":
        return "[M-H]-", "1-"
    else:
        raise ValueError("Polarity must be 'positive' or 'negative'")


def parse_out_file_from_content(lines, polarity):
    """Parse file content and extract metadata and spectra."""
    metadata = {}
    energy_spectra = {}
    current_energy = None
    adduct, charge = determine_adduct_and_charge(polarity)
    
    for line in lines:
        line = line.strip()
        
        # Stop parsing current energy level on a blank line
        if current_energy and line == "":
            current_energy = None
            continue
        
        # Extract metadata
        if line.startswith("#PREDICTED"):
            metadata["library"] = line.split("#", 1)[1]
        elif line.startswith("#ID="):
            metadata["id"] = line.split("=", 1)[1]
        elif line.startswith("#SMILES="):
            metadata["smiles"] = line.split("=", 1)[1]
        elif line.startswith("#InChiKey="):
            metadata["inchikey"] = line.split("=", 1)[1]
        elif line.startswith("#Formula="):
            metadata["formula"] = line.split("=", 1)[1]
        elif line.startswith("#PMass="):
            metadata["precursor_mz"] = float(line.split("=", 1)[1])
        elif line.startswith("energy"):
            current_energy = line
            energy_spectra[current_energy] = []
        elif current_energy and not line.startswith("#"):
            parts = line.split()
            if len(parts) >= 2:
                mz, intensity = map(float, parts[:2])
                energy_spectra[current_energy].append((mz, intensity))

    # Create matchms Spectrum objects for each energy level
    spectra_per_energy = {}
    for energy, peaks in energy_spectra.items():
        if peaks and "precursor_mz" in metadata:
            # Convert to numpy arrays
            mz_values = np.array([p[0] for p in peaks], dtype="float32")
            intensities = np.array([p[1] for p in peaks], dtype="float32")
            
            # Sort by m/z values
            sorted_indices = np.argsort(mz_values)
            mz_values = mz_values[sorted_indices]
            intensities = intensities[sorted_indices]
            
            spectrum_metadata = {
                "PEPMASS": metadata["precursor_mz"],
                "CHARGE": charge,
                "FILENAME": metadata["id"],
                "MOLECULAR_FORMULA": metadata["formula"],
                "IONMODE": polarity.upper(),
                "EXACTMASS": metadata["precursor_mz"] - 1.007276 if polarity == "positive" else metadata["precursor_mz"] + 1.007276,
                "NAME": metadata["id"],
                "SMILES": metadata["smiles"],
                "INCHI": Chem.MolToInchi(Chem.MolFromSmiles(metadata["smiles"])),
                "INCHIKEY": metadata["inchikey"],
                "LIBRARYQUALITY": metadata["library"],
                "SCANS": len(mz_values),
                "ADDUCT": adduct,
                "ENERGY": energy
            }
            
            spectrum = Spectrum(mz=mz_values, intensities=intensities, metadata=spectrum_metadata)
            spectra_per_energy[energy] = spectrum

    return spectra_per_energy


def merge_spectra(spectra):
    """Merge spectra per InChIKey by summing their intensities at matching m/z values."""
    merged_spectra = {}
    
    # Group spectra by InChIKey
    spectra_by_inchikey = defaultdict(list)
    for spec in spectra:
        inchikey = spec.metadata.get("inchikey", "unknown")
        spectra_by_inchikey[inchikey].append(spec)

    # Merge per InChIKey
    for inchikey, spec_list in spectra_by_inchikey.items():
        all_mz = np.concatenate([spec.peaks.mz for spec in spec_list])
        all_intensities = np.concatenate([spec.peaks.intensities for spec in spec_list])

        unique_mz = np.unique(all_mz)
        merged_intensities = np.array([np.sum(all_intensities[all_mz == mz]) for mz in unique_mz], dtype="float32")

        # Normalize intensities to a max of 100
        max_intensity = merged_intensities.max()
        if max_intensity > 0:
            merged_intensities = (merged_intensities / max_intensity) * 100

        # Preserve metadata (taking from the first spectrum)
        base_metadata = spec_list[0].metadata.copy()
        base_metadata["energy"] = "energySum"
        base_metadata["scans"] = len(unique_mz)

        merged_spectra[inchikey] = Spectrum(
            mz=np.array(unique_mz, dtype="float32"),
            intensities=np.array(merged_intensities, dtype="float32"),
            metadata=base_metadata
        )

    return list(merged_spectra.values())


def process_zip(zip_path, folder_in_zip, output_prefix, polarity, query_file, inchikey_file):
    """Process .log files from a directory inside a ZIP archive."""
    all_spectra = []
    spectra_by_energy = {}

    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        file_list = [f for f in zip_ref.namelist() if f.startswith(folder_in_zip) and f.endswith(".log")]
        
        for file in file_list:
            with zip_ref.open(file) as f:
                content = f.read().decode("utf-8").splitlines()
                spectra_per_energy = parse_out_file_from_content(content, polarity)
                for energy, spectrum in spectra_per_energy.items():
                    spectra_by_energy.setdefault(energy, []).append(spectrum)
                    all_spectra.append(spectrum)

    # Read query and short InChIKey files
    query_df = pl.read_csv(query_file, separator="\t").select("item")
    inchikey_df = pl.read_csv(inchikey_file, separator="\t")
    filtered_keys = query_df.join(inchikey_df, on="item", how="inner")["short_inchikey"].to_list()
    lotus_all_spectra = [s for s in all_spectra if s.metadata.get("name") in filtered_keys]

    polarity_short = polarity[:3].lower()

    # Save per-energy files for wikidata
    for energy, spectra in spectra_by_energy.items():
        save_as_json(spectra, f"{output_prefix}_wikidata_{polarity_short}_{energy}.json")
        save_as_mgf(spectra, f"{output_prefix}_wikidata_{polarity_short}_{energy}.mgf")
        save_as_msp(spectra, f"{output_prefix}_wikidata_{polarity_short}_{energy}.msp")
    
    # Save combined all-energy file for wikidata
    save_as_json(all_spectra, f"{output_prefix}_wikidata_{polarity_short}_energyAll.json")
    save_as_mgf(all_spectra, f"{output_prefix}_wikidata_{polarity_short}_energyAll.mgf")
    save_as_msp(all_spectra, f"{output_prefix}_wikidata_{polarity_short}_energyAll.msp")

    # Save merged summed spectrum for wikidata
    merged_spectra = merge_spectra(all_spectra)
    save_as_json(merged_spectra, f"{output_prefix}_wikidata_{polarity_short}_energySum.json")
    save_as_mgf(merged_spectra, f"{output_prefix}_wikidata_{polarity_short}_energySum.mgf")
    save_as_msp(merged_spectra, f"{output_prefix}_wikidata_{polarity_short}_energySum.msp")

    # Save per-energy files for lotus
    for energy, spectra in spectra_by_energy.items():
        filtered_spectra_energy = [s for s in spectra if s.metadata.get("name") in filtered_keys]
        save_as_json(filtered_spectra_energy, f"{output_prefix}_lotus_{polarity_short}_{energy}.json")
        save_as_mgf(filtered_spectra_energy, f"{output_prefix}_lotus_{polarity_short}_{energy}.mgf")
        save_as_msp(filtered_spectra_energy, f"{output_prefix}_lotus_{polarity_short}_{energy}.msp")
    
    # Save combined all-energy file for lotus
    save_as_json(lotus_all_spectra, f"{output_prefix}_lotus_{polarity_short}_energyAll.json")
    save_as_mgf(lotus_all_spectra, f"{output_prefix}_lotus_{polarity_short}_energyAll.mgf")
    save_as_msp(lotus_all_spectra, f"{output_prefix}_lotus_{polarity_short}_energyAll.msp")

    # Save merged summed spectrum for lotus
    merged_lotus_spectra = merge_spectra(lotus_all_spectra)
    save_as_json(merged_lotus_spectra, f"{output_prefix}_lotus_{polarity_short}_energySum.json")
    save_as_mgf(merged_lotus_spectra, f"{output_prefix}_lotus_{polarity_short}_energySum.mgf")
    save_as_msp(merged_lotus_spectra, f"{output_prefix}_lotus_{polarity_short}_energySum.msp")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process .log files from ZIP into MGF/MSP with filtering")
    parser.add_argument("zip_path", help="Path to the ZIP file")
    parser.add_argument("folder_in_zip", help="Folder inside ZIP containing .log files")
    parser.add_argument("output_prefix", help="Prefix for output files")
    parser.add_argument("--polarity", choices=["positive", "negative"], required=True, help="Polarity")
    parser.add_argument("--query_file", default="query_p703.tsv", help="TSV file with query items")
    parser.add_argument("--inchikey_file", default="short_inchikeys.tsv", help="TSV file with short InChIKeys")
    
    args = parser.parse_args()
    process_zip(args.zip_path, args.folder_in_zip, args.output_prefix, args.polarity, args.query_file, args.inchikey_file)
