import zipfile
import argparse
import numpy as np
import polars as pl
from collections import defaultdict
from matchms import Spectrum
from matchms.exporting import save_as_json, save_as_mgf, save_as_msp
from rdkit import Chem


def determine_adduct_and_charge(polarity):
    return ("[M+H]+", "1+") if polarity == "positive" else ("[M-H]-", "1-")


def parse_out_file_from_content(lines, polarity):
    metadata = {}
    energy_spectra = defaultdict(list)
    adduct, charge = determine_adduct_and_charge(polarity)
    current_energy = None
    
    for line in map(str.strip, lines):
        if not line:
            current_energy = None
            continue
        
        if line.startswith("#PREDICTED"):
            metadata["library"] = line.lstrip("#").strip()
        elif line.startswith("#"):
            key, _, value = line.partition("=")
            metadata[key.lstrip("#")] = value
        elif line.startswith("energy"):
            current_energy = line
        elif current_energy:
            try:
                mz, intensity = map(float, line.split()[:2])
                energy_spectra[current_energy].append((mz, intensity))
            except ValueError:
                pass
    
    precursor_mz = float(metadata.get("PMass", 0))
    exact_mass = precursor_mz - 1.007276 if polarity == "positive" else precursor_mz + 1.007276
    spectra = {}
    
    for energy, peaks in energy_spectra.items():
        if peaks and precursor_mz:
            mz_values, intensities = zip(*sorted(peaks))
            spectrum = Spectrum(
                mz=np.array(mz_values, dtype="float32"),
                intensities=np.array(intensities, dtype="float32"),
                metadata={
                    "PEPMASS": precursor_mz,
                    "CHARGE": charge,
                    "IONMODE": polarity.upper(),
                    "COLLISION_ENERGY": energy,
                    "ADDUCT": adduct,
                    "LIBRARYQUALITY": metadata.get("library"),
                    "FRAGMENTATION_MODE": "in silico",
                    "EXACTMASS": exact_mass,
                    "NAME": metadata.get("ID"),
                    "MOLECULAR_FORMULA": metadata.get("Formula"),
                    "SMILES": metadata.get("SMILES"),
                    "INCHI": Chem.MolToInchi(Chem.MolFromSmiles(metadata.get("SMILES", ""))) if metadata.get("SMILES") else "",
                    "INCHIKEY": metadata.get("InChiKey"),
                    "SCANS": 1,
                    "NUM_PEAKS": len(mz_values),
                }
            )
            spectra[energy] = spectrum
    return spectra


def merge_spectra(spectra):
    merged_spectra = {}
    
    spectra_by_id = defaultdict(list)
    for spec in spectra:
        _id = spec.metadata.get("compound_name", "unknown")
        spectra_by_id[_id].append(spec)

    for _id, spec_list in spectra_by_id.items():
        all_mz = np.concatenate([s.peaks.mz for s in spec_list])
        all_intensities = np.concatenate([s.peaks.intensities for s in spec_list])
        
        unique_mz = np.unique(all_mz)
        merged_intensities = np.array([np.sum(all_intensities[all_mz == mz]) for mz in unique_mz], dtype="float32")

        max_intensity = merged_intensities.max()
        if max_intensity > 0:
            merged_intensities = (merged_intensities / max_intensity) * 100
        
        base_metadata = spec_list[0].metadata.copy()
        base_metadata["collision_energy"] = "energySum"
        base_metadata["scans"] = len(spec_list)
        base_metadata["num_peaks"] = len(unique_mz)

        merged_spectra[_id] = Spectrum(
            mz=np.array(unique_mz, dtype="float32"),
            intensities=np.array(merged_intensities, dtype="float32"),
            metadata=base_metadata
        )
    
    return list(merged_spectra.values())


def process_zip(zip_path, folder_in_zip, output_prefix, polarity, query_file, inchikey_file):
    spectra, spectra_by_energy = [], defaultdict(list)
    
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        for file in filter(lambda f: f.startswith(folder_in_zip) and f.endswith(".log"), zip_ref.namelist()):
            with zip_ref.open(file) as f:
                spectra_per_energy = parse_out_file_from_content(f.read().decode().splitlines(), polarity)
                for energy, spectrum in spectra_per_energy.items():
                    spectra_by_energy[energy].append(spectrum)
                    spectra.append(spectrum)

    query_df, inchikey_df = pl.read_csv(query_file, separator="\t"), pl.read_csv(inchikey_file, separator="\t")
    filtered_keys = set(query_df.join(inchikey_df, on="item", how="inner")["short_inchikey"].to_list())
    
    lotus_spectra = [s for s in spectra if s.metadata.get("compound_name") in filtered_keys]
    merged_spectra, merged_lotus_spectra = merge_spectra(spectra), merge_spectra(lotus_spectra)

    polarity_short = polarity[:3].lower()


    def save_files(prefix, data):
        save_as_json(data, f"{prefix}.json")
        save_as_mgf(data, f"{prefix}.mgf")
        save_as_msp(data, f"{prefix}.msp")


    def sort_spectra(spectra_list, merged_spectra):
        merged_peaks = {s.metadata.get("compound_name", ""): s.metadata.get("num_peaks", 0) for s in merged_spectra}

        return sorted(spectra_list, key=lambda s: (
            merged_peaks.get(s.metadata.get("compound_name", ""), float('inf')),
            s.metadata.get("compound_name", ""),
            s.metadata.get("collision_energy", "")
        ))


    all_spectra = sort_spectra(spectra + merged_spectra, merged_spectra)
    merged_spectra = sort_spectra(merged_spectra, merged_spectra)
    
    all_lotus_spectra = sort_spectra(lotus_spectra + merged_lotus_spectra, merged_lotus_spectra)
    merged_lotus_spectra = sort_spectra(merged_lotus_spectra, merged_lotus_spectra)

    save_files(f"{output_prefix}_wikidata_{polarity_short}_energyAll", all_spectra)
    save_files(f"{output_prefix}_wikidata_{polarity_short}_energySum", merged_spectra)
    save_files(f"{output_prefix}_lotus_{polarity_short}_energyAll", all_lotus_spectra)
    save_files(f"{output_prefix}_lotus_{polarity_short}_energySum", merged_lotus_spectra)


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
