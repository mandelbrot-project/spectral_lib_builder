#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 16:08:34 2020

@author: pma
"""

import csv
import sys

# define help menu 
try:
    metadata_file_path = sys.argv[1]
    mass_files_folder_path = sys.argv[2]
    print('Parsing file'
        + metadata_file_path
        + ' and appending metadata to mass raw files loacted in '
        + mass_files_folder_path)
except:
    print('Please fill metadata file path followed by mass_files_folder_path')

# This will contain the tsv data
data = {}
header = False  # Was header read?

# Read the tsv file and store data in data
with open(metadata_file_path, 'r') as tsv_file:
    tsv_reader = csv.reader(tsv_file, delimiter='\t', quotechar='"')

    try:
        for row in tsv_reader:
            # Set the header if it wasn't yet
            if header is False:
                header = row
                continue  # Go to the next row
            if row != []:  # Take care of shitty windows format
                data[row[header.index('NPAID')]] = row
    except csv.Error as e:
        sys.exit('file {}, line {}: {}'.format(metadata_file_path,
                                               tsv_reader.line_num,
                                               e))

done = 0
skipped = 0

# Iterate over all files
for pos, row in data.items():
    # This will contain the mgf data
    content = ''
    print(row[header.index('NPAID')])
    filename = row[header.index('NPAID')]
    if ".mgf" not in filename:
        filename = "%s.mgf" % filename
    # Create content
    # To access a row of a specific name:
    #   row[header.index('NAMEOFTHEROW')]
    # Careful, as the program will crash if the row doesn't
    # exist (left as an exercise ;) )

    content += "BEGIN IONS\n"
    content += "PEPMASS={}\n".format(row[header.index('protonated_emass')])
    content += "CHARGE=1+\n"
    # MSLEVEL=xxx
    #content += "SOURCE_INSTRUMENT={}-{}\n".format(
        #row[header.index('INSTRUMENT')],
       # row[header.index('IONSOURCE')])
    content += "FILENAME={}\n".format(row[header.index('NPAID')])
    #content += "InChIKey={}\n".format(row[header.index('InChIKey')])
    content += "MOLECULAR_FORMULA={}\n".format(row[header.index('Molecular Formula')])
    content += "SEQ=*..*\n"
    #content += "NOTES={}:{}:{}:{}:{}:{}\n".format(
        #row[header.index('PI')],
        #row[header.index('DATACOLLECTOR')],
        # "N/A",  # Change that
        # "N/A",  # Change that
        # row[header.index('ACQUISITION')],
        # "N/A"  # Change that
        # )
    content += "IONMODE=POSITIVE\n"
    content += "EXACTMASS={}\n".format(row[header.index('EMW')])
    content += "NAME={}\n".format(row[header.index('NPAID')])
    content += "SMILES={}\n".format(row[header.index('SMILES')])
    content += "INCHI={}\n".format(row[header.index('InChI')])
    #content += "Synonyms={}\n".format(row[header.index('Syn_list')])
    #content += "LIB={}\n".format(row[header.index('LIB')])    
    #content += "INCHIAUX={}\n".format(row[header.index('INCHI_KEY')])
    content += "LIBRARYQUALITY=In Silico Fragmented CFM-ID\n"
    # Change that
    # ORGANISM=xxx
    # SPECTRUMID=xxx
    # ACTIVATION=xxx
    # INSTRUMENT=xxx
    # TITLE=xxx
    # content += "SCANS={}\n".format(row[header.index('EXTRACTSCAN')])
    content += "SCANS=1\n"

    # Try to open the file
    # Be careful given that the path below is relative be sure to check were you launch the script from   
    try:
        #with open('../npatlas_data/results_npatlas/{}'.format(filename), 'r') as raw_file:
        with open(str(mass_files_folder_path) + '{}'.format(filename), 'r') as raw_file:
            for line in raw_file:
                content += line.strip() + "\n"
            content += '\n'.join(list(raw_file))
            content += "END IONS\n"            
            # Write the output file
            # with open(filename.split('.')[0]+'.mgf', 'w') as mgf_file:
            with open(str(mass_files_folder_path) + '{}'.format(filename), 'w') as mgf_file:
                mgf_file.write(content)
            done += 1

    except IOError:
        print('File {} not found, skipped'.format(filename))
        skipped += 1
        continue

print('Treated {} files, skipped {}.'.format(done, skipped))
