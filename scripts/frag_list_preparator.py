#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 16:08:34 2020

@author: pma
"""

import pandas as pd
import sys


try:
    metadata_file_path = sys.argv[1]
    for_frag_file_path = sys.argv[2]
    duplicated_field_header = sys.argv[3]
    id_header = sys.argv[4]
    smiles_header = sys.argv[5]
    
    
    
    print('Parsing file ' 
          + metadata_file_path 
          + ' removing duplicated values in field: '
          + duplicated_field_header
          + ' and keeping unique ID column : '
          + id_header
          + ' and SMILES columns: '
          + smiles_header
          + ' as a space separated text file: ' 
          + for_frag_file_path)
except:
    print('Please add input and output file path followed by the header of i) the field to deduplicated, ii) the unique id field iii) the smiles code field')
    


df = pd.read_csv(metadata_file_path, sep=',', error_bad_lines=False, low_memory=False)


df = df.drop_duplicates(subset= duplicated_field_header )


df = df[[id_header, smiles_header]]


df.to_csv(for_frag_file_path, sep=' ', index=False, header=None)
