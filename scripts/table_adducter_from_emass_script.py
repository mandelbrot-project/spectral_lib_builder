#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 16:08:34 2020

@author: pma
"""

from __future__ import print_function
# from rdkit import Chem
import pandas as pd
# from rdkit.Chem import PandasTools
# from rdkit.Chem import Descriptors
import sys


try:
    metadata_file_path = sys.argv[1]
    delimiter_input = sys.argv[2]
    emass_field = sys.argv[3]
    adducted_metadata_file_path = sys.argv[4]
    delimiter_ouput = sys.argv[5]
    
    
    print('Parsing file' 
          + metadata_file_path 
          + ' and calculating exact mass, protonated and deprotonated adducts exact mass and appending to new metadatfile' 
          + adducted_metadata_file_path)
except:
    print('Please add input and output file path')
    


df = pd.read_csv(metadata_file_path, sep = delimiter_input) 


# df['ROMol'] = df['SMILES'].map(Chem.MolFromSmiles)

# df = df[~df['ROMol'].isnull()]


# df['EMW'] = df['ROMol'].map(Descriptors.ExactMolWt)


df['protonated_emass'] = df[emass_field] + 1.007276
df['deprotonated_emass'] = df[emass_field] - 1.007276
# df = df.drop('ROMol', axis=1)


df.to_csv(adducted_metadata_file_path, sep = delimiter_ouput, index=False)
