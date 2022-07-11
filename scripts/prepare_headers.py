#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 16:08:34 2020

@author: pma
"""

# Prepares headers for mgf header population
# from __future__ import print_function
from rdkit import Chem
import pandas as pd
from rdkit.Chem import PandasTools
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

import sys

try:
    metadata_file_path = sys.argv[1]
    delimiter_input = sys.argv[2]
    adducted_metadata_file_path = sys.argv[3]
    column_name_smiles = sys.argv[4]
    column_name_sik = sys.argv[5]
    print('Parsing file ' 
          + metadata_file_path 
          + ' and calculating exact mass, protonated and deprotonated adducts exact mass and appending to new metadata file ' 
          + adducted_metadata_file_path)
except:
    print('Please add input and output file path')    
    print(delimiter_ouput)

if delimiter_input == "s+":
    df = pd.read_csv(
        metadata_file_path, 
        names=[column_name_sik, column_name_smiles], 
        header=None,
        delim_whitespace=True) 
else:
    df = pd.read_csv(metadata_file_path, sep = delimiter_input) 

df['SMILES'] = df[column_name_smiles]
df['SHORT_IK'] = df[column_name_sik].replace('"', "")
df['ROMol'] = df['SMILES'].map(Chem.MolFromSmiles)

df = df[~df['ROMol'].isnull()]

df['EMW'] = df['ROMol'].map(Descriptors.ExactMolWt)
df['MF'] = df['ROMol'].map(rdMolDescriptors.CalcMolFormula)
df['InChI'] = df['ROMol'].map(Chem.MolToInchi)

df['protonated_emass'] = df['EMW'] + 1.007276
df['deprotonated_emass'] = df['EMW'] - 1.007276
df = df.drop('ROMol', axis=1)

df.to_csv(adducted_metadata_file_path, sep = '\t', index=False)
