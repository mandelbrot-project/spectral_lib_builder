#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 16:08:34 2020

@author: pma
"""

from __future__ import print_function
from rdkit import Chem
import pandas as pd
from rdkit.Chem import PandasTools
from rdkit.Chem import Descriptors
import sys


try:
    metadata_file_path = sys.argv[1]
    adducted_metadata_file_path = sys.argv[2]
    
    
    print('Parsing file' 
          + metadata_file_path 
          + ' and calculating exact mass, protonated and deprotonated adducts exact mass and appending to new metadatfile' 
          + adducted_metadata_file_path)
except:
    print('Please add input and output file path')
    


df = pd.read_csv(metadata_file_path, sep='\t') 


df['ROMol'] = df['SMILES'].map(Chem.MolFromSmiles)


df['EMW'] = df['ROMol'].map(Descriptors.ExactMolWt)


df['protonated_emass'] = df['EMW'] + 1.007276
df['deprotonated_emass'] = df['EMW'] - 1.007276
df = df.drop('ROMol', axis=1)


df.to_csv(adducted_metadata_file_path, sep='\t', index=False)
