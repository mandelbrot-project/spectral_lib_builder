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
import numpy as np

df = pd.read_csv("/Users/pma/is_fragmentation/npatlas_data/np_atlas_2019_12.tsv", sep='\t') 

df.head()
df.columns
df.info()

df['ROMol'] = df['SMILES'].map(Chem.MolFromSmiles)

#df['ROMol'] = df['InChI_DNP'].map(Chem.MolFromInchi)

#df['SMILES'] = df['ROMol'].map(Chem.MolToSmiles)


df['ROMol']
df.info()


df = df[~df['ROMol'].isnull()]

# if df['ROMol'].isnull() == True 
df['EMW'] = df['ROMol'].map(Descriptors.ExactMolWt)


df['EMW'].head()


df['protonated_emass'] = df['EMW'] + 1.007276
df['deprotonated_emass'] = df['EMW'] - 1.007276


df = df.drop('ROMol', axis=1)

df.to_csv("/Users/pma/is_fragmentation/npatlas_data/np_atlas_2019_12_adducted.tsv", sep='\t', index=False)

df.to_csv("/Users/pma/Desktop/DNP_SMILED.tsv", sep='\t', index=False)

df = df.fillna('NFound')

PandasTools.FrameToGridImage(df.head(8), legendsCol='NPAID', molsPerRow=4)


df.head()
df.columns


df = pd.read_csv("/Users/pma/Desktop/190602_DNP_TAXcof_CF.tsv", sep="\t")
