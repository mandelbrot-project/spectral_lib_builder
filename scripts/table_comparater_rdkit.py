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

df_tofrag = pd.read_csv("/Users/pma/is_fragmentation/npatlas_data/npatlas_for_frag.txt", sep=' ') 
df_fraged = pd.read_csv("/Users/pma/is_fragmentation/npatlas_data/list_mgf.txt", sep='\t', header=None) 

df_merged = pd.merge(df_tofrag, df_fraged, left_on='NPAID', right_on=0, how='left')

df_unfragged = df_merged[df_merged[0].isna() == True]
df_unfragged = df_unfragged.drop(0, axis = 1)




df_unfragged.to_csv("/Users/pma/is_fragmentation/npatlas_data/unfragged_compounds.txt", sep='\t', index=False)

df = df_unfragged


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

df.sort_values(by=['EMW'], ascending = False)

df.to_csv("/Users/pma/is_fragmentation/npatlas_data/unfragged_compounds.txt", sep='\t', index=False)





""" contribution from Hans de Winter """
from rdkit import Chem
from rdkit.Chem import AllChem

def _InitialiseNeutralisationReactions():
    patts= (
        # Imidazoles
        ('[n+;H]','n'),
        # Amines
        ('[N+;!H0]','N'),
        # Carboxylic acids and alcohols
        ('[$([O-]);!$([O-][#7])]','O'),
        # Thiols
        ('[S-;X1]','S'),
        # Sulfonamides
        ('[$([N-;X2]S(=O)=O)]','N'),
        # Enamines
        ('[$([N-;X2][C,N]=C)]','N'),
        # Tetrazoles
        ('[n-]','[nH]'),
        # Sulfoxides
        ('[$([S-]=O)]','S'),
        # Amides
        ('[$([N-]C=O)]','N'),
        )
    return [(Chem.MolFromSmarts(x),Chem.MolFromSmiles(y,False)) for x,y in patts]

_reactions=None
def NeutraliseCharges(smiles, reactions=None):
    global _reactions
    if reactions is None:
        if _reactions is None:
            _reactions=_InitialiseNeutralisationReactions()
        reactions=_reactions
    mol = Chem.MolFromSmiles(smiles)
    replaced = False
    for i,(reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant):
            replaced = True
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]
    if replaced:
        return (Chem.MolToSmiles(mol,True), True)
    else:
        return (smiles, False)
    
_reactions=None
def NeutraliseCharges_series(smiles, reactions=None):
    global _reactions
    if reactions is None:
        if _reactions is None:
            _reactions=_InitialiseNeutralisationReactions()
        reactions=_reactions
    mol = Chem.MolFromSmiles(smiles)
    replaced = False
    for i,(reactant, product) in enumerate(reactions):
        while mol.HasSubstructMatch(reactant):
            replaced = True
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]
    if replaced:
        return pd.Series([Chem.MolToSmiles(mol,True), True])
    else:
        return pd.Series([smiles, False])
    
    
smis=("c1cccc[nH+]1",
      "C[N+](C)(C)C","c1ccccc1[NH3+]",
      "CC(=O)[O-]","c1ccccc1[O-]",
      "CCS",
      "C[N-]S(=O)(=O)C",
      "C[N-]C=C","C[N-]N=C",
      "c1ccc[n-]1",
      "CC[N-]C(=O)CC")
for smi in smis:
    (molSmiles, neutralised) = NeutraliseCharges(smi)
    print(smi + "->" + molSmiles)
    
smis=("c1cccc[nH+]1",
      "C[N+](C)(C)C","c1ccccc1[NH3+]",
      "CC(=O)[O-]","c1ccccc1[O-]",
      "CCS",
      "C[N-]S(=O)(=O)C",
      "C[N-]C=C","C[N-]N=C",
      "c1ccc[n-]1",
      "CC[N-]C(=O)CC")
for smi in smis:
    (molSmiles, neutralised) = NeutraliseCharges(smi)
    print(smi + "->" + molSmiles)
    
    
for i in list(df['SMILES']):
    n = NeutraliseCharges(i)
    df['Nautralized_SMILE'] = n[0]
    df['Nautralization_status'] = n[1]
    print(n[0])




df['Neutralized_SMILE'] = df['SMILES'].apply(NeutraliseCharges)

df['SMILES'].apply(NeutraliseCharges)

df[['sum', 'difference']] = df.apply(
    lambda row: pd.Series(add_subtract(row['a'], row['b'])), axis=1)

Chem.SanitizeMol(m)

df['Sanitized'] = df['ROMol'].apply(Chem.SanitizeMol)


df[['Neutralized_SMILE', 'Neutralization_status']] = df['SMILES'].apply(NeutraliseCharges_series)


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
