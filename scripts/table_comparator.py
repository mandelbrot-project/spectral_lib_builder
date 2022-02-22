#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 16:08:34 2020

@author: pma
"""

from __future__ import print_function
import pandas as pd
import sys



try:
    tofrag_file_path = sys.argv[1]
    fragged_file_path = sys.argv[2]
    unfragged_file_path = sys.argv[3]
    
    
    print('Comparing to frag file' 
          + tofrag_file_path 
          + ' with fragged results ' 
          + fragged_file_path 
          + 'and outputting the difference in ' 
          + unfragged_file_path)
except:
    print('Please fill tofrag_file_path, fragged_file_path and unfragged_file_path')





df_tofrag = pd.read_csv(tofrag_file_path, sep=' ', header=None) 
df_fraged = pd.read_csv(fragged_file_path, sep='\t', header=None)
df_tofrag.rename(columns={0: 'NPAID', 1: 'SMILES'}, inplace=True)
df_merged = pd.merge(df_tofrag, df_fraged, left_on='NPAID', right_on=0, how='left')
df_unfragged = df_merged[df_merged[0].isna() == True]
df_unfragged = df_unfragged.drop(0, axis = 1)




df_unfragged.to_csv(unfragged_file_path, sep='\t', index=False)
