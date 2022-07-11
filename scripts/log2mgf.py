#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 16:08:34 2020

@author: pma
"""

# Cleans raw .log cfm-predict output to yield clean spectral files for mgf header population
import glob
import re
import pandas as pd
import os
import sys

# define help menu 
try:
    directory_name=sys.argv[1]
    extension_type=sys.argv[2]
    manipulation=sys.argv[3]
    print('Parsing directory ' + directory_name + ' with file extension: ' + extension_type)
except:
    print('Please pass parent directory, where the .log files are placed')

done = 0
skipped = 0

# this regex allows to suck only the two first columns of floats
p = re.compile(r'^^([0-9]*\.{1}[0-9]*) ([0-9]*\.{1}[0-9]*)', re.M)

# we iterate over each .log file of the following folder
for filename in glob.iglob(str(directory_name) + '*' + str(extension_type), recursive=True):
    
    try:
        # we then get a first block corresponding to masses only
        masses = open(filename,'r').read().split('\n\n')
        
        # and split this block according to te energy pattern
        split = re.split('energy.*', str(masses[0]))
        
        # since we want to keep all of them we append them to a unique df, order them by mass
        list_energy_all = pd.DataFrame(re.findall(p, split[1]))
        list_energy_all = pd.concat([list_energy_all, pd.DataFrame(re.findall(p, split[2]))])
        list_energy_all = pd.concat([list_energy_all, pd.DataFrame(re.findall(p, split[3]))])
        list_energy_all.rename(columns={0:'mass'},inplace=True)
        list_energy_all['mass'] = pd.to_numeric(list_energy_all['mass'], errors='coerce')
        list_energy_all = list_energy_all.astype(float)

        # manipulation
        if manipulation == 'mean':
            # must be one of ['mean','max','sum']
            list_energy_all_manipulated = list_energy_all.groupby(['mass']).mean()
        elif manipulation == 'max':
            # must be one of ['mean','max','sum']
            list_energy_all_manipulated = list_energy_all.groupby(['mass']).max()
        elif manipulation == 'sum':
            # must be one of ['mean','max','sum']
            list_energy_all_manipulated = list_energy_all.groupby(['mass']).sum()
        else:
            list_energy_all_manipulated = list_energy_all['mass']

        list_energy_all_manipulated = list_energy_all_manipulated.reset_index()

        list_energy_all = list_energy_all_manipulated.sort_values('mass', ascending = True)
    
        # uncomment if you want to output to a different filename
        mgf_filename = filename.split('.log')[0].replace('"', "")+'.mgf'
        
        list_energy_all.to_csv(mgf_filename,sep=' ', index=False, header=False)
        # be carefull if you remove the original files
        # os.remove(filename)
        done += 1
        
    except Exception as er:
        print(er, filename)
        skipped += 1
    
print('Treated {} files, skipped {}'.format(done, skipped))
