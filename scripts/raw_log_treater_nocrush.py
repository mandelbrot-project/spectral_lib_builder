
#%%

### This cleans raw .log cfm-predict output to yield clean spectral files for mgf header population

import glob
import re
import pandas as pd
import os
import sys

## define help menu 

try:
    directory_name=sys.argv[1]
    extension_type=sys.argv[2]
    print('Parsing directory' + directory_name + ' with file extension :' + extension_type)
except:
    print('Please pass parent directory, where the .log files are placed')


done = 0
skipped = 0

# this regex allows to suck only the two first columns of floats
p = re.compile(r'^^([0-9]*\.{1}[0-9]*) ([0-9]*\.{1}[0-9]*)', re.M)

## we iterate over each .log file of the following folder
#for filename in glob.iglob('../npatlas_data/results_npatlas/*.log', recursive=True):

for filename in glob.iglob(str(directory_name) + '/*/*' + str(extension_type), recursive=True):
    
    try:
        #we then get a first block corresponding to masses only
        coco_masses = open(filename,'r').read().split('\n\n')
        
        # and split this block according to te energy pattern
        coco_split = re.split('energy.*', str(coco_masses[0]))
        
        
        # since we want to keep all of them we append them to a unique df, order them by mass
        coco_list_energy_all = pd.DataFrame(re.findall(p, coco_split[1]))
        coco_list_energy_all = coco_list_energy_all.append(pd.DataFrame(re.findall(p, coco_split[2])))
        coco_list_energy_all = coco_list_energy_all.append(pd.DataFrame(re.findall(p, coco_split[3])))
        
        coco_list_energy_all.rename(columns={0:'mass'},inplace=True)
        
        coco_list_energy_all['mass'] = pd.to_numeric(coco_list_energy_all['mass'], errors='coerce')
        
        
        
        coco_list_energy_all = coco_list_energy_all.sort_values('mass', ascending = True)
        
    
        # uncomment if you want to output to a differnt filename
        mgf_filename = filename.split('.log')[0]+'.mgf'
        
        # be carefull here we'll crush the original files
        coco_list_energy_all.to_csv(mgf_filename,sep=' ', index=False, header=False)
        

        # os.remove(filename)
        
        done += 1
        print(filename, ' treated.')
        
    except Exception as er:
        print(er, filename)
        skipped += 1
    
print('Treated {} files, skipped {}'.format(done, skipped))


 
