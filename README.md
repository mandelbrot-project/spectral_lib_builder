# Tutorial for NPatlas - OpenNPDB in silico fragmentation data treatment


## Prior to fragmentation 

Complete metadate file is converted to list of Unique ID and smiles space separated (for cfm id input)

`python frag_list_preparator.py ../npatlas_data/np_atlas_2019_12.tsv ../npatlas_data/np_atlas_2019_12_for_frag.tsv NPAID SMILES`
`python scripts/frag_list_preparator.py ./open_np_db_data/open_NP_db.tsv ./open_np_db_data/open_NP_db_for_frag.txt shortinchikey shortinchikey smiles`
`python scripts/frag_list_preparator.py ../../../../210505_lotus_dnp_single_sik.csv ./lotus_data/lotus_data_for_frag.txt short_inchikey short_inchikey structure_smiles`


## Split the file

(Works on Linux based shells)

split -a 5 -l 500 -d ../open_np_db_data/open_NP_db_tofrag.txt ../open_np_db_data/opennpdb_tofrag/opennpdb_
split -a 5 -l 500 -d ./lotus_data/lotus_data_for_frag.txt ./lotus_data/lotus_data_for_frag/lotus_data_



## Prepare mutilple bash file to launch on baobab

Be careful the nodes partition names have changed.



## Fetching cfm-predict results

Download cfm-predict fragmentation results from the baob server using rsync command

`rsync -rvz -e 'ssh' --progress allardp@baobab2.unige.ch:/home/allardp/CFM_results/npatlas ./results`
`rsync -rvz -e 'ssh' --progress allardp@baobab2.unige.ch:/home/allardp/CFM_results/opennpdb ./results`


This code line : 

`find ./ -type f -name '*.mgf' | wc  `

allows to count all file in a folder. Here 25090 files for NPatlas


## Pruning the raw log files

The output of cfm-predict consist of .log file containing mass spectra, where each fragments are individually labelled and eventually linked to a substrcture. Such information might be usefull later but for now we only want to keep the raw ms data

(need to define a help function here)

`python raw_log_treater_npatlas.py ../npatlas_data/results_npatlas/npatlas/ .log`

At this step .log file should be pruned and contains only digits (m/z and intensities)




## Populating the mgf headers

### Preparation of the adducted metadata table

We need to prepare and adducted dataframe containing the protonated and deprotonated masses

This script recquire rdkit so we build a environment.yml file from a dedicated conda env

`conda env export -n conda-env -f /path/to/environment.yml`

`python table_adducter_npatlas_script.py ../npatlas_data/np_atlas_2019_12.tsv ../npatlas_data/np_atlas_2019_12_adducted.tsv`


### Addition of the metadata to the individual mgf headers

We can now populate each raw mgf with its corresponding metadata. For this we use the treat_npatlas.py script

`python treat_npatlas.py ../npatlas_data/np_atlas_2019_12_adducted.tsv ../npatlas_data/results_npatlas/npatlas/`


## Generating the final spectral file

We concatenate each documented mgf files to a single spectral mgf file.

`find ./ -type f -name '*.mgf' | while read F; do cat ${F} >> ../../npatlas_ISDB_pos.mgf; done`


## Outputting non-fragmented entries

For several reasons (charged compounds, some tautomers, structures too heavy to be fragmented in a reasonable amount of time) some entries might not have been fragmented. 

To find them we will first list all correctly converted mgf

`find ./ -type f -name '*.mgf' | sed 's!.*/!!' | sed 's!^!!' >  list_mgf.txt`

or here eventually without the extension

`find ./ -type f -name '*.mgf' | sed 's!.*/!!' | sed 's!.mgf!!' >  ../../list_mgf.txt`

And then the list is compared to the initial input using the table_comparator.py 

`python table_comparator.py ../npatlas_data/npatlas_for_frag.txt ../npatlas_data/list_mgf.txt ../npatlas_data/unfragged_list.txt`


### check molVS for structure standardization

https://molvs.readthedocs.io/en/latest/index.html


# Results

## NPAtlas

Fragmented NPAtlas spectral file (.mgf) and associated metadata are available here [NPAtlas_ISDB](https://www.dropbox.com/sh/rz9giwvzuhnvlpo/AABLJIu2EKo7pJrP-ALHFbfua?dl=0)




