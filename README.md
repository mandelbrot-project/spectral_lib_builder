# In Silico DataBase

## Foreword

All small python scripts require `environment.yml` to be installed to work.

## Prior to fragmentation 

Prepare a list of identified SMILES to fragment.
For the moment, we took the structures from [https://doi.org/10.5281/zenodo.6378223](https://doi.org/10.5281/zenodo.6378223) as starting point (a R solution is given at [scripts/lotus2cfm.R](scripts/lotus2cfm.R) .

//TODO Translate it for python afficionados (2022-07-07 AR)

## Prepare CFM on cluster

### Install

This installation procedure works on the [UniGE HPC](https://www.unige.ch/eresearch/en/services/hpc/). This does not mean it will work on another.


First, create a `cmf-4` directory:
```
mkdir cfm-4
```

Then load requested modules:
```
module load GCC/6.3.0-2.27 Singularity/2.4.2
```

Build cfm-4 from its last Docker image:
```
singularity build cfm-4/cfm.sif docker://wishartlab/cfmid
```

### Test

Pull the previsouly generated smiles list (this command is not generic, it needs to be adapted):
```
scp smiles4cfm.txt rutza@login2.baobab.hpc.unige.ch:smiles.txt
```

Also pull last bash commands:
```
scp scripts/run_cfm_test.sh rutza@login2.baobab.hpc.unige.ch:run_cfm_test.sh
scp scripts/run_cfm.sh rutza@login2.baobab.hpc.unige.ch:run_cfm.sh
scp scripts/run_cfm_neg.sh rutza@login2.baobab.hpc.unige.ch:run_cfm_neg.sh
```

Create a test file with 10 structures to check if everything works fine:
```
head smiles.txt -n 10 > test.txt
```

Create a `test` directory:
```
mkdir test
```

Split the test file asit would be split if real:
```
split --lines=1 --numeric-suffixes=1 --suffix-length=4 --additional-suffix=.txt test.txt test/test-
```

Create a `testout` directory:
```
mkdir testout
```

Run [run_cfm_test.sh](scripts/run_cfm_test.sh) in a sbatch array:
```
sbatch --array=1-10 run_cfm_test.sh
```

## Run CFM on cluster (full)

### Split in subfiles

Depending on the length of your SMILES list, you may want to split it:

First, create a `smiles` directory:
```
mkdir smiles
```

Split the big file:
```
split --lines=100 --numeric-suffixes=1 --suffix-length=4 --additional-suffix=.txt smiles.txt smiles/smiles-
```

Count how many entries it generated for the next step
```
ls smiles/ | wc -l
```

### Positive

Create the `posout` directory:
```
mkdir posout
```

Run [run_cfm.sh](scripts/run_cfm.sh) in a sbatch array (**adapt the array length**):
```
sbatch --array=1-2875 run_cfm.sh
```

### Negative

Create the `negout` directory:
```
mkdir negout
```

Run [run_cfm_neg.sh](scripts/run_cfm_neg.sh) in a sbatch array (**adapt the array length**):
```
sbatch --array=1-2875 run_cfm_neg.sh
```

## Fetch CFM results

Download CFM fragmentation results from the baob server (this command is not generic, it needs to be adapted):
```
scp -r rutza@login2.baobab.hpc.unige.ch:posout ./results
scp -r rutza@login2.baobab.hpc.unige.ch:negout ./results_neg
```

## Treating the raw log files

The output of cfm-predict consist of .log file containing mass spectra, where each fragments are individually labelled and eventually linked to a substrcture. 
Such information might be usefull later but for now we only want to keep the raw ms data:
If you want to merge the three different energies you can choose between 'max','mean', and 'sum' for the moment.
```
python scripts/log2mgf.py YOUR_RAW_RESULTS_DIR/ log sum
```

## Populating the mgf headers

### Preparation of the headers

We need to prepare and adducted table containing the protonated and deprotonated masses:
```
python scripts/prepare_headers.py YOUR_SMILES_LIST YOUR_DELIMITER YOUR_OUTPUT_PATH YOUR_SMILES_COLUMN_NAME YOUR_SHORT_INCHIKEY_COLUMN_NAME
```

### Addition of the metadata to the individual mgf headers

We can now populate each mgf with its corresponding metadata:
```
 python scripts/populate_headers.py YOUR_ADDUCTED_FILE_PATH YOUR_RAW_RESULTS_DIR_POS/ positive # or
 python scripts/populate_headers.py YOUR_ADDUCTED_FILE_PATH YOUR_RAW_RESULTS_DIR_NEG/ negative 
```

## Generating the final spectral file

We concatenate each documented mgf files to a single spectral mgf file:
```
bash scripts/concat.sh YOUR_RAW_RESULTS_DIR_POS/ isdb_pos.mgf # or
bash scripts/concat.sh YOUR_RAW_RESULTS_DIR_NEG/ isdb_neg.mgf
```

# Results

Can be found at:

//TODO Zenodo (2022-07-07 AR)
