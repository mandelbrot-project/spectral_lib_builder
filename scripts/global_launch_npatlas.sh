#!/bin/bash

for f in /home/allardp/bash_files/npatlas_bash/*
do
     echo "Launching bash Number $f"
        sbatch "$f"
done