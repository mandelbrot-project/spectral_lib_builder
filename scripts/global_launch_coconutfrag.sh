#!/bin/bash

for f in /home/allardp/coconut/sub_test/*
do
     echo "Launching bash Number $f"
        sbatch "$f"
done