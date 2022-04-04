#!/bin/sh
INPUT=$1
OUTPUT=$2
find $INPUT -type f -name '*.mgf' | while read F; do cat ${F} >> $OUTPUT; done
