#!/bin/bash
RAWFILES=$1
[[ $RAWFILES =~ /disco([^/]*)/ ]] && RANGE=${BASH_REMATCH[1]}
python3 classify_based_on_disco.py $RAWFILES $RANGE

rm -rf ./genomes
    
