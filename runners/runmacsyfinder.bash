#!/bin/bash

NUM_CORES=35
GENOMES_FILE="test"
CRISPRCASFINDER_FOLDER="/linuxhome/tmp/hielke/CRISPRCasFinder"
export PATH="$PATH:$CRISPRCASFINDER_FOLDER/bin"
export MACSY_HOME="$CRISPRCASFINDER_FOLDER/macsyfinder-1.0.5/"
export CUR_DIR=$PWD
export OUTPUT=$CUR_DIR/"output"
export GENOMES=$CUR_DIR/"genomes"
export PATRIC="/hosts/linuxhome/mgx/DB/PATRIC"
export GENOMES_DIR="$PATRIC/patricdb-201*/"
export EXTENSION="_crisprcasfinder"

mkdir -p "$OUTPUT"
mkdir -p "$GENOMES"

runmacsyfinder() {

    echo $1
    cp $GENOMES_DIR$1.fna $GENOMES
    mkdir $1$EXTENSTION && cd $_
    mkdir prodigal
    mkdir macsyfinder
    prodigal -i $GENOMES_DIR$1.fna a -c -m -g 11 -f gff -q -p single \
            -a prodigal/$1.faa -o prodigal/$1.gff 
    macsyfinder -w 1 all --out-dir macsyfinder --db-type ordered_replicon \
            -d CRISPRCASFINDER_FOLDER/CasFinder-2.0.2/DEF-SubTyping-2.0.2 \
            -p CRISPRCASFINDER_FOLDER/CasFinder-2.0.2/CASprofiles-2.0.2 \
            --sequence-db prodigal/$1.faa
    mv  prodigal/{$1.faa,$1.gff} $PATRIC/prodigal/data
    cd "$CUR_DIR"
    
    rm -f $GENOMES/$1.fna

}
export -f runmacsyfinder

xargs -a $GENOMES_FILE -d '\n' -I '{}' -P $NUM_CORES \
        bash -c 'runmacsyfinder "{}"' 
