#!/bin/bash

NUM_CORES=35
GENOMES_FILE="test"
CRISPRCASFINDER_FOLDER="/linuxhome/tmp/hielke/CRISPRCasFinder/"
export PATH="$PATH:$CRISPRCASFINDER_FOLDER/bin"
export MACSY_HOME="$CRISPRCASFINDER_FOLDER/macsyfinder-1.0.5/"
export CRISPRCASFINDER="$CRISPRCASFINDER_FOLDER/CRISPRCasFinder.pl"
export SOFILE="$CRISPRCASFINDER_FOLDER/sel392v2.so"
export CUR_DIR=$PWD
export OUTPUT=$CUR_DIR/"output6"
export GENOMES=$CUR_DIR/"genomes"
export GENOMES_DIR="/hosts/linuxhome/mgx/DB/PATRIC/patricdb-201*/"
export EXTENSION="_crisprcasfinder"

mkdir -p "$OUTPUT"
mkdir -p "$GENOMES"

runcrisprcasfinder() {

    echo $1
    cp -f $GENOMES_DIR$1.fna $GENOMES
    mkdir $1"_tmp_crisprcasfinder" && cd $_
    perl $CRISPRCASFINDER -q -cas -i $GENOMES/$1.fna \
            -so $SOFILE -out $OUTPUT/$1$EXTENSION -keep -log
    cd "$CUR_DIR"
    rm -f $GENOMES/$1.fna
    rm -rf $OUTPUT/$1$EXTENSION/CasFinder
    rm -rf $1"_tmp_crisprcasfinder"

}
export -f runcrisprcasfinder

xargs -a $GENOMES_FILE -d '\n' -I '{}' -P $NUM_CORES bash -c 'runcrisprcasfinder "{}"' 
