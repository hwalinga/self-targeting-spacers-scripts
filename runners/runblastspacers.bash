#!/bin/bash

NUM_CORES=35
SPACERS_FOLDER="/linuxhome/tmp/hielke/parsecrisprdetect/spacers/"
export PATRIC="/hosts/linuxhome/mgx/DB/PATRIC/patricdb-201*"
export CURDIR="$PWD"
export RESULTS="$CURDIR/results"
export TMP="$CURDIR/tmp"
export EXTENSION=".spacers.fna"

mkdir -p "$RESULTS"
mkdir -p "$TMP"

blastspacers() {

    spacerspath=$1
    spacerfile="${spacerspath##*/}"
    genome_id="${spacerfile%$EXTENSION}"
    TMP_DIR="$TMP/${genome_id}_selfhits"
    mkdir $TMP_DIR && cd $_
    cp ${PATRIC}/${genome_id}.fna . 
    makeblastdb -in "${genome_id}.fna" -dbtype nucl
    blastn -query $spacerspath -db $TMP_DIR/$genome_id.fna -task blastn-short -out $RESULTS/$genome_id.selfhits -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -dust no -gapopen 10 -gapextend 10 -penalty -1 -evalue 1
    rm -f *
    cd $CURDIR
    rmdir $TMP_DIR

}

export -f blastspacers

find $SPACERS_FOLDER -type f -print0 | 
        xargs -0 -I '{}' -P $NUM_CORES bash -c 'blastspacers "{}"'
