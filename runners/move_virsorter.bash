#!/bin/bash

SRC="./output/"
DST="/hosts/linuxhome/mgx/DB/PATRIC/patric_virsorter_new/data/"
CUR_DIR=$PWD

for NAME in $SRC*; do 
    DIR=${NAME##*/}
    GENOME=${DIR%*.virsorter}
    cd $NAME
    cp VIRSorter_global-phage-signal.csv "$DST""$GENOME"_VIRSorter_global-phage-signal.csv
    cp ./fasta/VIRSorter_mga_final.predict "$DST""$GENOME"_VIRSorter_mga_final.predict
    for i in ./Predicted_viral_sequences/*fasta; do 
        file=${i##*/}
        cp "$i" "$DST""$file"_"$GENOME"
    done;
    cd "$CUR_DIR"
done;
