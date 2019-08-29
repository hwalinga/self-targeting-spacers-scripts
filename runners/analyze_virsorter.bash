#!/bin/bash
NUM_CORES=24
GENOMES_FILE="test_list"
export HOME=/home/ksenia
export GENOMES_DIR="/hosts/linuxhome/mgx/DB/PATRIC/patricdb-201*/"
export VIRSORTER="/linuxhome/tmp/hielke/franklin/bin/wrapper_phage_contigs_sorter_iPlant.pl"
export CUR_DIR=$PWD
export OUTPUT=$CUR_DIR"/output/"
export GENOMES_DIR_LOCAL=$CUR_DIR"/genomes/"

mkdir -p $OUTPUT
mkdir -p $GENOMES_DIR_LOCAL
run_virsorter() {

    echo $1
    cp -f $GENOMES_DIR$1.fna $GENOMES_DIR_LOCAL
    virsorter_dir=$OUTPUT$1.virsorter
    mkdir $virsorter_dir && cd $_
    yes "" | perl $VIRSORTER -ncpu 1 -f $GENOMES_DIR_LOCAL$1.fna
    cd $CUR_DIR
    rm $GENOMES_DIR_LOCAL$1.fna

}
export -f run_virsorter
cat $GENOMES_FILE | xargs -d '\n' -I '{}' -P $NUM_CORES bash -c 'run_virsorter "{}"'
