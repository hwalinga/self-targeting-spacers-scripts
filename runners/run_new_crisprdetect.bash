NUM_CORES=22
export CRISPRDETECT="/home/hielke/CRISPRDetect_2.2_1/CRISPRDetect.pl"
export GENOMES_FILE="interesting_genomes.rest"
export GENOMES_DIR="/hosts/linuxhome/mgx/tmp/PATRIC/patricdb-201*/"
run_crisprdetect () {
    echo $1
    cp $GENOMES_DIR$1.fna . 
    perl $CRISPRDETECT -f $1.fna -o $1.crisprdetect -array_quality_score_cutoff 3 -T 10
    rm $1.fna
}
export -f run_crisprdetect

cat $GENOMES_FILE | xargs -d '\n' -I '{}' -P $NUM_CORES bash -c 'run_crisprdetect "{}"'
