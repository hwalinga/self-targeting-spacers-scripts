export GENOMES_DIR="/hosts/linuxhome/mgx/tmp/PATRIC/patricdb-201*/"
GENOMES_FILE="all_genomes_listed"

NUM_CORES=20
run_blast () {
    echo $1
    cp $GENOMES_DIR$1 . 
    blastx -query $1 -db faa-anticrispr2/all_anti_crispr.faa -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -evalue "1e-6" -out $1.blast
    rm $1
}
export -f run_blast

cat $GENOMES_FILE | xargs -d '\n' -I '{}' -P $NUM_CORES bash -c 'run_blast "{}"'
