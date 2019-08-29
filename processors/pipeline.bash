#!/bin/bash
NUMCORES=35
AWK=mawk
alias tawk="$AWK -F '\t' -v OFS='\t'"

find $CRISPRDETECTFILES -type f -name '*.crisprdetect' -print0 |
    parallel --xargs -q -0 -P $NUMCORES --line-buffer python3 \
    parsecrisprdetect/crisprdetectparser.py --spacers spacers \
    > spacers.info.tsv
find $CRISPRDETECTFILES -type f -name '*.crisprdetect.fp' -print0 |
    parallel --xargs -q -0 -P $NUMCORES --line-buffer python3 \
    parsecrisprdetect/crisprdetectparser.py --spacers deg_spacers \
    > deg_spacers.info.tsv

find /linuxhome/tmp/hielke/macsyfinder/output/ -mindepth 1 -maxdepth 1 -print0 |
    sed -z 's@$@/macsyfinder/macsyfinder.out@' |
    parallel --xargs -q --linebuffer -0 -P 35 \
    python3 ~/new_analysis/parsemacsyfinder.py \
    2> parsemacsyfinder.log |
    mtawk '{print $1"@"$2,$0}' | sort -k1,1 |
    join -t $'\t' -a{1,2} -e "NULL" -o 0 1.{2..8} 2.{2..11} - \
    <( mtawk '{print $1"@"$2,$0}' spacers.info.tsv | sort -k1,1) |
    python3 createtyping.py # genome.types / crispr.types

join -t $'\t' -a{1,2} -e "NULL" -o 0 1.{4..7} 2.{4..6} 2.9 combinetypes/types.tsv.keyed combinetypes/spacers.info.tsv.keyed | xcl

mtawk '{print $1"@"$2,$0}' ../parsemacsyfinder_SIXTH/types.tsv | sort -k1,1 --parallel=30 | join -a{1,2} -e "NULL" -t $'\t' -o 0 1.{4..7} 2.{4..6} 2.9 1.8 - <( mtawk '{print $1"@"$2,$0}' ../parsecrisprdetect/spacers.info.tsv | sort -k1,1 --parallel=30 ) | head


find ../blastspacers/results/ -type f -print0 | parallel --xargs -q -0 -P 70 --line-buffer mawk -F '\t' -v OFS='\t' -f ~/new_analysis/identfilter.awk |
    python3 ~/blast_selfhits/filter_out_arrays.py --arrayfiles ../parsecrisprdetect/{deg_,}spacers.info.tsv |
    python3 ~/blast_selfhits/get_flanks.py -P 70 |
    mawk -F '\t' -v OFS='\t' '{array=$3; sub("_[^_]*$", "", array); print $1"@"array,$0}' |
    sort --parallel=70 -k1,1 | join -t $'\t' - <(
    mtawk '{print $1"@"$3,$0}' ../parsecrisprdetect/spacers.info.tsv | sort -k1,1 --parallel=70 | join -t $'\t' - <( mtawk '{print $1"@"$3,$0}' ../combinetypes/crispr.types.tsv | sort -k1,1 --parallel=70 ) | sort -k2,2 --parallel=70 | join -t $'\t' -1 2 -2 1 - <(sort -k1,1 --parallel=70 ../combinetypes/genome.types.tsv) | mtawk '{print $2,$1,$3,$4,$8,$10,length($11),$13,$16}' | sort -k1,1 --parallel=70 ) | mtawk '{print $2,$3,$5,$6,$4,$7,$8,$9,$10,$11,$14,$15,$16,$17,$18,$19}' > flanked.merged


parallel --pipepart -a flanked.merged python3 ~/blast_selfhits/align_flanks.py >flanked.out 2>flanked.log
mtawk '$8 > 2000' flanked.out > flanked.contigsize.out
parallel --pipepart -a flanked.contigsize.out python3 ~/blast_selfhits/search_pam.py pam.out log.out


blastspacers/runblastspacers.bash
find ../blastspacers/results/ -type f -print0 |
    parallel --xargs -q -0 -P $NUMCORES \
    tawk -f identfilter.awk |
    python3 ~/blast_selfhits/filter_out_arrays.py \
    --arrayfiles parsecrisprdetect/{deg_,}spacers.info.tsv |
    python3 /home/hielke/blast_selfhits/get_flanks.py -P $NUMCORES |
    tawk
    '{array=$3; sub("_[^_]*$", "", array); print $1"@"array $0}' | sort -k1,1 |
    join -t $'\t' - <(tawk '{print $1"@"$3 $0}' spacers.info.tsv | sort -k1,1) |
    python3 ~/blast_selfhits/align_flanks.py | test

# Join CRISPR and phages.

distance_to_prophage=500
join -t $'\t' -a 1 -e NULL <( mtawk '{print $1"@"$2,$0}' ../parseselfhits_SEC/pam.out | sort -k1,1 --parallel=70 ) \
    <( mtawk '{print $1"@"$2,$0}' phages.coords.tsv | sort -k1,1 --parallel=70 ) |
    mtawk -v d=$distance_to_prophage -f ~/new_analysis/combinevirsorter.awk |
    python3 ~/new_analysis/print_unique.py

find data -type f -print0 | parallel -0 -P 30 python3 parse_interproscan.py {} ">" parsed/{/.}.features
find data -type f | parallel -P 70 --xargs python3 ~/new_analysis/parse_prodigal.py > gene_classification

join -t $'\t' <( mtawk '{print $1"@"$2,$0;}' ../combine_virsorter/virsorter.unique.combined.tsv |
    sort -k1,1 -S30% --parallel=70 )     <( mtawk '{gene=$2; sub("_[^_]*$","",gene); print $1"@"gene,$0;}' gene_classification |
    sort -k1,1 -S30% --parallel=70 ) |     python3 ~/new_analysis/print_gene_match.py > gene.matched.tsv.more
