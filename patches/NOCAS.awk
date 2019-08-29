#/usr/bin/awk -f
# USAGE
# awk -f NOCAS.awk genome.types.tsv gene.matched.tsv.3 > patched
FNR==NR{if($4==0){nocas[$1]=1;} next;}
$2 in nocas{$17 = "NOCAS"; print $0; next;}
1
