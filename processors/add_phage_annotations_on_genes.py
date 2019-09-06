import sys
from collections import defaultdict
from functools import partial

phage_file = "/home/hielke/bep/res/phages.coords.csv"
gene_file = "/home/hielke/bep/res/gene_classification.sel.90"

out_file = sys.stdout
out = partial(print, file=out_file, sep="\t")
# As also in processors/combinevirsorter.awk and processors/pipeline.bash
distance_to_prophage = 500

virsorter = defaultdict(lambda: defaultdict(set))
for line in map(lambda l: l.strip(), open(phage_file)):
    genome, contig, phage = line.split(maxsplit=2)
    virsorter[genome][contig].add(phage)

for F in map(lambda l: l.strip().split('\t'), open(gene_file)):
    genome_id, gene_id, start_gene, end_gene, *_ = F
    contig_id = "_".join(gene_id.split("_")[:-1])
    for phage in virsorter[genome_id][contig_id]:
        phage_id, start_phage, end_phage, *_ = phage.split()
        if start_phage == "NA":
            break
        start_gene, end_gene, start_phage, end_phage = map(
            int, [start_gene, end_gene, start_phage, end_phage]
        )
        if (start_gene <= end_phage + distance_to_prophage) and (
            end_gene >= start_phage - distance_to_prophage
        ):
            break
    else:
        phage_id = "0"
    out(*F, phage_id)
