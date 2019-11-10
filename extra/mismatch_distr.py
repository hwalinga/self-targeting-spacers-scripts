import argparse
import csv
import re
import sys
from collections import defaultdict

from Bio import SeqIO, pairwise2

ap = argparse.ArgumentParser()
ap.add_argument("inputfiles", nargs="*")
ap.add_argument("-t", "--sep", default="\t")
ap.add_argument("--spacers", type=argparse.FileType("r"))
args = ap.parse_args()

min_id = 90
out = sys.stdout
testing = False

# genome -> contig -> tuple
spacer_dict = defaultdict(lambda: defaultdict(list))

for F in csv.reader(args.spacers, delimiter=args.sep):
    if F[0] == "genome_id":  # header
        continue
    if float(F[8]) < min_id:
        continue
    genome = F[0]
    contig = F[4]
    array_id = F[13]
    start = int(F[2])
    end = int(F[3])
    spacer_seq = F[1]

    # from 1-indexed including end to 0-index excluding end
    if start < end:
        reverse = False
        coords = (start - 1, end)
    else:
        reverse = True
        coords = (end - 1, start)

    spacer_dict[genome][contig].append((coords, spacer_seq, array_id, reverse))

# if testing:
#     for genome, hits in spacer_dict.items():
#         print(genome, hits)

genome_mismatches = defaultdict(list)

for file in args.inputfiles:
    genome = re.sub(r"^.*\/", "", file).rstrip(".fna")
    for record in SeqIO.parse(file, "fasta"):
        contig = record.id
        genome_contig_seq = record.seq.upper()
        for coords, spacer_seq, array_id, reverse in spacer_dict[genome][contig]:
            hit_seq = genome_contig_seq[slice(*coords)]
            if reverse:
                hit_seq = hit_seq.reverse_complement()
            alignment = pairwise2.align.globalms(hit_seq, spacer_seq, 1, -1, -1, -1)
            if testing:
                score = alignment[0][2]
                hit_ident = score / len(spacer_seq)
                print(spacer_seq)
                print(hit_seq)
                print(genome, hit_ident, sep=args.sep, file=out)
                print(pairwise2.format_alignment(*alignment[0]))
                assert hit_ident > min_id / 100
            alignment_string = pairwise2.format_alignment(*alignment[0]).split("\n")[1]
            mismatches = [
                ind for ind, char in enumerate(alignment_string) if char != "|"
            ]
            mismatches_normalized = [ind / len(spacer_seq) for ind in mismatches]
            for m in mismatches_normalized:
                print(genome, array_id, m, sep=args.sep, file=out)
