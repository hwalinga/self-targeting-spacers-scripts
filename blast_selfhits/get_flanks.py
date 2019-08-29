#!/usr/bin/env python3

import fileinput
import argparse
import sys
import os
import csv
from collections import defaultdict
from glob import iglob
from multiprocessing import Pool, Lock
from functools import partial

from Bio import SeqIO

ap = argparse.ArgumentParser(description="print 5 print and 3 print resp")
ap.add_argument('selfhits', nargs='*')
ap.add_argument('-o', '--out', nargs='?', type=argparse.FileType('w'),
                default=sys.stdout)
ap.add_argument('-g', '--genome-path',
                default='/hosts/linuxhome/mgx/DB/PATRIC/patricdb-201*/')
ap.add_argument('-t', '--field-seperator', nargs='?', default='\t')
ap.add_argument('-P', '--max-procs', nargs='?', type=int, default=1, const=0)
args = ap.parse_args()

#lock = Lock()
#def synced_out(*inp):
#    with lock:
#        args.out.write(args.field_seperator.join(inp))
synced_out = partial(
        print, file=args.out, sep=args.field_seperator, flush=True)


flank = 100

# {genome: contig: {spacer_id, begin, start, etc.}}
selfhits = defaultdict(lambda: defaultdict(list))
for F in csv.reader(fileinput.input(args.selfhits),
                    delimiter=args.field_seperator):
    selfhits[F[0]][F[2]].append((F[1], *F[3:]))


def process_selfhit_item(selfhit_item, out=synced_out):
    k, v = selfhit_item
    genome_file = next(iglob(os.path.join(args.genome_path, k + '.fna')))
    for record in SeqIO.parse(genome_file, 'fasta'):
        for hit in v[record.id]:
            seq = record.seq.upper()
            # First assume Forward and adjust for Reverse.
            if hit[1] < hit[2]:
                direction = 'Forward'
                begin, end = int(hit[1]), int(hit[2])
            else:
                direction = 'Reverse'
                begin, end = int(hit[2]), int(hit[1])
            # Carefull: Blast range is 1-indexed and end inclusive,
            # but python range is 0-indexed and end exclusive.
            # This is changed with begin -= 1.

            begin_record = begin - flank - 1
            begin_record = 0 if begin_record < 0 else begin_record
            five_prime = seq[begin_record: begin - 1]

            end_record = end + flank
            end_record = len(seq) if end_record > len(seq) \
                else end_record
            three_prime = seq[end: end + flank]
            if direction == 'Reverse':
                five_prime, three_prime = three_prime.reverse_complement(), \
                                          five_prime.reverse_complement()
            out(k, record.id, *hit, five_prime, three_prime)


if args.max_procs == 1:
    for selfhit in selfhits.items():
        process_selfhit_item(selfhit, out=partial(print, file=args.out,
                             sep=args.field_seperator))
else:
    if args.max_procs == 0:
        args.max_procs = None  # This means maximum available processes
    with Pool(args.max_procs) as p:
        any(p.imap_unordered(process_selfhit_item, selfhits.items()))
