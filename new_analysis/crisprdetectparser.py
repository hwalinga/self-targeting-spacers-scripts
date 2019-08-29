#!/usr/bin/env python3

import argparse
import fileinput
from itertools import takewhile
import os
import sys
from functools import partial
from collections import defaultdict
import re

ap = argparse.ArgumentParser()
ap.add_argument('files', nargs='*')
ap.add_argument('-s', '--spacers', nargs='?', const='spacers_dir',
                default=None, help='Create spacers folder')
ap.add_argument('-o', '--out', nargs='?', type=argparse.FileType('a+'),
                default=sys.stdout)
ap.add_argument('-t', '--sep', nargs='?', default='\t')
args = ap.parse_args()

print = partial(print, file=args.out)
spacers_ext = ".spacers.fna"
filename2genome = lambda f: f.split('/')[-1].split(".crisprdetect")[0]

if args.spacers:
    if not os.path.isdir(args.spacers):
        os.mkdir(args.spacers)
    spacers_file = None

with fileinput.input(args.files) as f:
    for line in f:
        if f.isfirstline():
            genome = filename2genome(f.filename())
            contigs_counter = defaultdict(int)
            if args.spacers:
                if spacers_file:
                    spacers_file.close()
                spacers_file = open(
                    os.path.join(args.spacers, genome + spacers_ext), 'w')

        # Array definition starts now.
        if line.startswith(">"):
            F = line.split()
            orientation = F[-1].strip()
            contig = F[0].lstrip(">")
            contigs_counter[contig] += 1
            number = contigs_counter[contig]
            array_id = '{}_{}'.format(contig, number)

            all(takewhile(lambda l: not l.startswith("==="), f))
            # Loop over all spacers.
            for ind, line in enumerate(f):
                F = line.split()
                spacer_sequence = F[5]
                if spacer_sequence == "|":
                    end = int(F[0]) + (1 if orientation == 'Forward'
                                       else -1) * int(F[1])
                    num_spacers = ind
                    next(f)
                    line = next(f)
                    repeat = line.split()[4]
                    break
                if ind == 0:
                    begin = F[0]
                if args.spacers:
                    print('>{}_{}'.format(
                        array_id, ind + 1), file=spacers_file)
                    print(F[5], file=spacers_file)

            # Continue loop till Array family
            for line in f:
                if line.startswith("# Questionable array :"):
                    score = line.split(": ")[-1].strip()

                if line.startswith("# Array family :"):
                    family = re.search(
                        r'# Array family :\s(.*?)(?:\s|$)', line).group(1)
                    if orientation != 'Forward':
                        begin, end = end, begin
                    print(genome, contig, array_id, begin, end, orientation,
                          num_spacers, family, score, repeat, sep=args.sep)
                    break
