#!/usr/bin/env python3

import fileinput
import sys
import argparse
from collections import defaultdict
import csv
from functools import partial

ap = argparse.ArgumentParser()
ap.add_argument('-f', '--arrayfiles', nargs='+')
ap.add_argument('-o', '--out', nargs='?', default=sys.stdout,
                type=argparse.FileType('a+'))
ap.add_argument('inputfiles', nargs='*')
ap.add_argument('-t', '--sep', default='\t')

args = ap.parse_args()

print = partial(print, file=args.out, sep=args.sep)

# {(genome_id, contig_id): [(begin, end), ...]}
arrays = defaultdict(list)
for F in csv.reader(fileinput.input(args.arrayfiles), delimiter='\t'):
    arrays[F[0], F[1]].append((int(F[3]), int(F[4])))

# Find any overlapping range.

for F in csv.reader(fileinput.input(args.inputfiles), delimiter='\t'):
    if F[3] < F[4]:
        begin, end = int(F[3]), int(F[4])
    else:
        begin, end = int(F[4]), int(F[3])
    size = end - begin + 1
    for array in arrays[F[0], F[2]]:
        if size + array[1] - array[0] + 1 > max(array[1] - begin + 1,
                                                end - array[0] + 1):
            break  # This means overlapping array
    else:  # No break, so no overlapping array
        print(*F)
