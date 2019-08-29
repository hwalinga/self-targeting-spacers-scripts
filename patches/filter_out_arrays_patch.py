#!/usr/bin/env python3

import argparse
import csv
import fileinput
import sys
from collections import defaultdict
from functools import partial

ap = argparse.ArgumentParser()
ap.add_argument("-f", "--arrayfiles", nargs="+")
ap.add_argument(
    "-o", "--out", nargs="?", default=sys.stdout, type=argparse.FileType("a+")
)
ap.add_argument("inputfiles", nargs="*")
ap.add_argument("--flank", default=50)
ap.add_argument("-t", "--sep", default="\t")

args = ap.parse_args()

print = partial(print, file=args.out, sep=args.sep)

# {(genome_id, contig_id): [(begin, end), ...]}
arrays = defaultdict(list)
for F in csv.reader(fileinput.input(args.arrayfiles), delimiter="\t"):
    arrays[F[0], F[1]].append((int(F[3]), int(F[4])))

# Find any overlapping range or flanked by.

for F in csv.reader(fileinput.input(args.inputfiles), delimiter="\t"):
    if fileinput.isfirstline():  # There is now a header
        print(*F)
        continue

    if F[2] < F[3]:
        begin, end = int(F[2]), int(F[3])
    else:
        begin, end = int(F[3]), int(F[2])

    size = end - begin + 1
    for array in arrays[F[0], F[4]]:
        if size + array[1] - array[0] + 1 + args.flank > max(
            array[1] - begin + 1, end - array[0] + 1
        ):
            break  # This means overlapping array, or closer than args.flank
    else:  # No break, so no overlapping array
        print(*F)
