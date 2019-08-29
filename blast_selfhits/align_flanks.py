#!/usr/bin/env python3

import argparse
from Bio import pairwise2
from collections import defaultdict
import sys
from functools import partial
import os
import csv
import fileinput

ap = argparse.ArgumentParser()
ap.add_argument('inputfiles', nargs='*')
ap.add_argument('-o', '--out', nargs='?', type=argparse.FileType('w'),
                default=sys.stdout)
ap.add_argument('-t', '--sep', nargs='?', default='\t')
args = ap.parse_args()

out = partial(print, file=args.out, sep=args.sep)
log = partial(print, file=sys.stderr, sep=args.sep)

for F in csv.reader(fileinput.input(args.inputfiles), delimiter='\t'):
    repeat_len = int(F[13])
    if len(F[8]) < repeat_len or len(F[9]) < repeat_len:
        log(*F[0:4], '[WARN]', "possibly too short for correct alignment")
    if not F[8] or not F[9]:
        log(*F[0:4], '[SKIP]', "at the side")
        out(*F)
        continue
    score = int(pairwise2.align.globalms(
        F[8][-repeat_len:], F[9][:repeat_len], 1, -1, -1, -1)[0][2])

#    print(pairwise2.format_alignment(*pairwise2.align.globalms(F[8][-repeat_len:],
#                                                 F[9][:repeat_len])[0]))
    # out("%.2f" % (score / repeat_len))
    if score / repeat_len < 0.75:
        out(*F)
