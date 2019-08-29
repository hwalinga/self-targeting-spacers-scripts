#!/usr/bin/env python3

import argparse
import csv
import fileinput
import sys
from collections import Counter
from functools import partial
from itertools import chain, groupby, tee
from operator import itemgetter

ap = argparse.ArgumentParser()
ap.add_argument("inputfiles", nargs="*")
ap.add_argument(
    "-o", "--out", nargs="?", type=argparse.FileType("w"), default=sys.stdout
)
ap.add_argument("-t", "--sep", nargs="?", default="\t")
args = ap.parse_args()

out = partial(print, file=args.out, sep=args.sep)
log = partial(print, file=sys.stderr, sep=args.sep)

rows = csv.reader(fileinput.input(args.inputfiles), delimiter="\t")
for gene, assignments in groupby(rows, key=itemgetter(0)):

    desc, interpro, goterms = tee(filter(lambda row: len(row) >= 14, assignments), 3)
    desc = filter(bool, map(itemgetter(5), desc))
    interpro = filter(bool, map(itemgetter(12), interpro))
    goterms = filter(bool, map(itemgetter(13), goterms))
    desc = list(desc)

    c = Counter(chain(desc, interpro))
    try:
        m = map(itemgetter(0), next(groupby(c.most_common(), key=itemgetter(1)))[1])
    except StopIteration:
        m = [""]
    # If there are multiple description with the same count:
    # take longest string.
    d = max(m, key=len)
    if d:  # If there isn't any, take go-terms
        out(gene, d)
        continue
    item = Counter(goterms).most_common(1)
    if item:
        out(gene, item[0])
    else:
        out(gene, "UNKNOWN")
