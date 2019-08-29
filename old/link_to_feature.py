#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 27 22:18:18 2018

@author: hielke
"""

from collections import defaultdict
import sys

n = 0
test = True

if test:
    feature_file = "/home/hielke/bep/data/genomes/907.4.PATRIC.features.tab"
    gene_hits_file = "/home/hielke/bep/py/misc/gene_hits_907.4"
else:
    feature_file = sys.argv[1]
    gene_hits_file = sys.argv[2]

feature_dict = defaultdict(dict)
with open(feature_file) as f:
    next(f)  # header
    for line in f:
        splitted_line = line.split("\t")
        if splitted_line[4] == "source":
            continue
        contig = splitted_line[2]
        start, end = map(int, splitted_line[9:11])
        direction = splitted_line[11]
        start_gene = start if direction == "+" else end
        product = splitted_line[14]
        feature_dict[contig][start_gene] = splitted_line[14]


for line in open(gene_hits_file):
    splitted_orf = line.split()[1].split("_")
    contig = "_".join(splitted_orf[:-2])
    start_gene = min(map(int, splitted_orf[-1][1:-1].split("-")))
    feature = feature_dict[contig].get(start_gene)
    feature = feature if feature else ""
    print("\t".join([line.rstrip("\n"), feature]))
