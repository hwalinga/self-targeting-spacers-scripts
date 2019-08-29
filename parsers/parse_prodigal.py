#!/usr/bin/env python3

import argparse
import os
import sys
from functools import partial

from BCBio import GFF

ap = argparse.ArgumentParser()
ap.add_argument("-t", "--sep", default="\t")
ap.add_argument(
    "--prodigal",
    nargs="?",
    default=None,
    const="/hosts/linuxhome/mgx/DB/PATRIC/prodigal/data/",
)
ap.add_argument(
    "--interproscan-parsed",
    default="/hosts/linuxhome/mgx/DB/PATRIC/interproscan/parsed/",
)
ap.add_argument("inputfiles", nargs="*")
args = ap.parse_args()

log = partial(print, file=sys.stderr, sep=args.sep)
out = partial(print, file=sys.stdout, sep=args.sep)


def parse_parsed_interproscan(genome_id):
    """Return dictionary (gene_id -> desc)
    if no interproscan file can be find return None"""
    file = os.path.join(args.interproscan_parsed, genome_id + ".faa.features")
    if not os.path.isfile(file):
        return None
    return {
        item[0].strip(): item[1].strip()
        for item in map(lambda line: line.split("\t"), open(file))
    }


prodigal_files = os.listdir(args.prodigal) if args.prodigal else args.inputfiles
for prodigal_file in map(lambda l: l.strip(), prodigal_files):
    if prodigal_file.split(".")[-1] == "faa":
        continue
    genome_id = ".".join(os.path.basename(prodigal_file).split(".")[:-1])
    gene_id_dict = parse_parsed_interproscan(genome_id)
    for rec in GFF.parse(prodigal_file):
        contig_id = rec.id
        for ind, gene in enumerate(rec.features):
            gene_id = contig_id + "_" + str(ind + 1)
            if gene_id_dict is None:
                gene_classification = "UNDEFINED"
            else:
                gene_class = gene_id_dict.get(gene_id)
                gene_classification = gene_class if gene_class else "UNCLASSIFIED"
            out(
                genome_id,
                gene_id,
                gene.location.start,
                gene.location.end,
                gene.strand,
                gene_classification,
            )
