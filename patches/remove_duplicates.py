#!/usr/bin/env python3
"""
* Remove all genomes that have the same hit patterns
(thus indicating a double genome)
* Duplicates resulting from duplicates contigs
(same spacer_id with same coords on different contig.)
* Double arrays (same seq, same hit, same array location)
Also add information about spacer position.
"""

import sys
from itertools import islice
from operator import itemgetter

import pandas as pd

inp = sys.argv[1] if len(sys.argv) > 1 else sys.stdin

df = pd.read_csv(
    inp,
    sep="\t",
    names=[
        "genome_key",
        "genome_id",
        "contig_hit",
        "c_1",
        "c_2",
        "spacer_id",
        "hit_ident",
        "spacer_size",
        "contig_size",
        "3_prime",
        "5_prime",
        "array_id",
        "array_size",
        "array_confid",
        "repeat_size",
        "array_type",
        "genome_type",
        "PAM_confid",
        "poss_PAM",
        "PAM_side",
        "phage_class",
        "with_phage",
        "gene_hit",
        "gene_id",
        "gene_class",
    ],
    index_col=False,
    keep_default_na=False,
    dtype={"genome_id": str},
).set_index("genome_id")

double_genome = df.groupby("spacer_id").apply(
    lambda x: "-".join((islice(set(x.index), 0, None)))
)
double_array_ids = double_genome[double_genome.str.contains("-")].unique()
remove_these_genomes = set()
for i in double_array_ids:
    genomes = sorted(i.split("-"))
    genome_hits = [len(df.loc[genome_id]) for genome_id in genomes]
    max_index, _ = max(enumerate(genomes), key=itemgetter(1))
    del genomes[max_index]
    remove_these_genomes.update(genomes)
df = df[~df.index.isin(remove_these_genomes)]

# By the definition of the spacer ID,
# the last number indicates the position in the array
# The orientation of the array is defined by CRISPRDetect

df = df.assign(spacer_pos=lambda x: x["spacer_id"].str.split("_").str[-1].astype(int))


def seq_iter(filehandle):
    for seq_id in filehandle:
        yield seq_id.strip()[1:], next(filehandle).strip()


spacers = seq_iter(open("all_spacers.fna"))
all_seq = pd.DataFrame(spacers, columns=["spacer_id", "seq"])

# Merge main data frame with spacer sequences.
df = (
    df.reset_index()
    .merge(all_seq, how="left", on="spacer_id", validate="m:1")
    .set_index("genome_id")
)

df = (
    df.reset_index()
    .groupby(["spacer_id", "c_1", "c_2"])
    .first()
    .reset_index()
    .set_index("genome_id")
)
df = (
    df.reset_index()
    .groupby(["seq", "c_1", "c_2", "contig_hit", "spacer_pos"])
    .first()
    .reset_index()
    .set_index("genome_id")
)

# relative spacer position
# to change a number n range a-b with fixed a and variable b,
# to a relative position 0-100 use
# (n-a)/(b-a)
# 0: leader
# 1: tail
df["rel_spacer_pos"] = df["spacer_pos"].subtract(1).divide(df["array_size"].subtract(1))

df["name_spacer_pos"] = df["rel_spacer_pos"].apply(
    lambda x: {0: "leader", 1: "tail"}.get(x) or "middle"
)

df.to_csv(sys.stdout, sep="\t")
