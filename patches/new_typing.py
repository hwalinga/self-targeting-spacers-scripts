#!/usr/bin/env python3

import argparse
import sys

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("inputfile", nargs="?")
ap.add_argument("-f", "--typefile")
ap.add_argument(
    "-o", "--out", nargs="?", default=sys.stdout, type=argparse.FileType("a+")
)

args = ap.parse_args()

inp = args.inputfile or sys.stdin

df = pd.read_csv(
    inp,
    sep="\t",
    index_col=False,
    keep_default_na=False,
    dtype={"genome_id": str, "spacer_pos": int, "spacer_size": int},
).set_index("genome_id")

df_genome_class = pd.read_csv(
    args.typefile,
    sep="\t",
    names=["genome_id", "type", "class", "num_genes"],
    dtype={"genome_id": str, "num_genes": int},
).set_index("genome_id")

# Naming used for the groupings.
# Also see $PATCHES/new_typing.py
NOCAS = "No CAS genes"
Mixed = "Mixed CRISPR System"
Complete = "Complete CRISPR System"
Incomplete = "Incomplete CRISPR System"
# 'other'  # Those that exist in all, but not anymore in the selfhit group.
other = "Cas-TypeVC | CAS-TypeVIB2"
# 'CAS-TypeVC': 1, 'CAS-TypeVIB2': 2


def new_type_group(row):
    if row["num_genes"] == 0:
        return NOCAS
    if "Unclassified" in row["class"]:
        return Incomplete
    if "MultipleConfirmed" in row["class"]:
        # Confirmed Mixed:
        # The CRISPR array associated confirms that the type is a mixed type.
        return Mixed
    if "/" in row["type"] or row["type"] == "CAS":
        return Mixed
    return row["type"]


#     if "Confirmed" in row['class']:
#         return "Confirmed " + row['type']
#     if "Accompanied" in row['class']:
#         return "Accompanied " + row['type']


df_genome_class["new_type_group"] = df_genome_class.apply(new_type_group, axis=1)
df = df.join(df_genome_class)

df.to_csv(sys.stdout, sep="\t")
