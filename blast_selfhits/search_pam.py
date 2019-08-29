#!/usr/bin/env python3

import argparse
import csv
import fileinput
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from functools import partial

from Bio.Seq import Seq
from Bio.SeqUtils import nt_search

ap = argparse.ArgumentParser()
ap.add_argument(
    "--pam",
    nargs="?",
    default="/home/hielke/new_analysis/pam.xml",
    type=argparse.FileType(),
)
ap.add_argument(
    "-o", "--out", nargs="?", default=sys.stdout, type=argparse.FileType("a+")
)
ap.add_argument("inputfiles", nargs="*")
ap.add_argument("-t", "--sep", default="\t")

args = ap.parse_args()

out = partial(print, file=args.out, sep=args.sep)

pam = defaultdict(list)
side_dict = dict()
root = ET.parse(args.pam).getroot()
for child in root:
    if not child or child.attrib.get("unknown") or not child.attrib.get("side"):
        continue
    side_dict[child.attrib["name"]] = child.attrib["side"]
    pam[child.attrib["name"]] = [seq.text.strip() for seq in child]
args.pam.close()


class MatchTarget:
    def __init__(self, five_prime, three_prime, side):
        self.side = side
        # # test code
        # if side == "3":
        #     self.target = Seq(three_prime)
        # else:
        #     self.target = Seq(five_prime)
        # Since the blast hits directly at the spacer, we have to convert to
        # the other side (the protospacer).
        if side == "3":
            self.target = Seq(five_prime).reverse_complement()
        else:
            self.target = Seq(three_prime).reverse_complement()

    def get_target(self, size):
        if self.side == "3":
            return self.target[:size]
        else:
            return self.target[-size:]

    def match_target(self, pam):
        match = nt_search(str(self.get_target(len(pam))), pam)
        return len(match) > 1


"""
-1 : no PAM available
0  : PAM incorrect
1  : crispr_type correct
2  : genome_type correct
3  : both correct
"""
for F in csv.reader(fileinput.input(args.inputfiles), delimiter=args.sep):
    five_prime = F[8]
    three_prime = F[9]
    crispr_type = F[-2]
    genome_type = F[-1]

    s_c = side_dict.get(crispr_type)
    s_g = side_dict.get(genome_type)
    if not s_c and not s_g:
        out(*F, -1, "NULL", "NULL")
        continue
    res = 0
    if s_g:
        mt = MatchTarget(five_prime, three_prime, s_c)
        poss_pam = mt.get_target(10)
        side = s_g
        for p in pam[genome_type]:
            if mt.match_target(p):
                res += 2
                break
    if s_c:
        mt = MatchTarget(five_prime, three_prime, s_c)
        poss_pam = mt.get_target(10)
        side = s_c
        for p in pam[crispr_type]:
            if mt.match_target(p):
                res += 1
                break
    out(*F, res, poss_pam, side)
