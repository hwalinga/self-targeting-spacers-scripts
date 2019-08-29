#!/usr/bin/env python3

import sys
import argparse
import fileinput
from functools import partial
import os
import re
import ast
from operator import itemgetter
from uuid import uuid4
from collections import Counter, defaultdict
from BCBio import GFF

ap = argparse.ArgumentParser()
ap.add_argument('inputfiles', nargs='*',
                help='macsyfinder.out files')
ap.add_argument('--gene-coords', action='store_true')
ap.add_argument('-t', '--sep', default='\t')
ap.add_argument('--prodigal', default=None, nargs='?',
                const="/hosts/linuxhome/mgx/DB/PATRIC/prodigal/data/")
args = ap.parse_args()

log = partial(print, file=sys.stderr, sep=args.sep)
out = partial(print, file=sys.stdout, sep=args.sep)

output_ext = '_crisprcasfinder'


class ProdiagalPrinter:
    """ Save jobs so that they can be combined with the prodigal output in
    one go"""
    prodigal_path = args.prodigal

    def __init__(self):
        self.genome = ""
        self.jobs = []
        self.out = out

    def init(self):
        self.genome = ""
        self.jobs.clear()

    def save_job(self, *job):
        if not self.jobs:
            self.genome = job[0]
        if self.prodigal_path and self.genome != job[0]:
            log("[ERROR]", self.genome, job[0])
            return
        self.jobs.append(list(job))

    def match_genes(self):
        prodigal_file = os.path.join(self.prodigal_path, self.genome + ".gff")
        contig_dict = defaultdict(list)
        for job in self.jobs:
            contigs = ["_".join(contig_id.split("_")[:-1]) for contig_id
                       in job[3]]
            most_common_contig = Counter(contigs).most_common(1)[0][0]
            job[3] = [contig_id for contig_id, contig
                      in zip(job[3], contigs)
                      if contig == most_common_contig]
            job[4] = "|".join(job[4])
            contig_dict[most_common_contig].append(job)
        for rec in GFF.parse(prodigal_file):
            for job in contig_dict[rec.id]:
                start = sys.maxsize
                end = 0
                for contig_id in job[3]:
                    loc = rec.features[
                        int(contig_id.split("_")[-1]) - 1].location
                    start = min(start, int(loc.start))
                    end = max(end, int(loc.end))
                job[3:4] = [start, end]
                job[1:1] = [rec.id]

    def print_jobs(self):
        if not self.jobs:
            return
        if self.prodigal_path:
            self.match_genes()
        for job in self.jobs:
            self.out(*job)
        self.init()


if args.prodigal:
    p = ProdiagalPrinter()
    out = p.save_job


def analyse_clusters(f, disambiguation=False):
    """Analyse clusters of macsyfinder.out.

    Can also analyse the subsection of the Disambiguation step.
    In that case disambiguation links is None.
    That is also the reason why this is in a function,
    so that it can be called again when a disambiguation is found.
    (not recursive though)
    """
    if disambiguation:
        clusters = dict()
    for line in f:
        if line.startswith("MacSyFinder"):
            genome = re.search(r'/([^/]*).faa', next(f)).group(1)
            # genome = os.path.basename(os.path.dirname(
            #     fileinput.filename())).rstrip(output_ext)
            clusters = dict()
            disambiguation_links = dict()
            children_clusters = set()
            if args.prodigal:
                p.print_jobs()
            continue
        # This is for the disambiguation step.
        if line.startswith("=> none kept"):
            return None, False
        if line.startswith('--- Cluster'):
            # There is an optional (indicated by the regex modifier "?")
            # question mark. (Which is the same symbol!)
            crispr_type = [re.search(r'--- Cluster\s(.*)\s\??\s---',
                                     line).group(1)]
            genes_pos = ast.literal_eval(next(f))
            genes_type = ast.literal_eval(next(f))
            genes_number = frozenset(ast.literal_eval(next(f)))
            line = next(f)
            considering = re.search(r'Considering (.*):', line)
            crispr_condition = ["undefined"]
            ambiguous = None
            # None if nothing is considered, either disambiguous or undefined
            # True if after considering no decision is made
            # False if there is decision made
            while considering:
                crispr_condition.append(re.search(r'->\s(.*)(?:\s|$)',
                                                  next(f)).group(1))
                line = next(f)
                if line.startswith("===="):  # decision made
                    crispr_type.append(considering.group(1))
                    ambiguous = False
                    break
                if line.startswith('='):  # Putative
                    crispr_type.append(considering.group(1))
                    ambiguous = True
                    break

                crispr_type.append(considering.group(1))
                considering = re.search(r'Considering (.*):', line)

            if ambiguous is False:  # decision made
                clusters[genes_number] = (crispr_type[-1],
                                          crispr_condition[-1], genes_pos,
                                          genes_type)
            if ambiguous is True:
                crispr_types = set(crispr_type[1:])
                clusters[genes_number] = ("|".join(crispr_types),
                                          "ambiguous" if len(crispr_types) > 1
                                          else crispr_condition[-1],
                                          genes_pos, genes_type)

            if line.startswith("Disambiguation"):
                if disambiguation:
                    # This means disambiguation inception. Odd
                    log("Disambiguation inception:", genome)
                    continue
                # Do the same for all subsets, but link to parent:
                clusters[genes_number] = (crispr_type[0], crispr_condition[0],
                                          genes_pos, genes_type)
                children, flag = analyse_clusters(f, disambiguation=True)
                if children is None:
                    # Happens when "none" of the other options are viable
                    continue
                if flag:
                    # Unexpected next line
                    log("Odd disambiguation case (unexpected next line):",
                        genome)
                if not children:
                    # Somehow there were no disambiguations made.
                    log("No disambiguation found:", genome)
                    continue
                disambiguation_links[genes_number] = set(children.keys())
                children_clusters.update(children.keys())
                clusters.update(children)
                continue
            if disambiguation:
                # Make sure when this is the case we terminate early
                if line.startswith("Store:"):
                    clusters[genes_number] = (crispr_type, crispr_condition,
                                              genes_pos, genes_type)
                    continue
                if line.startswith("=> disambiguated"):
                    clusters[genes_number] = (crispr_type, crispr_condition,
                                              genes_pos, genes_type)
                    return clusters, False
                return clusters, True
            if ambiguous is None:  # Stopped early i.e. Undefined
                clusters[genes_number] = (crispr_type[0], crispr_condition[0],
                                          genes_pos, genes_type)

        # End of file is near indicator.
        # This does not always happen, but when it does not happen there is
        # no output anyway.
        if line.startswith("Counter("):
            for genes_num, cluster in clusters.items():
                if genes_num in children_clusters:
                    continue

                disambiguation_clusters_nums = \
                    disambiguation_links.get(genes_num)

                if not disambiguation_clusters_nums:
                    out(genome, *cluster)
                    continue

                disambiguation_clusters = {k: v for k, v
                                           in clusters.items()
                                           if k in
                                           disambiguation_clusters_nums}
                if len(disambiguation_clusters) == 1:
                    # Now just combine with parent.
                    children_cluster = disambiguation_clusters.popitem()[1]
                    cluster_types = {cluster[0], children_cluster[0]}
                    out(genome, "|".join(cluster_types),
                        "ambiguous" if len(cluster_types) > 1
                        else children_cluster[1],
                        cluster[2], cluster[3])
                    continue
                printed = set()
                total_printed = 0
                print_jobs = []
                for disambiguation_cluster_num, disambiguation_cluster \
                        in disambiguation_clusters.items():
                    if disambiguation_cluster_num in printed:
                        continue
                    # intersection test of any other one.
                    # if not intersection, just print this one.
                    # do not combine with parent, since it was splitted.
                    has_intersection = {disambiguation_cluster_num}
                    for other in disambiguation_clusters.keys():
                        if disambiguation_cluster_num == other:
                            continue
                        # intersection
                        if disambiguation_cluster_num & other:
                            has_intersection.add(other)

                    if len(has_intersection) == 1:
                        print_jobs.append([genome,
                                           *disambiguation_cluster])
                        printed.update(has_intersection)
                        total_printed += 1
                    else:
                        intersection_clusters = [
                            disambiguation_clusters[num] for
                            num in has_intersection]
                        crispr_type = "|".join(map(
                            itemgetter(0), intersection_clusters))
                        genes_pos = set().union(map(set, map(
                            itemgetter(2), intersection_clusters)))
                        genes_type = set().union(map(set, map(
                            itemgetter(3), intersection_clusters)))
                        print_jobs.append([genome, crispr_type,
                                           "ambiguous", genes_pos,
                                           genes_type])
                        printed.update(has_intersection)
                        total_printed += 1
                if total_printed == 1:  # Still combine with parent.
                    cluster_types = {cluster[0], print_jobs[0][1]}
                    out(genome, "|".join(cluster_types),
                        "ambiguous" if len(cluster_types) > 1
                        else print_jobs[0][2],
                        cluster[2], cluster[3])
                else:
                    for job in print_jobs:
                        out(*job)


ID = uuid4()
log("PROGRAM STARTED", ID)
with fileinput.input(args.inputfiles) as f:
    analyse_clusters(f)
if args.prodigal:
    p.print_jobs()
log("PROGRAM FINISHED", ID)
