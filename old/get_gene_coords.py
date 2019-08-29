#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 14:59:37 2018

@author: hielke
"""

from collections import defaultdict
import re
from shutil import copy 
from subprocess import Popen, PIPE
from glob import glob
from os import path
import os 
import pathlib
import sys
from operator import itemgetter
from functools import reduce

MAIN_GENOME_DIR = ""

genome_oi = "907.4"
#genome_oi = "108619.109"

test = True
if test:
    MAIN_GENOME_DIR = None
    GENOME_DIR = "/home/hielke/bep/data/genomes/"
    ORF_DIR = "/home/hielke/bep/data/genomes/orf/"
    LOCAL_GENOME_DIR = "./genomes"
    OUTPUT = ""
    INPUT_FILE = "/home/hielke/bep/jups/fp_combined_with_extra.tsv"
else:
    MAIN_GENOME_DIR = "FILL THIS ONE IN"
    GENOME_DIR = None
    ORF_DIR = "./orf/"
    LOCAL_GENOME_DIR = "./genomes/"
    OUTPUT = "./output/"
    INPUT_FILE = ""

for directory in [ORF_DIR, LOCAL_GENOME_DIR]:
    pathlib.Path(directory).mkdir(exist_ok=True) 

genomes = [None]
genomes = []
test_genome = GENOME_DIR + genome_oi + ".fna"
log = open('log', 'w')

gene_locations = defaultdict(dict)
prog = re.compile('^>(.*?)_(\d+)\s\[(\d+)\D+(\d+)\]') 
# >FMGY01000_001_1 [22 - 66] ^= FMGY01000_001 1 22 66

for genome in genomes or []:
    if not test:
        genome_location = GENOME_DIR + genome_oi + ".fna"
        genome_location = glob(genome_location)
        if len(genome_location) == 0:
            log.write("Not found:\t%s\n" % genome)
            continue
        genome_location = genome_location[0]
        if " " in genome_location:
            log.write("Has spaces!:\t%s\n" % genome)
            continue
        
        copy(genome_location, "./genomes/")
        genome_location = "./genomes/" + genome + ".fna"
    else:
        genome_location = test_genome
        genome = re.search(".*/(.*?).fna", test_genome).group(1)
#        genome_location = glob(genome_location)
    
    genome = genome_oi.replace('.', '-')
    
    orf_file = ORF_DIR + genome + ".orf"
    if not path.isfile(orf_file):
        print("creating orf file", file=sys.stderr)
        cmd = 'getorf -find 3 %s %s' % (genome_location, orf_file)
        res = Popen(cmd, shell=True, stdout=PIPE, universal_newlines=True)
        for line in res.stdout: print(line)
        
    
    for line in open(orf_file):
        m = prog.search(line)
        if not m:
            continue
        contig, gene_id, s, e = m.group(1), int(m.group(2)), \
                int(m.group(3)), int(m.group(4))
        entry = gene_locations[genome].get(contig)
        if not entry:
            # because gene_id is 1-indexed we initialize 0 to None.
            entry = [(0, 0)] 
            gene_locations[genome][contig] = entry
        if not gene_id == len(entry):
            print("Gene_id is incorrect here:\t(%s,%s)\n" 
                % (genome, contig))
        entry.append((s, e))
    
        
# COUNTING IN EACH BUCKET
        
gene_locations_bucket = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
#with open("/home/hielke/bep/jups/hits/" + genome_oi + ".hits.tsv") as f:
#    header = next(f)
#    for line in f:
#        splitted_line = line.split()
#        spacer = splitted_line[1]
#        g = int(splitted_line[3])
#        contig = splitted_line[4]
#        contig_coords = gene_locations[genome][contig]
#        for ind, (s, e) in enumerate(contig_coords):
#            if s <= g <= e or e <= g <= s:
#                gene_locations_bucket[genome][contig][ind].add(spacer)

def get_gene_locations(genome):
    if test:
        genome_location = path.join(GENOME_DIR, genome + ".fna")
    else:
        source_genome = glob(path.join(MAIN_GENOME_DIR, genome + ".fna"))
        genome_location = copy(source_genome[0], LOCAL_GENOME_DIR)
        
    
    if test:
        orf_file = ORF_DIR + genome + ".orf"
    if not path.isfile(orf_file):
        print("creating orf file")
        cmd = 'getorf -find 3 %s %s' % (genome_location, orf_file)
        res = Popen(cmd, shell=True, stdout=PIPE, universal_newlines=True)
        for line in res.stdout: print(line)
        
    gene_locations = dict()
    for line in open(orf_file):
        m = prog.search(line)
        if not m:
            continue
        contig, gene_id, s, e = m.group(1), int(m.group(2)), \
                int(m.group(3)), int(m.group(4))
        entry = gene_locations.get(contig)
        if not entry:
            # because gene_id is 1-indexed we initialize 0 to None.
            entry = [(0, 0)] 
            gene_locations[contig] = entry
        if not gene_id == len(entry):
            print("Gene_id is incorrect here:\t(%s,%s)\n" 
                % (genome, contig))
        entry.append((s, e))
    
    if not test:
        os.remove(orf_file)
        os.remove(genome_location)
    
    return gene_locations

def print_genome_results(genome_oi, gene_locations):
    for contig, contigs_dict in gene_locations_bucket[genome_oi].items():
        for ind, val in contigs_dict.items():
            if len(val) > 0:
                spacers = "|".join(map(itemgetter(0), val))
                inframe = "|".join(map(reduce(lambda f, g: lambda x: g(f(x)),
                                              [itemgetter(1), int, str]), val))
                coords = "-".join(map(str, gene_locations[contig][ind]))
                if test:
                    output_file = sys.stdout
                else:
                    output_file = open(path.join(OUTPUT, genome_oi + ".gene.hits"))
                print(genome_oi + '\t' + contig + '_' + str(ind)  
                        + '_(' + coords + ')\t' + spacers + '\t' + inframe, 
                        file=output_file)

printed = list()
current_genome = ''
with open(INPUT_FILE) as f:
    for line in f:
        splitted_line = line.split('\t')
        genome = splitted_line[0]
        if genome != current_genome:
            if not splitted_line[6]: # no selfhits
                continue
            if test and genome != genome_oi:
                continue
            if current_genome:
                # print current_genome if it exists
                printed.append(current_genome)
                print_genome_results(current_genome, gene_locations)
            gene_locations = get_gene_locations(genome)
            current_genome = genome
        contig = splitted_line[9]
        spacer_start, spacer_end = map(lambda x: int(x.split('.')[0]), 
                                       splitted_line[6:8])
        spacer = splitted_line[18]
        contig_coords = gene_locations[contig]
        
        for ind, (s, e) in enumerate(contig_coords):
            minimum, maximum = min(s, e), max(s, e)
            if minimum <= spacer_start <= maximum \
                    or minimum <= spacer_start <= maximum:
                inframe = (s - e) * (spacer_start - spacer_end) > 0
                gene_locations_bucket[genome][contig][ind].add((spacer, inframe))
   
if current_genome not in printed:
    print_genome_results(current_genome, gene_locations)      

log.close()
