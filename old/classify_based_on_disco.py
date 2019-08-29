#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 10:18:11 2018

@author: hielke
"""

from glob import glob
from collections import defaultdict
from operator import itemgetter
from subprocess import Popen, PIPE
from itertools import chain
from types import MethodType
from shutil import copy
from os import mkdir
from os import path
import re
import sys

from cluster import create_clusters_from_proximity_matrix, gene_proximity_matrix


def get_key_from_list_dict(list_item, list_dict):
    """If you have a dictionary with lists, this function will return the key 
    of the list that contains the `list_item`. It'll return the first it finds.
    """
    for key, alist in list_dict.items():
        if list_item in alist:
            return key


# /hosts/linuxhome/plasmid/tmp/hielke

UNDERSCORES_IN_GENOME = 0
FULL_REFS = "/home/hielke/CRISPRdisco/data/fullrefs/"
RAW_FILES = sys.argv[1]
RANGE = sys.argv[2]
GENOME_DIR = "/hosts/linuxhome/mgx/tmp/PATRIC/patricdb-201*/"
LOG_FILE = "log"
TARGET_FILE = "reclassified_crispr_systems_" + RANGE
target_file = open(TARGET_FILE, "w")
test = False
if not test:
    if not path.exists("./genomes"):
        mkdir("./genomes/")
    genomes = [
        line.split(",")[1].replace("-", ".")
        for line in open(glob(RAW_FILES + "*summary_list*")[0])
    ]
else:
    test_genome = "/home/hielke/bep/data/genomes/108619.109.fna"
    genomes = [0]

log = open(LOG_FILE, "w")

# ---- READ IN THE GENOMES ----


gene_locations = defaultdict(dict)
prog = re.compile("^>(.*?)_(\d+)\s\[(\d+)\D+(\d+)\]")
# >FMGY01000_001_1 [22 - 66] ^= FMGY01000_001 1 22 66

for genome in genomes:
    if not test:
        genome_location = GENOME_DIR + genome + ".fna"
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

    genome = genome.replace(".", "-")

    cmd = (
        'getorf -find 1 %s /dev/stdout | grep "^>" | cut -d" " -f1,2,3,4'
        % genome_location
    )
    res = Popen(cmd, shell=True, stdout=PIPE, universal_newlines=True)
    for line in res.stdout:
        m = prog.search(line)
        contig, gene_id, s, e = (
            m.group(1),
            int(m.group(2)),
            int(m.group(3)),
            int(m.group(4)),
        )
        entry = gene_locations[genome].get(contig)
        if not entry:
            # because gene_id is 1-indexed we initialize 0 to None.
            entry = [None]
            gene_locations[genome][contig] = entry
        if not gene_id == len(entry):
            log.write("Gene_id is incorrect here:\t(%s,%s)\n" % (genome, contig))
        entry.append((s, e))

# ---- PROCESS THE PROTEIN DATA

protein_crispr_type_dict = {}  # contains 'type gene'


def get_gene_group(file):
    file_name = file.split("/")[-1]
    gene_group = file_name.split("_")[0]
    if gene_group != "Ref":
        return gene_group
    return file_name.split("_")[5]


unidentified_map_Burstein = {"CasY": "TypeV-D", "CasX": "TypeV-E"}

for file in glob(FULL_REFS + "*"):
    for line in open(file):
        if not line.startswith(">"):
            continue
        anid = line.split()[0][1:]
        gene_group = get_gene_group(file)

        if "|" not in line:

            if "cas1_" in file:
                if "Cas4" in line:  # fusion
                    protein_crispr_type_dict[anid] = "CAS-I-U,TypeV-B*", gene_group
                    continue
                protein_crispr_type_dict[anid] = "PROT_CAS_1", gene_group
                continue

            items = [
                item for key, item in unidentified_map_Burstein.items() if key in line
            ]

            if not items:  # possibly Shmakov
                if "cas12c" in file:
                    protein_crispr_type_dict[anid] = "TypeV-C", gene_group
                    continue
                else:
                    print("ERROR: COULD NOT IDENTIFY")
                    print("    " + file)
                    print("    " + line)
                    continue

            protein_crispr_type_dict[anid] = items[0], gene_group

            assert len(items) == 1, "More items matched than expected"

            continue

        if "cas1_" in file and not "Makarova" in file:
            protein_crispr_type_dict[anid] = "PROT_CAS_1", gene_group
            continue

        if "Shmakov" in file:
            family = line.split("|")[-1].strip()
            if family.startswith(">"):
                print(file)
                print(line)
            protein_crispr_type_dict[anid] = family, gene_group
            continue

        family = line.split("|")[-3].strip()
        if family.startswith(">"):
            family = "PROT_CAS_2"
        if family == "cas12a":
            family = "TypeV-A"
        protein_crispr_type_dict[anid] = family, gene_group

        if False:
            print("ERROR: COULD NOT IDENTIFY")
            print("    " + file)
            print("    " + line)

for key, item in protein_crispr_type_dict.items():
    if item[0] == "TypeVB":
        protein_crispr_type_dict[key] = "TypeV-B", item[1]
    if item[0] == "CAS-V":  # cas12a(cpf1)
        protein_crispr_type_dict[key] = "TypeV-A", item[1]


# ---- PROCESS THE HITS ----


import csv


genome_crispr_hits = defaultdict(lambda: defaultdict(dict))

blast_results_file = glob(RAW_FILES + "*blast_results*")
# testfile = "/home/hielke/bep/data/disco/discodata/disco801,1400p/raw_files/csm5_blast_results_Seqs_2018-07-09.csv"
for file in blast_results_file:
    with open(file) as f:
        csvreader = csv.reader(f, delimiter=",", quotechar='"')
        next(csvreader)  # header
        for row in csvreader:
            if int(row[3].split(".")[0]) < 75 or int(row[-2].split(".")[0]) < 75:
                continue
            gene = row[1]
            hit_id = row[2].split("_")
            genome = hit_id[0 : UNDERSCORES_IN_GENOME + 1]
            genome = "_".join(genome)
            gene_id = int(hit_id[-1])
            contig = "_".join(hit_id[UNDERSCORES_IN_GENOME + 1 : -1])

            genome_crispr_hits[genome][contig][gene_id] = protein_crispr_type_dict[gene]

# genome_crispr_hits[gene] = genome_crispr_hits_gene

available_genes = {
    "2OG",
    "CasX",
    "CasY",
    "DEDDh",
    "DinG",
    "PD-DExK",
    "PrimPol",
    "c2c10",
    "c2c4",
    "c2c5",
    "c2c8",
    "c2c9",
    "cas1",
    "cas10",
    "cas10d",
    "cas11",
    "cas12a",
    "cas12b",
    "cas12c",
    "cas13a",
    "cas13b",
    "cas13c",
    "cas2",
    "cas3",
    "cas4",
    "cas5",
    "cas6",
    "cas7",
    "cas8a",
    "cas8b",
    "cas8c",
    "cas8e",
    "cas8f",
    "cas8u",
    "cas9",
    "casR",
    "cmr1gr7",
    "cmr3gr5",
    "cmr4gr7",
    "cmr5gr11",
    "cmr6gr7",
    "cmr7",
    "cmr8gr7",
    "csa3",
    "csa5gr11",
    "csb1gr7",
    "csb2gr5",
    "csb3",
    "csc1gr5",
    "csc2gr7",
    "cse2gr11",
    "csf1gr8",
    "csf2gr7",
    "csf3gr5",
    "csf4gr11",
    "csf5gr6",
    "csm2gr11",
    "csm3gr7",
    "csm4gr5",
    "csm5gr7",
    "csm6",
    "csn2",
    "csx1",
    "csx10gr5",
    "csx15",
    "csx16",
    "csx18",
    "csx19",
    "csx20",
    "csx21",
    "csx22",
    "csx23",
    "csx24",
    "csx25",
    "csx26",
    "csx3",
}


# signature_genes = {}


subtypes_dict = {
    "CAS-I": {
        "CAS-I-A",
        "CAS-I-B",
        "CAS-I-C",
        "CAS-I-D",
        "CAS-I-E",
        "CAS-I-F",
        "CAS-I-U",
    },
    "CAS-II": {"CAS-II-A", "CAS-II-B", "CAS-II-C"},
    "CAS-III": {"CAS-III-A", "CAS-III-B", "CAS-III-C", "CAS-III-D"},
    "CAS-IV": {"CAS-IV"},
    "CAS-V": {"TypeV-A", "TypeV-B", "TypeV-C", "TypeV-D", "TypeV-E"},
    "CAS-VI": {"Type6A", "Type6B", "Type6C"},
    "CAS-VU": {"TypeVU1", "TypeVU2", "TypeVU3", "TypeVU4", "TypeVU5"},
}

class_dict = {
    "class-1": {"CAS-I", "CAS-III", "CAS-IV"},
    "class-2": {"CAS-V", "CAS-VI", "CAS-VU", "CAS-II"},
}

translate_type_dict = {
    "CAS-I": "TypeI",
    "CAS-II": "TypeII",
    "CAS-III": "TypeIII",
    "CAS-IV": "TypeIV",
    "CAS-V": "TypeV",
    "CAS-VI": "TypeVI",
    "CAS-VU": "TypeV-U",
    "CAS-I-A": "TypeI-A",
    "CAS-I-B": "TypeI-B",
    "CAS-I-C": "TypeI-C",
    "CAS-I-D": "TypeI-D",
    "CAS-I-E": "TypeI-E",
    "CAS-I-F": "TypeI-F",
    "CAS-I-U": "TypeI-U",
    "CAS-II-A": "TypeII-A",
    "CAS-II-B": "TypeII-B",
    "CAS-II-C": "TypeII-C",
    "CAS-III-A": "TypeIII-A",
    "CAS-III-B": "TypeIII-B",
    "CAS-III-C": "TypeIII-C",
    "CAS-III-D": "TypeIII-D",
    "TypeV-A": "TypeV-A",
    "TypeV-B": "TypeV-B",
    "TypeV-C": "TypeV-C",
    "TypeV-D": "TypeV-D",
    "TypeV-E": "TypeV-E",
    "Type6A": "TypeVI-A",
    "Type6B": "TypeVI-B",
    "Type6C": "TypeVI-C",
    "TypeVU1": "TypeV-U1",
    "TypeVU2": "TypeV-U2",
    "TypeVU3": "TypeV-U3",
    "TypeVU4": "TypeV-U4",
    "TypeVU5": "TypeV-U5",
    "PROT_CAS_1": "",
    "PROT_CAS_2": "",
}


# if gene == "Ref_set_from_Makarova_2015_cse2gr11_generated_2017-05-15.fasta":
# return "cs2gr11"


VU_systems = {"TypeV-U1", "TypeV-U2", "TypeV-U3", "TypeV-U4", "TypeV-U5"}


def process_type(atype):
    res = {}
    if atype in {"PROT_CAS_1", "PROT_CAS_2"}:
        return {"subtype": "", "thetype": "", "theclass": "class-1"}
    res["subtype"] = translate_type_dict[atype]
    thetype = [key for key, item in subtypes_dict.items() if atype in item][0]
    res["type"] = translate_type_dict[thetype]
    res["theclass"] = [key for key, item in class_dict.items() if thetype in item][0]
    return res


class Classification:
    def __init__(self):
        self.type = set()
        self.subtype = set()
        self.theclass = set()
        self.VU = set()
        self.others = set()
        self.completeness = ""
        self.ambiguous = None
        self.contig = ""
        self.start = 0
        self.end = sys.maxsize
        self.genes = set()
        self.adaption_module = None

    def update_attr(self, attr_dict):
        for key, item in attr_dict.items():
            if key in dir(self):
                to_update_item = getattr(self, key)
                if type(to_update_item) is set:
                    to_update_item.add(item)

    def set_to_string(self):
        for key in dir(self):
            if key.startswith("__"):
                continue
            to_change_item = getattr(self, key)
            if type(to_change_item) is not set:
                continue
            to_change_item = filter(bool, to_change_item)
            setattr(self, key, ",".join(sorted(to_change_item)))

    def completeness_not_single(self):
        assert self.genes, "There are no genes added yet!"
        self.has_adaption()
        theclass = "class-1" if "class-1" in self.theclass else "class-2"
        if theclass == "class-1":
            return self.completeness_complex_eff_module()
        if "TypeV-B" not in self.subtype:
            self.completeness = "Complete"
            return
        # That one case that a CRISPR loci only has the cas1/cas4 fusion.
        self.completeness = "Complete" if "cas12b" in self.genes else "Partial"

    def completeness_complex_eff_module(self):
        assert self.genes, "There are no genes added yet!"
        gr5 = "cas5" in self.genes or any(["gr5" in g for g in self.genes])
        gr7 = "cas7" in self.genes or any(["gr7" in g for g in self.genes])
        hass_eff = gr5 and gr7
        if "TypeIII" in self.type:
            gr11 = "cas11" in self.genes or any(["gr11" in g for g in self.genes])
            gr10 = "cas10" in self.genes or "cas10d" in self.genes
            hass_eff = hass_eff and gr11 and gr10
        else:
            if "TypeI-F" in self.subtype or "TypeI-E" in self.subtype:
                hass_eff = hass_eff and (
                    "cas6" in self.genes or any(["gr6" in g for g in self.genes])
                )
                if "TypeI-E" in self.subtype:
                    hass_eff = hass_eff and any(["gr11" in g for g in self.genes])
            if not {"Type-IV"} == self.subtype:
                if not {"TypeI-D"} == self.subtype:
                    hass_eff = hass_eff and (
                        any(["cas8" in g or "gr8" in g for g in self.genes])
                    )
                else:
                    hass_eff = hass_eff and (
                        "cas10" in self.genes or "cas10d" in self.genes
                    )
        self.completeness = "Complete" if hass_eff else "Partial"

    def has_adaption(self):
        assert self.genes, "There are no genes added yet!"
        if (
            "TypeI" in self.type
            or "TypeII" in self.type
            or "TypeIII" in self.type
            or "TypeV" in self.type
        ):
            self.adaption_module = "cas1" in self.genes
        else:
            self.adaption_module = True

    def process_others(self, sorted_types, score):

        while sorted_types:
            atype, other_score = sorted_types.pop()
            if other_score == score:
                self.ambiguous = True
                self.update_attr(process_type(atype))
            else:
                translated_type = translate_type_dict[atype]
                if translated_type in self.subtype:
                    continue
                if translated_type in VU_systems:
                    self.VU.add(translated_type)
                else:
                    self.others.add(translated_type)

    def remove_VU(self):
        self.VU.update(self.subtype & VU_systems)
        self.subtype -= VU_systems
        self.type -= {"TypeV-U"}


# Classification = namedtuple('Classification',
#      'type subtype VU others completeness ambiguous external spacer_array')


def classify(types_and_genes):
    # returns {
    #   type: [RES(, RES, ...)? | None],
    #   subtype: [RES(, RES, ...)? | None],
    #   VU: [RES(, RES, ...)? | None],
    #   Others: [RES(, RES, ..)? | None]
    #   Completeness: [Complete | Partial | Single],
    #   Ambiguous: [True | False]
    #   external: [1-6] #(Virsorter),
    #   spacer_array: str
    #   genes: (gene(, gene, ...)?)
    #   }
    print(">>>>>>>")
    print(types_and_genes)

    num_genes = len(set(types_and_genes))
    candidate_types = []
    for type_and_gene in types_and_genes:
        types = type_and_gene[0]
        for atype in types.split(","):
            if atype in subtypes_dict:
                for subtype in subtypes_dict[atype]:
                    candidate_types.append(subtype)
            else:
                candidate_types.append(atype)

    count_types = defaultdict(int)
    for atype in candidate_types:
        if atype[-1] == "*":
            score = 1
            atype = atype[:-1]
        else:
            score = 2
        count_types[atype] += score
    sorted_types = sorted(count_types.items(), key=itemgetter(1))

    classification = Classification()
    genes = [type_and_gene[1] for type_and_gene in types_and_genes]
    classification.genes.update(genes)
    atype, score = sorted_types.pop()
    if num_genes == 1:
        print("1")
        classification.ambiguous = False
        classification.completeness = "Single"
        classification.update_attr(process_type(atype))
        classification.remove_VU()
    elif score >= num_genes * 2 * 2 / 3:
        print("2")
        classification.ambiguous = False
        classification.update_attr(process_type(atype))
        classification.remove_VU()
        classification.process_others(sorted_types, score)
        classification.completeness_not_single()
    else:
        print("3")
        classification.ambiguous = True
        classification.update_attr(process_type(atype))
        classification.remove_VU()
        classification.process_others(sorted_types, score)
        classification.completeness_not_single()

    return classification


# classifications = defaultdict(list)


def key_is_attr(key, inst):
    return not key.startswith("__") and not type(getattr(inst, key)) == MethodType


def write_to_file(file, aclassifcation, genome):
    aclassification.set_to_string()

    line = "\t".join(
        map(
            str,
            [
                (
                    lambda attr: str(
                        (0, 1)[attr]
                        if type(attr) is bool
                        else -1
                        if attr is None
                        else attr
                    )
                )(getattr(aclassification, key))
                for key in dir(aclassification)
                if key_is_attr(key, aclassification)
            ],
        )
    )
    target_file.write("%s\t%s\n" % (genome, line))


dummy_classification = Classification()
line = "\t".join(
    [key for key in dir(dummy_classification) if key_is_attr(key, dummy_classification)]
)
target_file.write("%s\t%s\n" % ("genome", line))

for genome, contig_hits in genome_crispr_hits.items():
    if test and genome != "108619-109":
        continue
    print("continue with the loop")
    print(contig_hits.items())
    for contig, hits in contig_hits.items():
        gene_ids = list(hits.keys())
        gene_coords = gene_locations[genome][contig]
        gene_coords_collected = [gene_coords[gene_id] for gene_id in gene_ids]
        proximity_matrix = gene_proximity_matrix(gene_coords_collected, 3000)
        clusters = create_clusters_from_proximity_matrix(proximity_matrix)
        for cluster in clusters:
            types_and_genes = [hits[gene_ids[c]] for c in cluster]
            aclassification = classify(types_and_genes)
            aclassification.contig = contig
            gene_coords_this_cluster = list(
                chain(*[gene_coords[gene_ids[c]] for c in cluster])
            )
            aclassification.start = min(gene_coords_this_cluster)
            aclassification.end = max(gene_coords_this_cluster)
            write_to_file(target_file, aclassification, genome)


log.close()
target_file.close()
