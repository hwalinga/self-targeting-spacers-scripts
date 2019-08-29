import argparse
import fileinput
import sys
from functools import partial
import csv
from uuid import uuid4
import re
from operator import itemgetter
from collections import defaultdict
from itertools import chain

ap = argparse.ArgumentParser()
ap.add_argument('inputfiles', nargs='*')
ap.add_argument('-t', '--sep', default='\t')
ap.add_argument('--found-repeats', type=argparse.FileType('r'),
                default="/home/hielke/new_analysis/found_repeat_families")
ap.add_argument('--genome', type=argparse.FileType('a+'),
                default='genome.types.tsv')
ap.add_argument('--crispr', type=argparse.FileType('a+'),
                default='crispr.types.tsv')
args = ap.parse_args()


def crisprrepeat2castype(crisprrepeat):
    if "NA" in crisprrepeat:
        return frozenset()
    majorType, subType = crisprrepeat.strip().split('-')
    return frozenset("CAS-Type" + majorType + sub
                     for sub in subType.split('/'))


classified_repeats = {item for line in args.found_repeats
                      for item in crisprrepeat2castype(line)}

genome_out = partial(print, file=args.genome, sep=args.sep)
crispr_out = partial(print, file=args.crispr, sep=args.sep)
log = partial(print, file=sys.stderr, sep='\n', end='\n---- ----\n')

ID = uuid4()
log("[INFO] Program %s started" % ID)

array2crispr_distance = 20000


def split_list_oncondition(condition, alist):
    false_list, true_list = [], []
    for i in alist:
        (false_list, true_list)[condition(i)].append(i)
    return false_list, true_list


def megasplit(pattern, string):
    """Split with zero lenght
    https://stackoverflow.com/questions/29988595/python-regex-splitting-on-pattern-match-that-is-an-empty-string
    """
    splits = list((m.start(), m.end()) for m in re.finditer(pattern, string))
    starts = [0] + [i[1] for i in splits]
    ends = [i[0] for i in splits] + [len(string)]
    return [string[start:end] for start, end in zip(starts, ends)]


casclassarray_possibilities = {
    1: "Loner",  # The crispr array does not have a cas locus near.
    2: "Wrong",  # The crispr array has a wrong family compared to the close locus.
    3: "Undefined",  # The crispr cas locus has a type with undefined repeat family.
    4: "Undetermined",  # The crispr repeat is of an undefined repeat family.
    5: "DistantlyConfirmed",  # There is a correct crispr but not near.
    6: "Confirmed",  # The crispr repeat and cas locus correct.
}

def determine_one_caslocus(cas_item, crispr, cas2crispr_links):
    """Returns (tuple, len == 2):
        casclasstype {'Single', 'Uncomplete', 'Ambiguous'}
        casclassarray_possibilities ( see above this function in source file )
    """
    k, v = cas_item
    castypeset = set(v[0].split("|"))
    casclasstype = v[1].split("_")[0].capitalize()
    match = set()
    casclassarray_poss = 1  # Loner
    for crisprlink in cas2crispr_links[k]:
        # this if only reached if there is a crisprlink
        if not castypeset & classified_repeats:
            # If the castypeset has no cas-typing that is a classifiable repeat
            return "|".join(castypeset), casclasstype + "Undefined"
        repeat_family = crisprrepeat2castype(crispr[crisprlink][0])
        if not repeat_family:
            # If the linked crispr has no classified repeat
            casclassarray_poss = max(4, casclassarray_poss)  # Undetermined
            continue
        repeat_combination = castypeset & repeat_family
        if repeat_combination:
            # If there is a match we collect that in match
            match.union(repeat_combination)
            casclassarray_poss = max(6, casclassarray_poss)
            continue
        # Else: There is no match between the castypeset and the repeat family
        casclassarray_poss = max(2, casclassarray_poss)  # Wrong
    if match:
        return "|".join(match), \
                casclasstype + "Confirmed"
    all_repeat_families = {crisprrepeat2castype(c[0]) for c in crispr.values()}
    match = all_repeat_families & castypeset
    if match:
        return "|".join(match), \
            casclasstype + "DistantlyConfirmed"
    return "|".join(castypeset), \
        casclasstype + casclassarray_possibilities[casclassarray_poss]


def say_one_type(list_of_types):
    """say_one_type( collect all types by split on | and rsplit (with re.sub)
    the subtype to find compromis type."""
    types_list = {s for t in list_of_types for s in t.split("|")}
    if not list_of_types:
        return "NA"
    if len(types_list) == 1:
        return types_list.pop()
    no_single_cas = types_list - {"CAS"}
    if not no_single_cas:
        return "CAS"
    typed, subtypes = map(set, zip(*[megasplit(r'(?<=V|I)(?!V|I)', s) for s in
                                     no_single_cas]))
    if len(typed) == 1:
        if len(subtypes) > 3:
            return typed.pop()
        ans = typed.pop() + "/".join(sorted(subtypes))
        return ans
    return "CAS"


def combine_types(list_of_types):
    """combine multipe types as one with a /"""
    return say_one_type(list_of_types)
    # This is probably duplicate code.
    if not list_of_types:
        return "NA"
    types_list = {s for t in list_of_types for s in t.split("|")}
    types_list_iter = map(lambda x: megasplit(r'(?<=V|I)(?!V|I)', x),
                          types_list)
    first = next(types_list_iter)
    return first[0] + "/".join(sorted(set([first[1]]
                               + list(map(itemgetter(1), types_list_iter)))))


def split_type_back(possible_combined_type):
    """Return types back as set and split again for combined types"""
    if possible_combined_type == 'CAS':
        return {'CAS'}
    pre, post = megasplit(r'(?<=V|I)(?!V|I)', possible_combined_type)
    return {pre + p for p in post.split('/')}


def determine_cas(cas, crispr, cas2crispr_links):
    amm_genes = sum([v[-1] for k, v in cas.items()])
    cas_loci = [determine_one_caslocus(c, crispr, cas2crispr_links)
                for c in cas.items()]

    cas_loci_extra, cas_loci_single = split_list_oncondition(
        lambda c: "Single" in c[1], cas_loci)

    if not cas_loci_single:
        non_confirmed, confirmed = split_list_oncondition(
            lambda c: "Confirmed" in c[1], cas_loci_extra)
        if confirmed:
            Extra = "Extra" if non_confirmed else ""
            return say_one_type(map(itemgetter(0), confirmed)), \
                "UnclassifiedConfirmed" + Extra, amm_genes
        return say_one_type(map(itemgetter(0), cas_loci_extra)), \
            "Unclassified", amm_genes
#        cas_loci_other, cas_loci_uncomplete = split_list_oncondition(
#                lambda c: "Uncomplete" in c[1], cas_loci_extra)
#        if cas_loci_uncomplete:
#            if cas_loci_other:
#                return say_one_type(cas_loci_uncomplete + cas_loci_other), \
#                    "MultipleUncomplete"
#            else:
    Extra = "Extra" if cas_loci_extra else ""

    accompanied, lonely = split_list_oncondition(
        lambda c: "Loner" in c[1], cas_loci_single)

    if accompanied:
        Accompanied = "Accompanied"
        if lonely:
            Accompanied = "AdditionalAccompanied"
        cas_loci_single = accompanied
    else:
        Accompanied = ""
        cas_loci_single = lonely

    non_confirmed, confirmed = split_list_oncondition(
        lambda c: "Confirmed" in c[1], cas_loci_single)

    if not confirmed:
        the_type = say_one_type(map(itemgetter(0), non_confirmed))
        if len(non_confirmed) == 1 or (
                "/" not in the_type and the_type != "CAS"):
            return non_confirmed.pop()[0], \
                Accompanied + "SingleUnknown" + Extra, amm_genes
        return the_type, Accompanied + "MultipleUnknown" + Extra, amm_genes

    close, distantly = split_list_oncondition(
        lambda c: "DistantlyConfirmed" in c[1], confirmed)
    if close:
        the_type = combine_types(map(itemgetter(0), confirmed))
        if len(confirmed) == 1 or ("/" not in the_type and the_type != "CAS"):
            consensus_type = combine_types(map(
                itemgetter(0), chain(confirmed, non_confirmed)))
            if len(confirmed) + len(non_confirmed) == 1 or (
                    "/" not in consensus_type and consensus_type != "CAS"):
                return close.pop()[0], \
                    Accompanied + "SingleConfirmed" + Extra, amm_genes
            return close.pop()[0], \
                Accompanied + "MultipleSingleConfirmed" + Extra, amm_genes
        return the_type, Accompanied + "MultipleConfirmed" + Extra, amm_genes
    if len(distantly) == 1:
        if non_confirmed:
            return close.pop()[0], \
                Accompanied + "MultipleSingleDistantlyConfirmed" + Extra, \
                amm_genes
        return close.pop()[0], \
            Accompanied + "SingleDistantlyConfirmed" + Extra, amm_genes

    # implicit `len(comfirmed) > 1`
    return combine_types(map(itemgetter(0), confirmed)), \
        Accompanied + "MultipleDistantlyConfirmed" + Extra, amm_genes


def determine_one_crispr(crispr_item, cas, crispr2cas_links, genome_type):
    k, v = crispr_item
    compatible_type = set(crisprrepeat2castype(v[0]))
    casclassarray_poss = 1  # Loner
    for caslink in crispr2cas_links[k]:
        # this is only reached if there is a crisprlink
        if not compatible_type:
            return (*v, "UndeterminedLinked")
        castypeset = set(cas[caslink][0].split('|'))
        if not castypeset & classified_repeats:
            casclassarray_poss = max(3, casclassarray_poss)  # Undefined
        if castypeset & compatible_type:
            return (*v, "Confirmed")

    if not compatible_type:
        return (*v, "UndeterminedLoner")
    if split_type_back(genome_type[0]) & compatible_type:
        return (*v, "DistantlyConfirmed")
    return (*v, casclassarray_possibilities[casclassarray_poss])


def determine_crispr(crispr, cas, crispr2cas_links, genome_type):
    res = []
    for crispr_item in crispr.items():
        res.append(determine_one_crispr(crispr_item, cas, crispr2cas_links,
                                        genome_type))
    return res


def print_results(res, genome):
    # (key, begin, end) -> (type, class, number of genes)
    cas = {(F[0], F[3], F[4]): (F[1], F[2], F[-1].count('|') + 1) for F in res
           if F[1] != "NULL"}
    # (key, begin, end) -> (type, ID)
    crispr = {(F[0], F[6], F[7]): (F[8], F[5]) for F in res if F[-2] != "NULL"}
    # cas_key -> crispr_key if distance < 20.000
    cas2crispr_links = defaultdict(list)
    for F in res:
        if "NULL" not in (F[2], F[-2]) \
                and min(F[6] - F[4], F[3] - F[7]) < array2crispr_distance:
            cas2crispr_links[F[0], F[3], F[4]].append((F[0], F[6], F[7]))

    crispr2cas_links = defaultdict(list)
    for k, v in cas2crispr_links.items():
        for d in v:
            crispr2cas_links[d].append(k)

    castype, casclass, amm_genes = determine_cas(cas, crispr, cas2crispr_links)
    genome_out(genome, castype, casclass, amm_genes)
    for c in determine_crispr(crispr, cas, crispr2cas_links,
                              (castype, casclass)):
        crispr_out(genome, combine_types(crisprrepeat2castype(c[0])), c[1],
                   c[2])





last_genome = ""
results = []
for F in csv.reader(fileinput.input(args.inputfiles), delimiter=args.sep):
    genome = F[0].split('@')[0]
    if genome != last_genome:
        if last_genome:
            print_results(results, last_genome)
        results = []
        last_genome = genome
    F[3], F[4], F[6], F[7] = map(lambda x: 'NULL' if x == 'NULL' else int(x),
                                 [F[3], F[4], F[6], F[7]])
    results.append(F)

print_results(results, last_genome)

log("[INFO] Program %s finished" % ID)
