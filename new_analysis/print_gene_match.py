import sys
from collections import defaultdict
from functools import partial
from itertools import groupby
from operator import itemgetter

sep = "\t"
out = partial(print, sep=sep, file=sys.stdout)

for key, F_it in groupby(
    map(lambda l: l.strip().split(sep), sys.stdin), key=lambda x: (x[0], *x[3:6])
):
    protospacer_hit_ori = "-1" if int(key[1]) < int(key[2]) else "1"
    for F in F_it:
        if int(F[-4]) < int(key[1]) < int(F[-3]) and int(F[-4]) < int(key[2]) < int(
            F[-3]
        ):
            if F[-2] == protospacer_hit_ori:
                hit_class = "RNA+"
            else:
                hit_class = "RNA-"
            out(*F[:-6], hit_class, F[-5], F[-1])
            break
    else:
        out(*F[:-6], "INTERGENIC", "NULL", "NULL")
