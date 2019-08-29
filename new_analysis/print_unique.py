import sys
from collections import defaultdict
from functools import partial
from itertools import groupby
from operator import itemgetter

sep = "\t"
out = partial(print, sep=sep, file=sys.stdout)

for key, F_it in groupby(
    map(lambda l: l.strip().split(sep), sys.stdin), key=itemgetter(0)
):
    A = defaultdict(list)
    for F in F_it:
        A[F[0], F[1], F[2], F[3], F[4]].append(F)
    for l in A.values():
        if len(l) == 1:
            out(*l[0])
            continue
        for i in l:
            if i[-1] == "2":
                out(*i)
                break
        else:
            out(*l[0])
