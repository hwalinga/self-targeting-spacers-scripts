from collections import defaultdict
dist_to_array = 20000 #(CRISPRdisco)

OUT_FILE = 're_systems_comb'

crispr_arrays = defaultdict(lambda: defaultdict(set))
for line in open("all_meta_info"):
    genome, array, contig, system = line.strip().split()
    crispr_arrays[genome][contig].add('_'.join([array, system]))

virsorter = defaultdict(lambda: defaultdict(str))
for line in open("all.phages.coords.tsv"):
    genome, contig, phage = line.strip().split(maxsplit=2)
    virsorter[genome][contig] = phage

out = open(OUT_FILE, 'w')
with open('reclassified_crispr_systems') as f:
    header = next(f) 
    new_header = '%s\t%s\t%s\n' % (header[:-1], "on_phage", "spacer_arrays")
    out.write(new_header)
    for line in f:
        no_linebreak_line = line[:-1]
        line = line.strip()
        splitted_line = line.split('\t')
        genome = splitted_line[0]
        contig = splitted_line[5]
        s = int(splitted_line[9])
        e = int(splitted_line[6])
        arrays = set()
        for array in crispr_arrays[genome][contig]:
            _, _, coords, _, system =  array.split('_')
            ss, se = map(int, coords.split('-', maxsplit=1))
            if abs(min(ss, se) - max(s, e)) < dist_to_array or \
                abs(max(ss, se) - min(s, e)) < dist_to_array:
                arrays.add(array) 

        phage = virsorter[genome].get(contig)
        if not phage:
            phage_id = "0"
        elif len(phage) == 1:
            phage_id = phage
        else:
            phage_id, coords, _ = phage.split()
            vs, ve = map(int, coords.split('-', maxsplit=1))
            if not(vs <= s <= ve or vs <= e <= ve):
                phage_id = "0"
        new_line = '%s\t%s\t%s\n' % (no_linebreak_line, phage_id, ",".join(arrays))
        out.write(new_line)


