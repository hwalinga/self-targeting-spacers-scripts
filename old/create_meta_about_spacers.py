import sys
import re

DOTS_IN_EXTENSION = 0
GENOME_PATH = sys.argv[1]
# GENOME_PATH = "output/ESC_AA8779AA_AS.crispr"

base_file = GENOME_PATH.split('/')[-1]
genome = ".".join(base_file.split(".")[:-(DOTS_IN_EXTENSION + 1)])

family_reg_prog = re.compile(':\s(.*)\s\[')
with open(GENOME_PATH) as f:
    for line in f: 
        if line.startswith('Array'):
            spacer_array = line.split("*", maxsplit=1)[0].strip()
            spacer_array = "_".join(spacer_array.split()) #EG: Array_1_583-401
            line_splitted = next(f).split()
            direction = line_splitted[-1][0]
            spacer_array += '_' + direction
            print(genome, end='\t')
            print(spacer_array, end='\t')
            contig = line_splitted[0][1:]
            print(contig, end='\t')
            while True:
                line = next(f)
                if line.startswith('# Array'):
                    if line.strip() == '# Array family : NA':
                        fam = 'NA'
                    else:
                        fam = family_reg_prog.search(line).group(1)
                    print(fam, end='\n')
                    break



