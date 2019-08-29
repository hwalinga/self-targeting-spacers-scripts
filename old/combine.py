from glob import glob
from collections import defaultdict

spacer_count_dict = defaultdict(dict)
for spacer_count_file in glob("output/*.spacer_count"):
    with open(spacer_count_file) as f:
        for line in f:
            genome, spacer_array, spacer_count = line.strip().split()
            genome = genome.split(".")[0]
            spacer_array = spacer_array[1:]  # remove leading '>'
            spacer_count = int(spacer_count.split("_")[-1])
            spacer_count += 1  # 0 indexed
            spacer_count_dict[genome][spacer_array] = spacer_count

crispr_repeat_family_dict = defaultdict(dict)
with open("meta.info") as f:
    for line in f:
        if line.startswith("output/"):  # new genome
            genome = line.split("/")[-1].split(".")[0]
            continue
        spacer_array, crispr_repeat_family = line.split(":\t")
        spacer_array = "_".join(spacer_array.split("_")[:-1])
        crispr_repeat_family_dict[genome][spacer_array] = crispr_repeat_family

IN_file = "all_franklin_spacers"
OUT_file = "all_franklin_spacers.combined"

OUT_handle = open(OUT_file, "w")
with open(IN_file) as f:
    for line in f:
        genome, spacer_id, rest = line.strip().split(maxsplit=2)
        genome = genome.split("/")[-1].split(".")[0]
        spacer_array = "_".join(spacer_id.split("_")[:-3])
        spacer_count = spacer_count_dict[genome][spacer_array]
        crispr_repeat_family = crispr_repeat_family_dict[genome][spacer_array]
        OUT_handle.write(
            "\t".join(
                [
                    genome,
                    spacer_array,
                    spacer_id,
                    rest,
                    str(spacer_count),
                    crispr_repeat_family,
                ]
            )
        )
OUT_handle.close()
