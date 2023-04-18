from tqdm import tqdm
import pandas as pd
from read_configs import read_titles, read_erros


def filter_run(
    input_file="./Result/ct_analizer_clustered.csv",
    output_file="./Result/result_level1_filter.csv",
    chunksize=10**5,
):
    titles = read_titles()
    header = True
    for chunk in tqdm(pd.read_csv(input_file, chunksize=chunksize)):
        level1 = chunk[chunk[titles.msg] == "-"]
        level1 = level1[level1[titles.star_branching] != "yes"]
        level1 = level1[level1[titles.domain_inter_struc] != "yes"]
        level1 = level1[level1[titles.domain_star_inter_struc] != "yes"]
        level1 = level1[
            level1[titles.n_term_struc] != "no residues between miR and miR*"
        ]

        level1[titles.num_of_lnk_res] = level1[titles.num_of_lnk_res].apply(
            lambda x: int(x)
        )
        level1 = level1[level1[titles.num_of_lnk_res] >= 0]

        level1.to_csv(output_file, header=header, mode="a", index=False)
        header = False
