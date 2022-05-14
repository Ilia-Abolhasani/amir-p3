import sys
import tqdm
import json
import numpy as np
import pandas as pd


class dotdict(dict):
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


with open('./src/config/titles.json') as json_file:
    titles = dotdict(json.load(json_file))

header = True
for chunk in tqdm(pd.read_csv("./Result/ct_analizer.csv", chunksize=10 ** 5)):
    level1 = chunk[chunk[titles.msg] == '-']
    level1 = level1[level1[titles.star_branching] != 'yes']
    level1 = level1[level1[titles.domain_inter_struc] != 'yes']
    level1 = level1[level1[titles.domain_star_inter_struc] != 'yes']
    level1 = level1[level1[titles.n_term_struc] != "no residues between miR and miR*"]

    level1[titles.num_of_lnk_res] = level1[titles.num_of_lnk_res].apply(lambda x: int(x))
    level1 = level1[level1[titles.num_of_lnk_res] >= 0]

    level1.to_csv("./Result/result_level1_filter.csv", header=header, mode='a', index=False)
    header = False
