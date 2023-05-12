#!/usr/bin/env python
import os
import pandas as pd
from read_configs import *
from utils import *
from filter import filter_run
from postprocess import postprocess
from ct_analizer import get_row
from subprocess import Popen, PIPE, STDOUT
import multiprocessing as mp

from tqdm.contrib.concurrent import process_map
from tqdm.notebook import tqdm


def start(
        input_genome_path,
        temp_path
):
    # Extention
    os.system(
        "bedtools getfasta -fi {input_genome_path} -fo {temp_path}/extended_original.txt -s -bed {temp_path}/extension_index.bed"
    )
    os.system("rm input_genome.fna.fai")

    # Convert hit region to upper case and other region to lower case
    ext = fasta_to_df(f"{temp_path}/extended_original.txt")
    info = pd.read_csv(f"{temp_path}/hit_index_info.csv")
    info["tag"] = info["tag"].apply(lambda x: x[1:])

    ext = ext.sort_values(by=["tag"]).reset_index()
    ext["help_tag"] = ext.apply(lambda r: r["tag"] + str(r.name), axis=1)
    del ext["tag"]

    info = info.sort_values(by=["tag"]).reset_index()
    info["help_tag"] = info.apply(
        lambda row: row["tag"] + str(row.name), axis=1)

    def redefined_tag(row):
        tag = row["tag"]
        [sstart, send] = tag.split(":")[-1].split("(")[0].split("-")
        sstart = int(sstart) + 1
        sign = tag.split("(")[-1].split(")")[0]
        return f"{tag.split(':')[0]}|{sign}|{sstart}-{send}|{row['hit_start']+1}-{row['hit_end']}"

    info["tag"] = info.apply(lambda row: redefined_tag(row), axis=1)
    ext = pd.merge(ext, info, how="inner", on="help_tag")

    def emphasis_hit(row):
        seq = list(row["data"].lower())
        s = row["hit_start"]
        e = row["hit_end"]
        seq[s:e] = list("".join(seq[s:e]).upper())
        return "".join(seq)

    ext["data"] = ext.apply(lambda row: emphasis_hit(row), axis=1)
    df_to_fasta(ext[["tag", "data"]], f"{temp_path}/extended_modified.txt")
