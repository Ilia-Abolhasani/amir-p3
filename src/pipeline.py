#!/usr/bin/env python
import os
import sys
import json
import time
import glob
import math
import shutil
import argparse
import numpy as np
import pandas as pd
from read_configs import *
import mirbase
from secondary_structure import secondary_structure
from utils import *
from filter import filter_run
from postprocess import postprocess
from ct_analizer import get_row
import blastN
import extension
import dimond
from networkx.algorithms.clique import find_cliques as maximal_cliques
import networkx
from subprocess import Popen, PIPE, STDOUT
import urllib.parse
import multiprocessing as mp

from tqdm.contrib.concurrent import process_map
from tqdm.notebook import tqdm


def start(
        input_genome_name,
        experiment,
        experiment_dir,
        mirbase_dir,
        nonconformity,
        flanking_value,
        seed_start,
        seed_end,
        hit_threshold,
        precursor_threshold,
        boi_threshold,
        num_cpus,
        secondary_structure_method
):
    input_genome_path = f"{experiment_dir}/{experiment}/{input_genome_name}"
    temp_path = f"{experiment_dir}/{experiment}/Temp"
    result_path = f"{experiment_dir}/{experiment}/Result"
    current_path = os.getcwd()

    if not os.path.exists(temp_path):
        os.mkdir(temp_path)

    if not os.path.exists(result_path):
        os.mkdir(result_path)

    # download mirbase if not exist
    if not os.path.exists(mirbase_dir):
        mirbase.download(mirbase_dir)
    selected = mirbase.select_mirs(mirbase_dir)
    df_to_fasta(selected, f"{temp_path}/mature_microRNA_queries.fasta")

    # Remove redundant cdhit-est
    os.system(
        "./software/cdhit/cd-hit-est -i ./{temp_path_f}/mature_microRNA_queries.fasta  -o ./{temp_path_f}/NR_mature_microRNA_queries.fasta     -c 1 -r 0 -G 1 -g 1 -b 30 -l 10 -aL 0 -AL 99999999 -aS 0     -AS 99999999 -s 0 -S 0"
    )

    # reformat
    with open(f"{temp_path}/NR_mature_microRNA_queries.fasta.clstr", "r") as file:
        text = file.read()
    lines = [line for line in text.split("\n") if len(line) > 0]
    cluster = []
    seqid = []
    last_cluster = ""
    for l in lines:
        if l[0] == ">":
            last_cluster = l.replace(">Cluster ", "C")
        else:
            cluster.append(last_cluster)
            seqid.append(l.split(", >")[1].split("...")[0])
    seq2cluster = pd.DataFrame({"seqid": seqid, "cluster": cluster})

    df = fasta_to_df(f"{temp_path}/mature_microRNA_queries.fasta")
    df["accession"] = df["tag"].apply(lambda x: x.split(" ")[0])
    seq2cluster = pd.merge(
        df, seq2cluster, how="inner", left_on="accession", right_on="seqid"
    )

    seq2cluster = pd.merge(seq2cluster, mature, how="inner", left_on="tag", right_on="tag")[
        ["cluster", "seqid", "tag", "confidence"]
    ]
    seq2cluster.to_csv(f"{temp_path}/seq2cluster.csv", index=False)

    df = fasta_to_df(f"{temp_path}/NR_mature_microRNA_queries.fasta")
    df["tag"] = df["tag"].apply(lambda x: x.split(" ")[0])
    df = pd.merge(df, seq2cluster, how="inner", left_on="tag", right_on="seqid")[
        ["cluster", "data"]
    ]
    lines = []
    df.apply(lambda row: lines.append(
        f">{row['cluster']}\n{row['data']}\n"), axis=1)
    with open(f"{temp_path}/BLASTn_queries.fasta", "w") as file:
        file.write("".join(lines))

    # BlastN
    blastN.start(input_genome_path,
                 temp_path,
                 nonconformity,
                 flanking_value,
                 num_cpus)
    # Extention with bed tools
    extension.start(input_genome_path,
                    temp_path)

    # Dimond
    dimond.start(temp_path)

    # Secondary structure prediction
    secondary_structure(secondary_structure_method, result_path, num_cpus)

    # CTAnalizer
    base = f"{result_path}/secondary_structure/{secondary_structure_method}/"
    df = fasta_to_df(f"{temp_path}/extended_modified_non_coding.txt")
    index_list = []
    for index, row in df.iterrows():
        tag = reformat(row["tag"])
        if len(glob.glob(f"{base + tag}/*.ct")) != 0:
            index_list.append(index)
    df = df.iloc[index_list, :]
    print(df.shape)

    def run(tag, path, extra):
        try:
            return get_row(tag, path, extra)
        except Exception as e:
            print(str(e), tag)
            return pd.Series()

    def get_df_by_tag(tag, extra=0):
        ct_files = glob.glob(f"{base}{reformat(tag)}/*.ct")
        return pd.Series(ct_files).apply(lambda path: run(tag, path, extra))

    # Apply on current data
    seq2cluster = pd.read_csv(f"{temp_path}/seq2cluster.csv")
    seq2cluster["tag"] = seq2cluster.groupby(["cluster"])["tag"].transform(
        lambda x: ",".join(x)
    )
    seq2cluster["seqid"] = seq2cluster.groupby(["cluster"])["seqid"].transform(
        lambda x: ",".join(x)
    )
    seq2cluster = seq2cluster.drop_duplicates()
    tag2cluster = pd.read_csv(
        f"./{temp_path}/pipe_seprated_location_list.csv", sep="\t")
    tag2cluster["location_tag"] = tag2cluster["location_tag"].apply(
        lambda x: x[1:])
    data = pd.merge(
        seq2cluster, tag2cluster, how="inner", left_on="cluster", right_on="qseqid"
    )
    data["Reference miRNA cluster"] = data["cluster"]
    data["Reference miRNA IDs"] = data["seqid"]
    data["Reference miRNA IDs and species"] = data["tag"]
    data = data[
        [
            "location_tag",
            "Reference miRNA cluster",
            "Reference miRNA IDs",
            "Reference miRNA IDs and species",
            "confidence",
        ]
    ]

    rcols_ref = [
        "Reference miRNA cluster",
        "Reference miRNA IDs",
        "Reference miRNA IDs and species",
    ]

    rcols_boi = ["boi seq", "boi name", "boi dotbracket"]

    rcols = [*rcols_ref, *rcols_boi]

    rcols_dg = [*rcols, "delta G"]

    def selection(row):
        global repeted
        tuple_row = tuple(row)
        if tuple_row not in repeted:
            repeted[tuple_row] = row.name
            return True
        return False

    def boi_selection(row):
        global repeted_boi
        tuple_row_boi = tuple(row[rcols_boi])
        dg = row["delta G"]
        if tuple_row_boi not in repeted_boi:
            repeted_boi[tuple_row_boi] = {
                "counter": 1,
                "dg": dg,
                "ref clusters": row["Reference miRNA cluster"],
                "ref ids": row["Reference miRNA IDs"],
                "ref species": row["Reference miRNA IDs and species"],
                "lock": False,
            }
        else:
            value = repeted_boi[tuple_row_boi]
            value["counter"] += 1
            value["dg"] = min(value["dg"], dg)
            value["ref clusters"] += "," + row["Reference miRNA cluster"]
            value["ref ids"] += "," + row["Reference miRNA IDs"]
            value["ref species"] += "," + \
                row["Reference miRNA IDs and species"]
            repeted_boi[tuple_row_boi] = value

    os.system("rm {result_path}/ct_analizer.csv")
    chunksize = 1 * (10**4)
    max_workers = mp.cpu_count() - 4
    num_terminal = 5  # acceptable_terminal_structures

    repeted = {}
    repeted_boi = {}
    header = True
    orders = None
    arr = np.array_split(df["tag"], max(df["tag"].shape[0] // chunksize, 1))
    for chunk in tqdm(arr):
        dfs = []
        for row in process_map(
            get_df_by_tag, chunk, tqdm_class=tqdm, max_workers=max_workers, chunksize=5
        ):
            dfs.append(row)
        chunk = pd.concat(dfs, axis=0)
        chunk = pd.merge(
            data, chunk, how="right", left_on="location_tag", right_on="seq name"
        )
        del chunk["location_tag"]
        if header:
            orders = chunk.columns
        for col in orders:
            if col not in chunk.columns:
                chunk[col] = np.nan

        for col in chunk.columns:
            if col not in orders:
                print(f"Error in {col}")
        chunk = chunk.reindex(columns=orders)
        chunk = chunk.replace(np.nan, "-").replace("", "-")
        # delete repeated
        selected = chunk[rcols].apply(lambda row: selection(row), axis=1)
        chunk = chunk[selected]

        # cluster refs
        chunk[rcols_dg].apply(lambda row: boi_selection(row), axis=1)
        chunk.to_csv(f"{result_path}/ct_analizer.csv",
                     header=header, mode="a", index=False)
        header = False

    os.system("rm {result_path}/ct_analizer_clustered.csv")

    def isKeepCluster(row):
        global repeted_boi
        dg = row["delta G"]
        if row["boi name"] == "-":
            return False
        tuple_row_boi = tuple(row[rcols_boi])
        value = repeted_boi[tuple_row_boi]
        if value["counter"] == 1:
            return True
        if value["dg"] != dg:
            return False
        if value["lock"]:
            return False
        value["lock"] = True
        repeted_boi[tuple_row_boi] = value
        return True

    def makeCluster(row):
        tuple_row_boi = tuple(row[rcols_boi])
        value = repeted_boi[tuple_row_boi]
        if value["counter"] != 1:
            for ref_c in ["ref clusters", "ref ids", "ref species"]:
                value[ref_c] = value[ref_c].replace(
                    " ,", ",").replace(", ", ",")
                value[ref_c] = value[ref_c].split(",")
                value[ref_c] = set(value[ref_c])
                value[ref_c] = ", ".join(value[ref_c])
            row["Reference miRNA cluster"] = value["ref clusters"]
            row["Reference miRNA IDs"] = value["ref ids"]
            row["Reference miRNA IDs and species"] = value["ref species"]
        return row

    header = True
    for chunk in tqdm(pd.read_csv(f"{result_path}/ct_analizer.csv", chunksize=10**5)):
        chunk = chunk[chunk[rcols_dg].apply(
            lambda row: isKeepCluster(row), axis=1)]
        chunk = chunk.apply(lambda row: makeCluster(row), axis=1)
        chunk.to_csv(f"{result_path}/ct_analizer_clustered.csv",
                     mode="a", index=False)
        header = False

    # Filters
    os.system("rm {result_path}/result_level1_filter.csv")
    filter1_run(
        input_file=f"{result_path}/ct_analizer_clustered.csv",
        output_file=f"{result_path}/result_level1_filter.csv",
    )

    config = {
        "delta_g_min": -999,
        "delta_g_max": 1,
        "hit_len_min": 21,
        "hit_len_max": 21,
        "hit_complementarity_percentage_min": 0.5,
        "hit_complementarity_percentage_max": 1.0,
        "number_of_terminal_structure_min": 0,
        "number_of_terminal_structure_max": 5,
        "boi_gc_content_min": 45,
        "boi_gc_content_max": 94,
        "num_of_linking_residues_min": 5,
        "num_of_linking_residues_max": 159,
        "hit_gc_content_percentage_min": 37,
        "hit_gc_content_percentage_max": 86,
        "precursor_mfei_min": 0.87,
        "precursor_mfei_max": 1.3695556794836854,
        "border_line_mismatch_max": 0,
        "border_line_bulge_max": 0,
        "border_line_internal_max": 0,
        "total_num_of_nonmatching_positions": 5,
        "total_num_of_mismached_positions": 5,
        "total_num_of_positions_in_bulges_and_loops": 2,
        "max_allowed_mismatch_size_in_hit_region": 2,
        "max_allowed_bulge_size_in_hit_region": 1,
        "max_allowed_internal_loop_size_in_hit_region": 3,
        "max_allowed_hsbl_ssbl_size": 2,
        "minimum_required_clear_region": 0,
        "acceptable_num_for_hit_locations_in_bulges_or_loops": 2,
        "acceptable_num_for_unmatched_locations_in_hit_region": 5,
        "delete_if_mature_duplex_involvement_in_apical_loop": "YES",
        "border_line_structure_allowance": "NOT ACCEPTED",
    }

    filter2(
        input_file=f"{result_path}/result_level1_filter.csv",
        output_file=f"{result_path}/result_level2_filter.csv",
        config=config,
    )

    # Cluster JSC
    result = pd.read_csv(f"{result_path}/result_level2_filter.csv")
    print(result.shape)
    result.head(2)

    def jaccard(A, B):
        a1 = int(A.split("-")[0])
        a2 = int(A.split("-")[1])
        b1 = int(B.split("-")[0])
        b2 = int(B.split("-")[1])
        s1 = set([i for i in range(min(a1, a2), max(a1, a2) + 1)])
        s2 = set([i for i in range(min(b1, b2), max(b1, b2) + 1)])
        intersection = len(s1.intersection(s2))
        union = (len(s1) + len(s2)) - intersection
        return float(intersection) / union

    def same2cluster(same_dict, threshold=0.8):
        counter = 1
        for key in same_dict:
            item2cluster = {}
            SET = list(set(same_dict[key]))
            G = networkx.Graph()
            # add nodes
            for s in SET:
                G.add_node(s)
            # add edges
            for i in range(0, len(SET)):
                for j in range(i + 1, len(SET)):
                    if jaccard(SET[i], SET[j]) >= threshold:
                        G.add_edge(SET[i], SET[j],
                                   weight=jaccard(SET[i], SET[j]))

            # get maximal
            for clique in maximal_cliques(G):
                unique = str(counter).zfill(4)
                counter += 1
                for item in clique:
                    if item in item2cluster:
                        item2cluster[item].append(unique)
                    else:
                        item2cluster[item] = [unique]
            same_dict[key] = item2cluster
        return same_dict

    # hit jaccard similarity
    same_strand_hit = {}

    def same_strand(row):
        global same_strand_hit
        chrom = row["chromosome"]
        sign = row["sign"]
        hit = row["hit position on chromosome"]
        key = f"{chrom}{sign}"
        if key in same_strand_hit:
            same_strand_hit[key].append(hit)
        else:
            same_strand_hit[key] = [hit]

    rcols_hit = ["chromosome", "sign", "hit position on chromosome"]
    result[rcols_hit].apply(lambda row: same_strand(row), axis=1)
    same_strand_hit = same2cluster(same_strand_hit, threshold=hit_threshold)

    def _f(row):
        key = f"{row['chromosome']}{row['sign']}"
        return same_strand_hit[f"{key}"][row["hit position on chromosome"]]

    result["hit cluster number"] = result[rcols_hit].apply(
        lambda row: _f(row), axis=1)

    rcols_boi = ["boi seq", "boi name", "boi dotbracket"]
    same_strand_boi = {}

    def same_strand(row):
        global same_strand_boi
        seq = row["boi seq"].lower()
        name = row["boi name"]
        dotbracket = row["boi dotbracket"].lower()
        if pd.isna(name) or name == "-":
            return
        name = name.split("|")
        loction = name[0] + name[1] + name[2]
        hit = name[3]
        key = f"{seq}{dotbracket}{loction}"
        if key in same_strand_boi:
            same_strand_boi[key].append(hit)
        else:
            same_strand_boi[key] = [hit]

    result[rcols_boi].apply(lambda row: same_strand(row), axis=1)

    same_strand_boi = same2cluster(same_strand_boi, threshold=boi_threshold)

    def _f(row):
        seq = row["boi seq"].lower()
        name = row["boi name"]
        dotbracket = row["boi dotbracket"].lower()
        if pd.isna(name) or name == "-":
            return
        name = name.split("|")
        loction = name[0] + name[1] + name[2]
        hit = name[3]
        key = f"{seq}{dotbracket}{loction}"
        return same_strand_boi[f"{key}"][hit]

    result["boi cluster number"] = result[rcols_boi].apply(
        lambda row: _f(row), axis=1)

    same_strand_precursor = {}
    rcols_pre = ["precursor seq", "precursor name", "precursor dotbracket"]

    def same_strand(row):
        global same_strand_precursor
        seq = row["precursor seq"].lower()
        name = row["precursor name"]
        dotbracket = row["precursor dotbracket"].lower()
        if pd.isna(name) or name == "-":
            return
        name = name.split("|")
        loction = name[0] + name[1] + name[2]
        hit = name[3]
        key = f"{seq}{dotbracket}{loction}"
        if key in same_strand_precursor:
            same_strand_precursor[key].append(hit)
        else:
            same_strand_precursor[key] = [hit]

    result[rcols_pre].apply(lambda row: same_strand(row), axis=1)
    same_strand_precursor = same2cluster(
        same_strand_precursor, threshold=precursor_threshold
    )

    def _f(row):
        seq = row["precursor seq"].lower()
        name = row["precursor name"]
        dotbracket = row["precursor dotbracket"].lower()
        if pd.isna(name) or name == "-":
            return
        name = name.split("|")
        loction = name[0] + name[1] + name[2]
        hit = name[3]
        key = f"{seq}{dotbracket}{loction}"
        return same_strand_precursor[f"{key}"][hit]

    result["precursor cluster number"] = result[rcols_pre].apply(
        lambda row: _f(row), axis=1
    )

    hit2cluster = {}
    hit_unique = result["hit seq"].unique()
    for i in range(0, hit_unique.shape[0]):
        hit2cluster[hit_unique[i]] = str(i + 1).zfill(4)
    result["identical hit cluster"] = result["hit seq"].apply(
        lambda hit: hit2cluster[hit])

    result["seed region"] = result["hit seq"].apply(
        lambda hit: hit[seed_start - 1: seed_end]
    )

    result.to_csv(
        f"{result_path}/result_level2_filter_clustered.csv", index=False)
    os.system(
        "zip -r {result_path}/result_level2_filter_clustered.zip {result_path}/result_level2_filter_clustered.csv"
    )

    # BlastX
    os.system("makeblastdb -in ./NR/nr -dbtype prot -out ./NR/nr_database")
    os.ipython().system(
        'blastx -query ./input_blastx.txt         -db ./NR/nr_database         -out ./Temp/BlastX/blastx         -num_threads 20         -evalue 1e-3         -outfmt "6 qseqid sseqid qstart qend evalue bitscore score length frames qframe qcovs qcovhsp staxids"'
    )
    blx = pd.read_csv("./Temp/BlastX/blastx", sep="\t", header=None)
    blx.columns = "qseqid sseqid qstart qend evalue bitscore score length frames qframe qcovs qcovhsp staxids".split(
        " "
    )
    coding_seq = blx["qseqid"].unique()
