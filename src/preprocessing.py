import os
import json
import numpy as np
import pandas as pd
import multiprocessing as mp
from numpy.random import randint, rand
from sklearn.preprocessing import OneHotEncoder
from keras.utils import np_utils


def reverse_complement(dna):
    out = ""
    for d in dna[::-1]:
        if d == "A":
            out += "T"
        if d == "T":
            out += "A"
        if d == "C":
            out += "G"
        if d == "G":
            out += "C"
    return out


k = 4
nuc = ["A", "C", "G", "T"]
np.random.seed(0)
tnf = set(
    ["".join(np.random.choice(nuc, size=k, replace=True))
     for i in range(0, 2000)]
)
tnf = sorted(list(tnf))
_tnf = {}
for kmer in tnf:
    if reverse_complement(kmer) not in _tnf:
        _tnf[kmer] = 0


def tnf_calc(dna):
    dna = dna.replace("N", "")
    counter = 0
    out = _tnf.copy()
    for i in range(0, len(dna) - (k - 1)):        
        kmer = dna[i: i + k]        
        if kmer in out:
            out[kmer] += 1
        else:
            out[reverse_complement(kmer)] += 1
        counter += 1
    return pd.Series(out) / counter


def nuc_one_hot(n):
    n = n.upper()
    out = {"A": 0, "C": 0, "G": 0, "T": 0}
    if n in out:
        out[n.upper()] = 1
    return pd.Series(out)


def com_one_hot(n):
    n = n.upper()
    if n == "-":
        n = "-/-"
    out = np.zeros(8)
    if n[0] == "A":
        out[0] = 1
    if n[0] == "C":
        out[1] = 1
    if n[0] == "G":
        out[2] = 1
    if n[0] == "T":
        out[3] = 1
    if n[2] == "A":
        out[4] = 1
    if n[2] == "C":
        out[5] = 1
    if n[2] == "G":
        out[6] = 1
    if n[2] == "T":
        out[7] = 1
    return pd.Series(out)


def get_max_in_region(row, type_str, size_str, region):
    out = 0
    mismatch_type = row[type_str]
    for i in range(len(mismatch_type)):
        if mismatch_type[i] == region:
            if row[size_str][i] > out:
                out = row[size_str][i]
    return out


def get_number_in_region(row, type_str, size_str, region):
    counter = 0
    mismatch_type = row[type_str]
    for i in range(len(mismatch_type)):
        if mismatch_type[i] == region:
            counter += 1
    return counter


def get_sum_in_region(row, type_str, size_str, region):
    out = 0
    mismatch_type = row[type_str]
    for i in range(len(mismatch_type)):
        if mismatch_type[i] == region:
            out += row[size_str][i]
    return out


def preprocessing(df, mu=None, std=None, make_standard=True):    
    result = df.copy()
    result = result.reset_index(drop=True)

    cols = [
        "hit len",
        "flanking GC content",
        "flanking MFEI",
        "hit GC content",
        "hit complementarity percentage",
        "num of linking residues",
        "boi GC content",
        "boi delta G",
        "boi AMFE",
        "boi MFEI",
        "branch#1 apical loop size",
        "precursor gc content",
        "precursor delta G",
        "precursor AMFE",
        "precursor MFEI",
        "number of terminal structures",
        "primary stem length",
        "mismatch",
        "bulge",
        "internal loop",
        "base structure corrected length",
        "primary stem corrected length",
        "Loop distal junction distance",
        "Loop proximal junction distance",
        "sum of residue in terminal loop",
        "acceptable num for hit locations in bulges or loops",
        "acceptable num for unmatched locations in hit region * 2",
        "acceptable num for unmatched locations in hit region",
        "acceptable num for hit locations in bulges or loops mayers",
        "acceptable num for unmatched locations in hit region mayers * 2",
        "acceptable num for unmatched locations in hit region mayers",
        "total num of mismached positions",
        "total num of nonmatching positions",
        "total num of positions in bulges and loops",
        "mature duplex involvement in apical loop",
    ]
    for region in [
        "loop distal",
        "hit region",
        "distal border line",
        "loop proximal",
        "proximal border line",
    ]:
        _col = ["mismatch type", "mismatch size"]
        result[f"max mismatch in {region}"] = result[_col].apply(
            lambda row: get_max_in_region(
                row, "mismatch type", "mismatch size", region
            ),
            axis=1,
        )
        cols.append(f"max mismatch in {region}")
        result[f"number mismatch in {region}"] = result[_col].apply(
            lambda row: get_number_in_region(
                row, "mismatch type", "mismatch size", region
            ),
            axis=1,
        )
        cols.append(f"number mismatch in {region}")
        result[f"sum mismatch in {region}"] = result[_col].apply(
            lambda row: get_sum_in_region(
                row, "mismatch type", "mismatch size", region
            ),
            axis=1,
        )
        cols.append(f"sum mismatch in {region}")
        _col = ["bulge type", "bulge size"]
        result[f"max bulge in {region}"] = result[_col].apply(
            lambda row: get_max_in_region(
                row, "bulge type", "bulge size", region),
            axis=1,
        )
        cols.append(f"max bulge in {region}")
        result[f"number bulge in {region}"] = result[_col].apply(
            lambda row: get_number_in_region(
                row, "bulge type", "bulge size", region),
            axis=1,
        )
        cols.append(f"number bulge in {region}")
        result[f"sum bulge in {region}"] = result[_col].apply(
            lambda row: get_sum_in_region(
                row, "bulge type", "bulge size", region),
            axis=1,
        )
        cols.append(f"sum bulge in {region}")
        _col = ["internal type", "internal loop total size"]
        result[f"max loop in {region}"] = result[_col].apply(
            lambda row: get_max_in_region(
                row, "internal type", "internal loop total size", region
            ),
            axis=1,
        )
        cols.append(f"max loop in {region}")
        result[f"number loop in {region}"] = result[_col].apply(
            lambda row: get_number_in_region(
                row, "internal type", "internal loop total size", region
            ),
            axis=1,
        )
        cols.append(f"number loop in {region}")
        result[f"sum loop in {region}"] = result[_col].apply(
            lambda row: get_sum_in_region(
                row, "internal type", "internal loop total size", region
            ),
            axis=1,
        )
        cols.append(f"sum loop in {region}")
    mirs = result[["mir type"]]
    X = result[cols]
    # replace inf with max
    for c in X.columns:
        m = X[X[c] != np.inf][c].max()
        X[c].replace([np.inf], m, inplace=True)
    # standardization
    if mu is None or std is None:
        mu = X.mean()
        std = X.std()
    if(make_standard):
        X = (X - mu) / std
    # mir type
    encoder = OneHotEncoder(handle_unknown="ignore")
    encoder_df = pd.DataFrame(encoder.fit_transform(mirs).toarray())
    X = X.join(encoder_df)
    # TNF
    _tnf = pd.DataFrame(df["hit seq"].apply(lambda x: tnf_calc(x)))
    X = X.join(_tnf).reset_index(drop=True)
    # start end nucleotide
    nuc_list = []
    for c in ["-3", "-2", "-1", "", "+1", "+2"]:
        col_name = f"hit start{c} composition"
        _nuc = pd.DataFrame(
            result[col_name].apply(lambda x: com_one_hot(x))
        ).reset_index(drop=True)
        _nuc.columns = [f"start_{c}_{i}" for i in _nuc.columns]
        nuc_list.append(_nuc)
        col_name = f"hit end{c} composition"
        _nuc = pd.DataFrame(
            result[col_name].apply(lambda x: com_one_hot(x))
        ).reset_index(drop=True)
        _nuc.columns = [f"end_{c}_{i}" for i in _nuc.columns]
        nuc_list.append(_nuc)

    X = pd.concat([X, *nuc_list], axis=1)
    # connectivity
    cols = []
    for c in ["-3", "-2", "-1", "", "+1", "+2"]:
        cols.append(f"connectivity hit start{c}")
        cols.append(f"connectivity hit end{c}")
    #for c in range(2,13):
        #cols.append(f"seed connectivity{c}")
    X = X.join(result[cols])
    X = X.astype("float32")
    return [X, mu, std]


def get_target(df, reference):
    Y = df["hit seq"].isin(reference["data"])
    Y = Y.apply(lambda x: 1 if x else 0)
    Y = np_utils.to_categorical(Y, 2)
    Y = Y.astype("float32")
    return Y
