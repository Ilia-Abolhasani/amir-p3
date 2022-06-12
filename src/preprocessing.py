import os
import json
import numpy as np
import pandas as pd
import multiprocessing as mp
from numpy.random import randint, rand
from sklearn.preprocessing import OneHotEncoder
from keras.utils import np_utils


def reverse_complement(dna):
    dna = dna[::-1]
    dna = dna.replace("A", "1").replace("C", "2").replace("G", "3").replace("T", "4")
    dna = dna.replace("1", "T").replace("2", "G").replace("3", "C").replace("4", "T")
    return dna

k = 4
nuc = ["A", "C", "G", "T"]
tnf = set([''.join(np.random.choice(nuc, size=k, replace=True)) for i in range(0, 20000)])
_tnf = {}
for kmer in tnf:
    if(reverse_complement(kmer) not in _tnf):
        _tnf[kmer] = 0


def tnf_calc(dna):
    counter = 0
    out = _tnf.copy()        
    for i in range(0, len(dna)-(k-1)):        
        kmer = dna[i:i+k]        
        if(kmer in out):
            out[kmer] += 1    
        else:
            out[reverse_complement(kmer)] += 1    
        counter += 1        
    return pd.Series(out)/counter    


def nuc_one_hot(n):    
    n = n.upper() 
    out = {"A":0, "C":0, "G":0, "T":0}
    if(n in out):
        out[n.upper()] = 1
    return pd.Series(out)


# groups = [['A/A', 'G/A', 'A/G', 'G/G'],
#           ['C/C', 'T/T',  'T/C', 'C/T'],
#           ['A/C', 'C/A'],
#           ['T/A', 'A/T'],
#           ['G/C', 'C/G'],
#           ['G/T', 'T/G'],
#           ['-']]

          
# def com_one_hot(n):            
#     out = np.zeros(len(groups))
#     n = n.upper()        
#     for i in range(len(groups)):
#         if n in groups[i]:
#             out[i] = 1
#     return pd.Series(out)

def com_one_hot(n):      
    if(n == "-"):
        n = "-/-"
    n = n.upper()        
    out1 = nuc_one_hot(n[0]) 
    out1.index = [ "1" + col for col in out1.index]
    out2 = nuc_one_hot(n[2])         
    out2.index = [ "2" + col for col in out2.index]
    return pd.concat([out1, out2])


def sum_of_size_in_hit(row, type_str, size_str):
    _sum = 0
    mismatch_type = row[type_str]
    for i in range(len(mismatch_type)):
        if mismatch_type[i] == "hit region":
            _sum += row[size_str][i]
    return _sum


def sum_of_size_in_hit_only_zero(row):
    _sum = 0
    bulge_type = row["bulge type"]
    bulge_strand = row["bulge strand"]
    for i in range(len(bulge_type)):
        if bulge_type[i] == "hit region" and bulge_strand[i] == "zero":
            _sum += row["bulge size"][i]
    return _sum


def number_of_residue(row):
    hit_end = row['hit end']    
    hit_start = row['hit start']
    psep = row['psep']
    if(psep == '-'):
        return 0
    mir_type = row['mir type']
    if mir_type == '5p':
        if psep < hit_end:
            return hit_end - psep
    if mir_type == '3p':
        if psep > hit_start:
            return psep - hit_start
    return 0


def sum_of_size_in_border_line(row, border_type, type_str, size_str, start, end):
    _sum = 0
    _size = row[size_str]    
    _start = row[start]    
    _end = row[end]
    mir_type = row['mir type']     
    mismatch_type = row[type_str]    
    for i in range(len(mismatch_type)):
        if mismatch_type[i] == border_type:
            if border_type == "distal border line":
                if mir_type  == '5p':
                    _sum += _size[i] - _start[i]
                if mir_type == '3p':
                    _sum += _size[i] - _end[i]
            if border_type == "proximal border line":
                if mir_type == '5p':
                    _sum += _size[i] - _end[i]
                if mir_type == '3p':
                    _sum += _size[i] - _start[i]
    return _sum


def check_involvement(row):
    if(row['number of terminal structures'] == "-"):
        return None
    if(row['number of terminal structures'] > 1):
        return True
    start = row['branch#1 apical loop start']
    end = row['branch#1 apical loop end']
    for col in ['hit start', 'hit end', 'star start', 'star end']:
        if(start < row[col] < end ):
            return False                    
    return True


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


def preprocessing(df):        
    result = df.copy()
    cols = ['mismatch type', 'mismatch size']
    sum_missmatch = result[cols].apply(lambda row: sum_of_size_in_hit(row, 'mismatch type', 'mismatch size'), axis=1)

    cols = ['bulge type', 'bulge size']
    sum_bulge = result[cols].apply(lambda row: sum_of_size_in_hit(row, 'bulge type', 'bulge size'), axis=1)
    
    cols = ['internal type', 'internal loop total size']
    sum_internal = result[cols].apply(lambda row: sum_of_size_in_hit(row, 'internal type', 'internal loop total size'), axis=1)

    cols = ['internal type', 'internal loop HSBL']
    sum_internal_hsbl = result[cols].apply(lambda row: sum_of_size_in_hit(row, 'internal type', 'internal loop HSBL'),axis=1)

    cols = ["mir type", 'mismatch type', 'mismatch size', 'mismatch start', 'mismatch end'] 
    sum_missmatch_border_proximal = result[cols].apply(lambda row: sum_of_size_in_border_line(row, 'proximal border line', 'mismatch type', 'mismatch size', 'mismatch start', 'mismatch end'), axis=1)
    sum_missmatch_border_distal = result[cols].apply(lambda row: sum_of_size_in_border_line(row, 'distal border line', 'mismatch type', 'mismatch size', 'mismatch start', 'mismatch end'), axis=1)

    cols = ["mir type", 'bulge type', 'bulge size', 'bulge start', 'bulge end']
    sum_bulge_border_proximal = result[cols].apply(lambda row: sum_of_size_in_border_line(row, 'proximal border line', 'bulge type', 'bulge size', 'bulge start', 'bulge end'), axis=1)
    sum_bulge_border_distal = result[cols].apply(lambda row: sum_of_size_in_border_line(row, 'distal border line', 'bulge type', 'bulge size', 'bulge start', 'bulge end'), axis=1)

    cols = ["mir type", 'internal type', 'internal loop HSBL', 'internal start', 'internal end']
    sum_internal_border_proximal = result[cols].apply(lambda row: sum_of_size_in_border_line(row, 'proximal border line', 'internal type', 'internal loop HSBL', 'internal start', 'internal end'), axis=1)
    sum_internal_border_distal = result[cols].apply(lambda row: sum_of_size_in_border_line(row, 'distal border line', 'internal type', 'internal loop HSBL', 'internal start', 'internal end'), axis=1)

    cols = ["hit start", "hit end", "psep", "mir type"]
    sum_of_residue = result[cols].apply(lambda row: number_of_residue(row), axis=1)
    
    cols = ["bulge type", "bulge strand", "bulge size"]
    sum_bulge_zero = result[cols].apply(lambda row: sum_of_size_in_hit_only_zero(row), axis=1)

    result['sum of residue in terminal loop'] = sum_of_residue

    _sum = sum_bulge + sum_internal + sum_bulge_border_proximal + sum_bulge_border_distal + sum_internal_border_proximal + sum_internal_border_distal + sum_of_residue
    result["ACCEPTABLE_NUM_FOR_HIT_LOCATIONS_IN_BULGES_OR_LOOPS"] = _sum

    result["ACCEPTABLE_NUM_FOR_UNMATCHED_LOCATIONS_IN_HIT_REGION * 2"] = (_sum + (sum_missmatch + sum_missmatch_border_proximal + sum_missmatch_border_distal) * 2)
    result["ACCEPTABLE_NUM_FOR_UNMATCHED_LOCATIONS_IN_HIT_REGION"] = (_sum + (sum_missmatch + sum_missmatch_border_proximal + sum_missmatch_border_distal))

    _sum = sum_bulge_zero + sum_internal_hsbl + sum_bulge_border_proximal + sum_bulge_border_distal + sum_internal_border_proximal + sum_internal_border_distal + sum_of_residue
    result["ACCEPTABLE_NUM_FOR_HIT_LOCATIONS_IN_BULGES_OR_LOOPS_mayers"] = _sum

    result["ACCEPTABLE_NUM_FOR_UNMATCHED_LOCATIONS_IN_HIT_REGION_mayers * 2"] = (_sum + (sum_missmatch + sum_missmatch_border_proximal + sum_missmatch_border_distal) * 2)
    result["ACCEPTABLE_NUM_FOR_UNMATCHED_LOCATIONS_IN_HIT_REGION_mayers"] = (_sum + (sum_missmatch + sum_missmatch_border_proximal + sum_missmatch_border_distal))

    result["TOTAL_NUM_OF_MISMACHED_POSITIONS"] = sum_missmatch

    result["TOTAL_NUM_OF_NONMATCHING_POSITIONS"] = (sum_missmatch + sum_bulge + sum_internal)

    result["TOTAL_NUM_OF_POSITIONS_IN_BULGES_AND_LOOPS"] = (sum_bulge + sum_internal)
    
    result["MATURE_DUPLEX_INVOLVEMENT_IN_APICAL_LOOP"] = result.apply(lambda row: check_involvement(row), axis=1)

    cols = ['hit len',
        'flanking GC content',
        'flanking MFEI',
        'hit GC content',        
        'hit complementarity percentage',        
        'num of linking residues',
        'boi GC content',
        'boi delta G',
        'boi AMFE',
        'boi MFEI',        
        'branch#1 apical loop size',        
        'precursor gc content',
        'precursor delta G',
        'precursor AMFE',
        'precursor MFEI',
        'number of terminal structures',
        'primary stem length',
        'mismatch',
        'bulge',
        'internal loop',
        'base structure corrected length',
        'primary stem corrected length',
        'Loop distal junction distance',
        'Loop proximal junction distance',
        'sum of residue in terminal loop',
        'ACCEPTABLE_NUM_FOR_HIT_LOCATIONS_IN_BULGES_OR_LOOPS',
        'ACCEPTABLE_NUM_FOR_UNMATCHED_LOCATIONS_IN_HIT_REGION * 2',
        'ACCEPTABLE_NUM_FOR_UNMATCHED_LOCATIONS_IN_HIT_REGION',
        'ACCEPTABLE_NUM_FOR_HIT_LOCATIONS_IN_BULGES_OR_LOOPS_mayers',
        'ACCEPTABLE_NUM_FOR_UNMATCHED_LOCATIONS_IN_HIT_REGION_mayers * 2',
        'ACCEPTABLE_NUM_FOR_UNMATCHED_LOCATIONS_IN_HIT_REGION_mayers',
        'TOTAL_NUM_OF_MISMACHED_POSITIONS',
        'TOTAL_NUM_OF_NONMATCHING_POSITIONS',
        'TOTAL_NUM_OF_POSITIONS_IN_BULGES_AND_LOOPS',
        'MATURE_DUPLEX_INVOLVEMENT_IN_APICAL_LOOP',
        ]

    for region in ["loop distal", "hit region", "distal border line", "loop proximal", "proximal border line"]:
        _col = ["mismatch type", "mismatch size"]
        result[f'max mismatch in {region}'] = result[_col].apply(lambda row: get_max_in_region(row, "mismatch type", "mismatch size", region), axis=1)
        cols.append(f'max mismatch in {region}')        
        result[f'number mismatch in {region}'] = result[_col].apply(lambda row: get_number_in_region(row, "mismatch type", "mismatch size", region), axis=1)
        cols.append(f'number mismatch in {region}')        
        result[f'sum mismatch in {region}'] = result[_col].apply(lambda row: get_sum_in_region(row, "mismatch type", "mismatch size", region), axis=1)    
        cols.append(f'sum mismatch in {region}')
        _col = ["bulge type", "bulge size"]
        result[f'max bulge in {region}'] = result[_col].apply(lambda row: get_max_in_region(row, "bulge type", "bulge size", region), axis=1)
        cols.append(f'max bulge in {region}')        
        result[f'number bulge in {region}'] = result[_col].apply(lambda row: get_number_in_region(row, "bulge type", "bulge size", region), axis=1)
        cols.append(f'number bulge in {region}')        
        result[f'sum bulge in {region}'] = result[_col].apply(lambda row: get_sum_in_region(row, "bulge type", "bulge size", region), axis=1)
        cols.append(f'sum bulge in {region}')
        _col = ["internal type", "internal loop total size"]
        result[f'max loop in {region}'] = result[_col].apply(lambda row: get_max_in_region(row, "internal type", "internal loop total size", region), axis=1)
        cols.append(f'max loop in {region}')        
        result[f'number loop in {region}'] = result[_col].apply(lambda row: get_number_in_region(row, "internal type", "internal loop total size", region), axis=1)
        cols.append(f'number loop in {region}')        
        result[f'sum loop in {region}'] = result[_col].apply(lambda row: get_sum_in_region(row, "internal type", "internal loop total size", region), axis=1)
        cols.append(f'sum loop in {region}')
        
    X = result[cols]    

    # replace inf with max
    for c in X.columns:
        m = X[X[c] != np.inf][c].max()
        X[c].replace([np.inf], m, inplace=True)    
    # standardization    
    mu = X.mean()
    std = X.std()
    X = (X - mu) / std    
    # mir type
    encoder = OneHotEncoder(handle_unknown='ignore')
    encoder_df = pd.DataFrame(encoder.fit_transform(result[['mir type']]).toarray())
    X = X.join(encoder_df)
    # TNF    
    _tnf = pd.DataFrame(df['hit seq'].apply(lambda x: tnf_calc(x)))
    X = X.join(_tnf)
    # start end nucleotide
    for c in ["-3", "-2", "-1", "", "+1", "+2"]:
        col_name = f"hit start{c} composition"
        _nuc = pd.DataFrame(df[col_name].apply(lambda x: com_one_hot(x)))
        _nuc.columns = [f"start_{c}_{i}" for i in _nuc.columns]
        X = X.join(_nuc)
        colname = f"hit end{c} composition"
        _nuc = pd.DataFrame(df[col_name].apply(lambda x: com_one_hot(x)))
        _nuc.columns = [f"end_{c}_{i}" for i in _nuc.columns ]
        X = X.join(_nuc)    
    # connectivity
    cols = []
    for c in  ["-3", "-2", "-1", "", "+1", "+2"]:
        cols.append(f"connectivity hit start{c}")
        cols.append(f"connectivity hit end{c}")
    #for c in range(2,13):
    #    cols.append(f"seed connectivity{c}")
    X = X.join(result[cols])    
    X = X.astype("float32")    
    return X


def get_target(df, reference):
    Y = df['hit seq'].isin(reference['data'])
    Y = Y.apply(lambda x: 1 if x else 0)
    Y = np_utils.to_categorical(Y, 2)
    Y = Y.astype("float32")
    return Y

