import numpy as np
import tqdm
import pandas as pd
import sys


# parameter
DELTA_G_MIN = -1000
DELTA_G_MAX = 0
HIT_LEN_MIN = 16
HIT_LEN_MAX = 25
HIT_COMPLEMENTARITY_PERCENTAGE_MIN = 0.3
HIT_COMPLEMENTARITY_PERCENTAGE_MAX = 1
NUMBER_OF_TERMINAL_STRUCTURE_MIN = 1
NUMBER_OF_TERMINAL_STRUCTURE_MAX = 5
BOI_GC_CONTENT_MIN = 0.01
BOI_GC_CONTENT_MAX = 1
BORDER_LINE_MISMATCH_MAX = 3
BORDER_LINE_BULGE_MAX = 3
BORDER_LINE_INTERNAL_MAX = 3
TOTAL_NUM_OF_NONMATCHING_POSITIONS = 5
TOTAL_NUM_OF_MISMACHED_POSITIONS = 4
TOTAL_NUM_OF_POSITIONS_IN_BULGES_AND_LOOPS = 3
MAX_ALLOWED_BULGE_SIZE_IN_HIT_REGION = 2
MAX_ALLOWED_INTERNAL_LOOP_SIZE_IN_HIT_REGION = 3
MINIMUM_REQUIRED_CLEAR_REGION = 8
ACCEPTABLE_NUM_FOR_HIT_LOCATIONS_IN_BULGES_OR_LOOPS = 3
ACCEPTABLE_NUM_FOR_UNMATCHED_LOCATIONS_IN_HIT_REGION = 5

MIN_NUM_OF_LINKING_RESIDUES = 0
MAX_NUM_OF_LINKING_RESIDUES = 1000
MIN_HIT_GC_CONTENT_PERCENTAGE = 25
MAX_HIT_GC_CONTENT_PERCENTAGE = 75
DELETE_IF_MATURE_DUPLEX_INVOLVEMENT_IN_APICAL_LOOP = "YES"

def is_allowed(row, type_str, size_str, limmit):
    mismatch_type = row[type_str]
    for i in range(len(mismatch_type)):
        if mismatch_type[i] == "hit region":
            if row[size_str][i] > limmit:
                return False
    return True


def is_allowed_clear(row):
    row['Loop distal junction distance'] = float(row['Loop distal junction distance'])
    row['Loop proximal junction distance'] = float(row['Loop proximal junction distance'])
    if(row['Loop distal junction distance'] >= MINIMUM_REQUIRED_CLEAR_REGION or
            row['Loop proximal junction distance'] >= MINIMUM_REQUIRED_CLEAR_REGION):
        return True
    return False


def check_border_line(row,type_str, size_str, limmit):
    valid = True
    mismatch_type = row[type_str]
    for i in range(len(mismatch_type)):
        if mismatch_type[i] == "distal border line" or mismatch_type[i] == "proximal border line":
            if row[size_str][i] > limmit:
                valid = False
    return valid


def sum_of_size_in_hit(row, type_str, size_str):
    _sum = 0
    mismatch_type = row[type_str]
    for i in range(len(mismatch_type)):
        if mismatch_type[i] == "hit region":
            _sum += row[size_str][i]
    return _sum


def number_of_residue(row):
    hit_end = row['hit end']
    hit_end = float(hit_end)
    hit_start = row['hit start']
    psep = row['psep']
    psep = float(psep)
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
    if(row['number of terminal structures'] > 1):
        return True
    start = row['branch#1 apical loop start']
    end = row['branch#1 apical loop end']
    for col in ['hit start', 'hit end', 'star start', 'star end']:
        if(start < row[col] < end ):
            return False                    
    return True


def convert(row):
    row['delta G'] = float(row['delta G'])
    row['number of terminal structures'] = float(row['number of terminal structures'])
    row['num of linking residues'] = float(row['num of linking residues'])
    row['boi GC content'] = float(row['boi GC content'])
    for item in ['mismatch type', 'mismatch size', 'mismatch start', 'mismatch end',
                 'bulge type', 'bulge size', 'bulge start', 'bulge end',
                 'internal type', 'internal loop total size', 'internal start', 'internal end']:
        row[item] = eval(row[item])
    return row


header = True
for chunk in tqdm.tqdm(pd.read_csv("../Result/result_level1_filter.csv", chunksize=10**4)):
    level2 = chunk.apply(lambda row: convert(row), axis=1)
    level2 = level2[level2['delta G'] >= DELTA_G_MIN]
    level2 = level2[level2['delta G'] <= DELTA_G_MAX]
    level2 = level2[level2['hit len'] >= HIT_LEN_MIN]
    level2 = level2[level2['hit len'] <= HIT_LEN_MAX]
    level2 = level2[level2['hit complementarity percentage'] >= HIT_COMPLEMENTARITY_PERCENTAGE_MIN]
    level2 = level2[level2['hit complementarity percentage'] <= HIT_COMPLEMENTARITY_PERCENTAGE_MAX]
    level2 = level2[level2['number of terminal structures'] >= NUMBER_OF_TERMINAL_STRUCTURE_MIN]
    level2 = level2[level2['number of terminal structures'] <= NUMBER_OF_TERMINAL_STRUCTURE_MAX]
    level2 = level2[level2['boi GC content'] >= BOI_GC_CONTENT_MIN]
    level2 = level2[level2['boi GC content'] <= BOI_GC_CONTENT_MAX]
    level2 = level2[level2['num of linking residues'] >= MIN_NUM_OF_LINKING_RESIDUES]
    level2 = level2[level2['num of linking residues'] <= MAX_NUM_OF_LINKING_RESIDUES]
    level2 = level2[level2['hit GC content'] >= MIN_HIT_GC_CONTENT_PERCENTAGE]
    level2 = level2[level2['hit GC content'] <= MAX_HIT_GC_CONTENT_PERCENTAGE]

    level2 = level2[level2.apply(lambda row: is_allowed(row, "bulge type", "bulge size", MAX_ALLOWED_BULGE_SIZE_IN_HIT_REGION), axis=1)]
    level2 = level2[level2.apply(lambda row: is_allowed(row, "internal type", "internal loop total size", MAX_ALLOWED_INTERNAL_LOOP_SIZE_IN_HIT_REGION), axis=1)]
    level2 = level2[level2.apply(lambda row: is_allowed_clear(row), axis=1)]
    level2 = level2[level2.apply(lambda row: check_border_line(row, "mismatch type", "mismatch size", BORDER_LINE_MISMATCH_MAX), axis=1)]
    level2 = level2[level2.apply(lambda row: check_border_line(row, "bulge type", "bulge size", BORDER_LINE_BULGE_MAX), axis=1)]
    level2 = level2[level2.apply(lambda row: check_border_line(row, "internal type", "internal loop total size", BORDER_LINE_INTERNAL_MAX), axis=1)]

    sum_missmatch = level2.apply(lambda row: sum_of_size_in_hit(row, 'mismatch type', 'mismatch size'), axis=1)
    sum_bulge = level2.apply(lambda row: sum_of_size_in_hit(row, 'bulge type', 'bulge size'), axis=1)
    sum_internal = level2.apply(lambda row: sum_of_size_in_hit(row, 'internal type', 'internal loop total size'), axis=1)
    sum_missmatch_border_proximal = level2.apply(lambda row: sum_of_size_in_border_line(row, 'proximal border line', 'mismatch type', 'mismatch size', 'mismatch start', 'mismatch end'), axis=1)
    sum_missmatch_border_distal = level2.apply(lambda row: sum_of_size_in_border_line(row, 'distal border line', 'mismatch type', 'mismatch size', 'mismatch start', 'mismatch end'), axis=1)
    sum_bulge_border_proximal = level2.apply(lambda row: sum_of_size_in_border_line(row, 'proximal border line', 'bulge type', 'bulge size', 'bulge start', 'bulge end'), axis=1)
    sum_bulge_border_distal = level2.apply(lambda row: sum_of_size_in_border_line(row, 'distal border line', 'bulge type', 'bulge size', 'bulge start', 'bulge end'), axis=1)
    sum_internal_border_proximal = level2.apply(lambda row: sum_of_size_in_border_line(row, 'proximal border line', 'internal type', 'internal loop HSBL', 'internal start', 'internal end'), axis=1)
    sum_internal_border_distal = level2.apply(lambda row: sum_of_size_in_border_line(row, 'distal border line', 'internal type', 'internal loop HSBL', 'internal start', 'internal end'), axis=1)
    sum_of_residue = level2.apply(lambda row: number_of_residue(row), axis=1)

    _sum = sum_bulge + sum_internal + sum_bulge_border_proximal + sum_bulge_border_distal + sum_internal_border_proximal + sum_internal_border_distal + sum_of_residue
    level2 = level2[_sum <= ACCEPTABLE_NUM_FOR_HIT_LOCATIONS_IN_BULGES_OR_LOOPS]

    _sum = (_sum + (sum_missmatch + sum_missmatch_border_proximal + sum_missmatch_border_distal) * 2)[level2.index]
    level2 = level2[_sum <= ACCEPTABLE_NUM_FOR_UNMATCHED_LOCATIONS_IN_HIT_REGION]

    _sum = (sum_missmatch + sum_bulge + sum_internal)[level2.index]
    level2 = level2[_sum <= TOTAL_NUM_OF_NONMATCHING_POSITIONS]

    _sum = sum_missmatch[level2.index]
    level2 = level2[_sum <= TOTAL_NUM_OF_MISMACHED_POSITIONS]

    _sum = (sum_bulge + sum_internal)[level2.index]
    level2 = level2[_sum <= TOTAL_NUM_OF_POSITIONS_IN_BULGES_AND_LOOPS]

    if(DELETE_IF_MATURE_DUPLEX_INVOLVEMENT_IN_APICAL_LOOP == "YES"):
        level2 = level2[level2.apply(lambda row: check_involvement(row), axis=1)]

    level2.to_csv("../Result/result_level2_filter.csv", header=header, mode='a', index=False)
    header = False