import sys
import json
import tqdm
import numpy as np
import pandas as pd


# parameter
class DotDict(dict):
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

with open('./src/config/titles.json') as json_file:
    titles = DotDict(json.load(json_file))

def _is_allowed(row, type_str, size_str, limit):
    mismatch_type = row[type_str]
    for i in range(len(mismatch_type)):
        if mismatch_type[i] == "hit region":
            if row[size_str][i] > limit:
                return False
    return True


def _is_allowed_clear(row, limit):
    if float(row[titles.distal_junc_dist]) >= limit:
        return True
    if float(row[titles.prx_junc_dist]) >= limit:
        return True
    return False


def _check_border_line(row,type_str, size_str, limit):
    valid = True
    mismatch_type = row[type_str]
    for i in range(len(mismatch_type)):
        if mismatch_type[i] == "distal border line" or mismatch_type[i] == "proximal border line":
            if row[size_str][i] > limit:
                valid = False
    return valid


def _sum_of_size_in_hit(row, type_str, size_str):
    _sum = 0
    mismatch_type = row[type_str]
    for i in range(len(mismatch_type)):
        if mismatch_type[i] == "hit region":
            _sum += row[size_str][i]
    return _sum


def _sum_of_size_in_hit_only_zero(row):
    _sum = 0
    bulge_type = row[titles.bulge_type]
    bulge_strand = row[titles.bulge_strand]
    for i in range(len(bulge_type)):
        if bulge_type[i] == "hit region" and bulge_strand[i] == "zero":
            _sum += row[titles.bulge_size][i]
    return _sum


def _number_of_residue(row):
    hit_end = float(row[titles.hit_end])
    hit_start = row[titles.hit_start]
    psep = row[titles.psep]
    psep = float(psep)
    mir_type = row[titles.mir_type]
    if mir_type == '5p':
        if psep < hit_end:
            return hit_end - psep
    if mir_type == '3p':
        if psep > hit_start:
            return psep - hit_start
    return 0


def _sum_of_size_in_border_line(row, border_type, type_str, size_str, start, end):
    _sum = 0
    _size = row[size_str]
    _start = row[start]
    _end = row[end]
    mir_type = row[titles.mir_type]
    mismatch_type = row[type_str]
    for i in range(len(mismatch_type)):
        if mismatch_type[i] == border_type:
            if border_type == "distal border line":
                if mir_type == '5p':
                    _sum += _size[i] - _start[i]
                if mir_type == '3p':
                    _sum += _size[i] - _end[i]
            if border_type == "proximal border line":
                if mir_type == '5p':
                    _sum += _size[i] - _end[i]
                if mir_type == '3p':
                    _sum += _size[i] - _start[i]
    return _sum


def _check_involvement(row):
    if row[titles.n_term_struc] > 1:
        return True
    start = row['branch#1 apical loop start']
    end = row['branch#1 apical loop end']
    for col in [titles.hit_start, titles.hit_end, titles.star_start, titles.star_end]:
        if start < row[col] < end:
            return False                    
    return True


def _loopHSBL_SSBL_check(row, limit):
    loopType = row['internal type']
    hsbl = row['internal loop HSBL']
    ssbl = row['internal loop SSBL']
    for i in range(len(loopType)):
        if loopType[i] != "loop proximal" and loopType[i] != "loop distal":
            if hsbl[i] > limit:
                return False
            if ssbl[i] > limit:
                return False
    return True


def convert(row):
    row[titles.dg] = float(row[titles.dg])
    row[titles.n_term_struc] = float(row[titles.n_term_struc])
    row[titles.num_of_lnk_res] = float(row[titles.num_of_lnk_res])
    row[titles.boi_gc] = float(row[titles.boi_gc])
    if row[titles.pre_mfei] == "-":
        row[titles.pre_mfei] = 10*10
    row[titles.pre_mfei] = float(row[titles.pre_mfei])
    for item in ['mismatch type', 'mismatch size', 'mismatch start', 'mismatch end',
                 'bulge type', 'bulge size', 'bulge start', 'bulge end', "bulge strand",
                 'internal type', 'internal loop total size', 'internal start', 'internal end', 'internal loop HSBL', 'internal loop SSBL']:
        row[item] = eval(row[item])
    return row


def filter2_run(inp, config):
    level2 = inp
    level2 = level2[level2[titles.dg] >= config.delta_g_min]
    level2 = level2[level2[titles.dg] <= config.delta_g_max]
    level2 = level2[level2['hit len'] >= config.hit_len_min]
    level2 = level2[level2['hit len'] <= config.hit_len_max]
    level2 = level2[level2['hit complementarity percentage'] >= config.hit_complementarity_percentage_min]
    level2 = level2[level2['hit complementarity percentage'] <= config.hit_complementarity_percentage_max]
    level2 = level2[level2['number of terminal structures'] >= config.number_of_terminal_structure_min]
    level2 = level2[level2['number of terminal structures'] <= config.number_of_terminal_structure_max]
    level2['hit GC content'] = level2['hit GC content']
    level2['boi GC content'] = level2['boi GC content']
    level2 = level2[level2['boi GC content'] >= config.boi_gc_content_min]
    level2 = level2[level2['boi GC content'] <= config.boi_gc_content_max]
    level2 = level2[level2['num of linking residues'] >= config.num_of_linking_residues_min]
    level2 = level2[level2['num of linking residues'] <= config.num_of_linking_residues_max]

    level2 = level2[level2['hit GC content'] >= config.hit_gc_content_percentage_min]
    level2 = level2[level2['hit GC content'] <= config.hit_gc_content_percentage_max]

    level2 = level2[level2['precursor MFEI'] >= config.precursor_mfei_min]
    level2 = level2[level2['precursor MFEI'] <= config.precursor_mfei_max]

    cols = ["mismatch type", "mismatch size"]
    level2 = level2[level2[cols].apply(lambda row: _is_allowed(row, "mismatch type", "mismatch size", config.max_allowed_mismatch_size_in_hit_region),axis=1)]

    cols = ["bulge type", "bulge size"]
    level2 = level2[level2[cols].apply(lambda row: _is_allowed(row, "bulge type", "bulge size", config.max_allowed_bulge_size_in_hit_region), axis=1)]

    cols = ["internal type", "internal loop total size"]
    level2 = level2[level2[cols].apply(lambda row: _is_allowed(row, "internal type", "internal loop total size", config.max_allowed_internal_loop_size_in_hit_region), axis=1)]

    cols = [titles.distal_junc_dist, titles.prx_junc_dist]
    level2 = level2[level2[cols].apply(lambda row: _is_allowed_clear(row, config.minimum_required_clear_region), axis=1)]

    cols = ["mismatch type", "mismatch size"]
    level2 = level2[level2[cols].apply(lambda row: _check_border_line(row, "mismatch type", "mismatch size", config.border_line_mismatch_max), axis=1)]

    cols = ["bulge type", "bulge size"]
    level2 = level2[level2[cols].apply(lambda row: _check_border_line(row, "bulge type", "bulge size", config.border_line_bulge_max), axis=1)]

    cols = ["internal type", "internal loop total size"]
    level2 = level2[level2[cols].apply(lambda row: _check_border_line(row, "internal type", "internal loop total size", config.border_line_internal_max), axis=1)]

    cols = ["internal type", "internal loop HSBL", "internal loop SSBL"]
    level2 = level2[level2[cols].apply(lambda row: _loopHSBL_SSBL_check(row, config.max_allowed_hsbl_ssbl_size), axis=1)]

    cols = ['mismatch type', 'mismatch size']
    sum_missmatch = level2[cols].apply(lambda row: _sum_of_size_in_hit(row, 'mismatch type', 'mismatch size'), axis=1)

    cols = [ 'bulge type', 'bulge size']
    sum_bulge = level2[cols].apply(lambda row: _sum_of_size_in_hit(row, 'bulge type', 'bulge size'), axis=1)

    cols = ['internal type', 'internal loop total size']
    sum_internal = level2[cols].apply(lambda row: _sum_of_size_in_hit(row, 'internal type', 'internal loop total size'), axis=1)

    cols = ['internal type', 'internal loop HSBL']
    sum_internal_hsbl = level2[cols].apply(lambda row: _sum_of_size_in_hit(row, 'internal type', 'internal loop HSBL'),axis=1)

    cols = [titles.mir_type, 'mismatch type', 'mismatch size', 'mismatch start', 'mismatch end']
    sum_missmatch_border_proximal = level2[cols].apply(lambda row: _sum_of_size_in_border_line(row, 'proximal border line', 'mismatch type', 'mismatch size', 'mismatch start', 'mismatch end'), axis=1)
    sum_missmatch_border_distal = level2[cols].apply(lambda row: _sum_of_size_in_border_line(row, 'distal border line', 'mismatch type', 'mismatch size', 'mismatch start', 'mismatch end'), axis=1)

    cols = [titles.mir_type, 'bulge type', 'bulge size', 'bulge start', 'bulge end']
    sum_bulge_border_proximal = level2[cols].apply(lambda row: _sum_of_size_in_border_line(row, 'proximal border line', 'bulge type', 'bulge size', 'bulge start', 'bulge end'), axis=1)
    sum_bulge_border_distal = level2[cols].apply(lambda row: _sum_of_size_in_border_line(row, 'distal border line', 'bulge type', 'bulge size', 'bulge start', 'bulge end'), axis=1)

    cols = [titles.mir_type, 'internal type', 'internal loop HSBL', 'internal start', 'internal end']
    sum_internal_border_proximal = level2[cols].apply(lambda row: _sum_of_size_in_border_line(row, 'proximal border line', 'internal type', 'internal loop HSBL', 'internal start', 'internal end'), axis=1)
    sum_internal_border_distal = level2[cols].apply(lambda row: _sum_of_size_in_border_line(row, 'distal border line', 'internal type', 'internal loop HSBL', 'internal start', 'internal end'), axis=1)

    cols = [titles.bulge_type, titles.bulge_strand, titles.bulge_size]
    sum_bulge_zero = level2[cols].apply(lambda row: _sum_of_size_in_hit_only_zero(row), axis=1)

    cols = [titles.hit_start, titles.hit_end, titles.psep, titles.mir_type]
    sum_of_residue = level2[cols].apply(lambda row: _number_of_residue(row), axis=1)

    _sum = sum_bulge_zero + sum_internal_hsbl + sum_bulge_border_proximal + sum_bulge_border_distal + sum_internal_border_proximal + sum_internal_border_distal + sum_of_residue
    level2 = level2[_sum <= config.acceptable_num_for_hit_locations_in_bulges_or_loops]

    _sum = (_sum + (sum_missmatch + sum_missmatch_border_proximal + sum_missmatch_border_distal) * 1)[level2.index] ###### Changed from 2 to 1
    level2 = level2[_sum <= config.acceptable_num_for_unmatched_locations_in_hit_region]

    _sum = (sum_missmatch + sum_bulge + sum_internal)[level2.index]
    level2 = level2[_sum <= config.total_num_of_nonmatching_positions]

    _sum = sum_missmatch[level2.index]
    level2 = level2[_sum <= config.total_num_of_mismached_positions]

    _sum = (sum_bulge + sum_internal)[level2.index]
    level2 = level2[_sum <= config.total_num_of_positions_in_bulges_and_loops]

    if config.delete_if_mature_duplex_involvement_in_apical_loop == "YES":
        cols = [titles.n_term_struc, titles.hit_start, titles.hit_end, titles.star_start, titles.star_end,
                'branch#1 apical loop start', 'branch#1 apical loop end']
        level2 = level2[level2[cols].apply(lambda row: _check_involvement(row), axis=1)]

    border_proximal = (sum_missmatch_border_proximal + sum_bulge_border_proximal + sum_internal_border_proximal)[level2.index] == 0
    border_distal = (sum_missmatch_border_distal + sum_bulge_border_distal + sum_internal_border_distal)[level2.index] == 0
    if config.border_line_structure_allowance == "NOT ACCEPTED":
        level2 = level2[border_proximal & border_distal]
    if config.border_line_structure_allowance == "1 END ONLY":
        level2 = level2[border_proximal | border_distal]
    return level2


#BORDER_LINE_STRUCTURE_ALLOWANCE = "1 END ONLY"  # "NOT ACCEPTED"   "1 END ONLY" "BOTH END"
def filter2(input_file, output_file, config=None, chunksize=10 ** 5):
    if config is None:  # read from filter_level2.json
        with open('./src/config/filter_level2.json') as json_file:
            config = json.load(json_file)
    config = DotDict(config)

    header = True
    for chunk in tqdm.tqdm(pd.read_csv(input_file, chunksize=chunksize)):
        chunk = chunk.apply(lambda row: convert(row), axis=1)
        level2 = filter2_run(chunk, config)
        level2.to_csv(output_file, header=header, mode='a', index=False)
        header = False
