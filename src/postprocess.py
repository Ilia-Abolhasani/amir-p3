import json
import tqdm
import numpy as np
import pandas as pd
from read_configs import DotDict, read_titles, read_erros

titles = read_titles()

def _is_allowed_clear(row, limit):
    if float(row[titles.distal_junc_dist]) >= limit:
        return True
    if float(row[titles.prx_junc_dist]) >= limit:
        return True
    return False


def _check_border_line(row, type_str, size_str, limit):
    valid = True
    mismatch_type = row[type_str]
    for i in range(len(mismatch_type)):
        if (
            mismatch_type[i] == "distal border line"
            or mismatch_type[i] == "proximal border line"
        ):
            if row[size_str][i] > limit:
                valid = False
    return valid



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
                if mir_type == "5p":
                    _sum += _size[i] - _start[i]
                if mir_type == "3p":
                    _sum += _size[i] - _end[i]
            if border_type == "proximal border line":
                if mir_type == "5p":
                    _sum += _size[i] - _end[i]
                if mir_type == "3p":
                    _sum += _size[i] - _start[i]
    return _sum
                                                              
            
def postprocess(df, config=None):
    level2 = df.copy()
    if config is None:  # read from filter_level2.json
        with open("./src/config/filter_level2.json") as json_file:
            config = json.load(json_file)
    config = DotDict(config)                           
    level2 = level2[level2["hit len"] >= config.hit_len_min]
    level2 = level2[level2["hit len"] <= config.hit_len_max]    
    level2 = level2[
        level2["number of terminal structures"]
        >= config.number_of_terminal_structure_min
    ]
    level2 = level2[
        level2["number of terminal structures"]
        <= config.number_of_terminal_structure_max
    ]            
    level2 = level2[level2["boi GC content"] >= config.boi_gc_content_min]
    level2 = level2[level2["boi GC content"] <= config.boi_gc_content_max]
                
    cols = [titles.distal_junc_dist, titles.prx_junc_dist]
    level2 = level2[
        level2[cols].apply(
            lambda row: _is_allowed_clear(row, config.minimum_required_clear_region),
            axis=1,
        )
    ]

    cols = ["mismatch type", "mismatch size"]
    level2 = level2[
        level2[cols].apply(
            lambda row: _check_border_line(
                row, "mismatch type", "mismatch size", config.border_line_mismatch_max
            ),
            axis=1,
        )
    ]

    cols = ["bulge type", "bulge size"]
    level2 = level2[
        level2[cols].apply(
            lambda row: _check_border_line(
                row, "bulge type", "bulge size", config.border_line_bulge_max
            ),
            axis=1,
        )
    ]

    cols = ["internal type", "internal loop total size"]
    level2 = level2[
        level2[cols].apply(
            lambda row: _check_border_line(
                row,
                "internal type",
                "internal loop total size",
                config.border_line_internal_max,
            ),
            axis=1,
        )
    ]
    

    cols = [
        titles.mir_type,
        "mismatch type",
        "mismatch size",
        "mismatch start",
        "mismatch end",
    ]
    sum_missmatch_border_proximal = level2[cols].apply(
        lambda row: _sum_of_size_in_border_line(
            row,
            "proximal border line",
            "mismatch type",
            "mismatch size",
            "mismatch start",
            "mismatch end",
        ),
        axis=1,
    )
    sum_missmatch_border_distal = level2[cols].apply(
        lambda row: _sum_of_size_in_border_line(
            row,
            "distal border line",
            "mismatch type",
            "mismatch size",
            "mismatch start",
            "mismatch end",
        ),
        axis=1,
    )

    cols = [titles.mir_type, "bulge type", "bulge size", "bulge start", "bulge end"]
    sum_bulge_border_proximal = level2[cols].apply(
        lambda row: _sum_of_size_in_border_line(
            row,
            "proximal border line",
            "bulge type",
            "bulge size",
            "bulge start",
            "bulge end",
        ),
        axis=1,
    )
    sum_bulge_border_distal = level2[cols].apply(
        lambda row: _sum_of_size_in_border_line(
            row,
            "distal border line",
            "bulge type",
            "bulge size",
            "bulge start",
            "bulge end",
        ),
        axis=1,
    )

    cols = [
        titles.mir_type,
        "internal type",
        "internal loop HSBL",
        "internal start",
        "internal end",
    ]
    sum_internal_border_proximal = level2[cols].apply(
        lambda row: _sum_of_size_in_border_line(
            row,
            "proximal border line",
            "internal type",
            "internal loop HSBL",
            "internal start",
            "internal end",
        ),
        axis=1,
    )
    sum_internal_border_distal = level2[cols].apply(
        lambda row: _sum_of_size_in_border_line(
            row,
            "distal border line",
            "internal type",
            "internal loop HSBL",
            "internal start",
            "internal end",
        ),
        axis=1,
    )               

    border_proximal = (
        sum_missmatch_border_proximal
        + sum_bulge_border_proximal
        + sum_internal_border_proximal
    )[level2.index] == 0
    border_distal = (
        sum_missmatch_border_distal
        + sum_bulge_border_distal
        + sum_internal_border_distal
    )[level2.index] == 0
    if config.border_line_structure_allowance == "NOT ACCEPTED":
        level2 = level2[border_proximal & border_distal]
    if config.border_line_structure_allowance == "1 END ONLY":
        level2 = level2[border_proximal | border_distal]
    return df.index[~df.index.isin(level2.index)]