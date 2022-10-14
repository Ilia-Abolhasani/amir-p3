import os
import json
import tqdm
import numpy as np
import pandas as pd
from read_configs import DotDict, read_titles, read_erros
from convertor import convert

titles = read_titles()


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


def _check_involvement(row):
    if row[titles.n_term_struc] > 1:
        return True
    start = row["branch#1 apical loop start"]
    end = row["branch#1 apical loop end"]
    for col in [titles.hit_start, titles.hit_end, titles.star_start, titles.star_end]:
        if start < row[col] < end:
            return False
    return True


def _loopHSBL_SSBL_check(row, limit):
    loopType = row["internal type"]
    hsbl = row["internal loop HSBL"]
    ssbl = row["internal loop SSBL"]
    for i in range(len(loopType)):
        if loopType[i] != "loop proximal" and loopType[i] != "loop distal":
            if hsbl[i] > limit:
                return False
            if ssbl[i] > limit:
                return False
    return True


def _postprocess_run(inp, config):
    level2 = inp
    if (config.delta_g_min in config):
        level2 = level2[level2[titles.dg] >= config.delta_g_min]
    if (config.delta_g_max in config):
        level2 = level2[level2[titles.dg] <= config.delta_g_max]
    if (config.hit_len_min in config):
        level2 = level2[level2[titles.hit_len] >= config.hit_len_min]
    if (config.hit_len_max in config):
        level2 = level2[level2[titles.hit_len] <= config.hit_len_max]
    if (config.hit_complementarity_percentage_min in config):
        level2 = level2[level2[titles.hit_cmplt_pct]
                        >= config.hit_complementarity_percentage_min]
    if (config.hit_complementarity_percentage_max in config):
        level2 = level2[level2[titles.hit_cmplt_pct] <=
                        config.hit_complementarity_percentage_max]
    if (config.number_of_terminal_structure_min in config):
        level2 = level2[level2[titles.n_term_struc]
                        >= config.number_of_terminal_structure_min]
    if (config.number_of_terminal_structure_max in config):
        level2 = level2[level2[titles.n_term_struc]
                        <= config.number_of_terminal_structure_max]

    level2[titles.hit_gc_content] = level2[titles.hit_gc_content]
    level2[titles.boi_gc] = level2[titles.boi_gc]
    if (config.boi_gc_content_min in config):
        level2 = level2[level2[titles.boi_gc] >= config.boi_gc_content_min]
    if (config.boi_gc_content_max in config):
        level2 = level2[level2[titles.boi_gc] <= config.boi_gc_content_max]
    if (config.num_of_linking_residues_min in config):
        level2 = level2[level2[titles.num_of_lnk_res]
                        >= config.num_of_linking_residues_min]
    if (config.num_of_linking_residues_max in config):
        level2 = level2[level2[titles.num_of_lnk_res]
                        <= config.num_of_linking_residues_max]
    if (config.hit_gc_content_percentage_min in config):
        level2 = level2[level2[titles.hit_gc_content] >=
                        config.hit_gc_content_percentage_min]
    if (config.hit_gc_content_percentage_max in config):
        level2 = level2[level2[titles.hit_gc_content] <=
                        config.hit_gc_content_percentage_max]
    if (config.precursor_mfei_min in config):
        level2 = level2[level2[titles.pre_mfei] >= config.precursor_mfei_min]
    if (config.precursor_mfei_max in config):
        level2 = level2[level2[titles.pre_mfei] <= config.precursor_mfei_max]

    if (config.max_allowed_mismatch_size_in_hit_region in config):
        level2 = level2[
            level2[["mismatch type", "mismatch size"]].apply(
                lambda row: _is_allowed(
                    row,
                    "mismatch type",
                    "mismatch size",
                    config.max_allowed_mismatch_size_in_hit_region,
                ),
                axis=1,
            )
        ]

    if (config.max_allowed_bulge_size_in_hit_region in config):
        level2 = level2[
            level2[["bulge type", "bulge size"]].apply(
                lambda row: _is_allowed(
                    row,
                    "bulge type",
                    "bulge size",
                    config.max_allowed_bulge_size_in_hit_region,
                ),
                axis=1,
            )
        ]
    if (config.max_allowed_internal_loop_size_in_hit_region in config):
        level2 = level2[
            level2[["internal type", "internal loop total size"]].apply(
                lambda row: _is_allowed(
                    row,
                    "internal type",
                    "internal loop total size",
                    config.max_allowed_internal_loop_size_in_hit_region,
                ),
                axis=1,
            )
        ]
    if (config.minimum_required_clear_region in config):
        level2 = level2[
            level2[[titles.distal_junc_dist, titles.prx_junc_dist]].apply(
                lambda row: _is_allowed_clear(
                    row, config.minimum_required_clear_region),
                axis=1,
            )
        ]
    if (config.border_line_mismatch_max in config):
        level2 = level2[
            level2[["mismatch type", "mismatch size"]].apply(
                lambda row: _check_border_line(
                    row, "mismatch type", "mismatch size", config.border_line_mismatch_max
                ),
                axis=1,
            )
        ]
    if (config.border_line_bulge_max in config):
        level2 = level2[
            level2[["bulge type", "bulge size"]].apply(
                lambda row: _check_border_line(
                    row, "bulge type", "bulge size", config.border_line_bulge_max
                ),
                axis=1,
            )
        ]
    if (config.border_line_internal_max in config):
        level2 = level2[
            level2[["internal type", "internal loop total size"]].apply(
                lambda row: _check_border_line(
                    row,
                    "internal type",
                    "internal loop total size",
                    config.border_line_internal_max,
                ),
                axis=1,
            )
        ]
    if (config.max_allowed_hsbl_ssbl_size in config):
        level2 = level2[
            level2[["internal type", "internal loop HSBL", "internal loop SSBL"]].apply(
                lambda row: _loopHSBL_SSBL_check(
                    row, config.max_allowed_hsbl_ssbl_size),
                axis=1,
            )
        ]
    if (config.acceptable_num_for_hit_locations_in_bulges_or_loops in config):
        level2 = level2[level2[titles.accept_num_hit_loc_bulg_loop] <=
                        config.acceptable_num_for_hit_locations_in_bulges_or_loops]
    if (config.acceptable_num_for_unmatched_locations_in_hit_region in config):
        level2 = level2[level2[titles.accept_num_unmatched_loc_hit_mayers] <=
                        config.acceptable_num_for_unmatched_locations_in_hit_region]
    if (config.total_num_of_nonmatching_positions in config):
        level2 = level2[level2[titles.total_num_nonmatching_pos]
                        <= config.total_num_of_nonmatching_positions]
    if (config.total_num_of_mismached_positions in config):
        level2 = level2[level2[titles.total_num_mismached_pos]
                        <= config.total_num_of_mismached_positions]
    if (config.total_num_of_positions_in_bulges_and_loops in config):
        level2 = level2[level2[titles.total_num_positions_bulges_loops] <=
                        config.total_num_of_positions_in_bulges_and_loops]
    if (config.delete_if_mature_duplex_involvement_in_apical_loop in config):
        if config.delete_if_mature_duplex_involvement_in_apical_loop == "YES":
            cols = [
                titles.n_term_struc,
                titles.hit_start,
                titles.hit_end,
                titles.star_start,
                titles.star_end,
                "branch#1 apical loop start",
                "branch#1 apical loop end",
            ]
            level2 = level2[level2[cols].apply(
                lambda row: _check_involvement(row), axis=1)]

    if (config.border_line_structure_allowance in config):
        border_proximal = level2[titles.sum_border_proximal] == 0
        border_distal = level2[titles.sum_border_distal] == 0
        if config.border_line_structure_allowance == "NOT ACCEPTED":
            level2 = level2[border_proximal & border_distal]
        if config.border_line_structure_allowance == "1 END ONLY":
            level2 = level2[border_proximal | border_distal]
    return level2


# BORDER_LINE_STRUCTURE_ALLOWANCE = "1 END ONLY"  # "NOT ACCEPTED"   "1 END ONLY" "BOTH END"
def postprocess(input_file, output_file, config=None, chunksize=10**5):
    if config is None:  # read from filter_level2.json
        with open(os.path.dirname(__file__) + "/config/filter_postprocess.json") as json_file:
            config = json.load(json_file)
    config = DotDict(config)

    header = True
    for chunk in tqdm.tqdm(pd.read_csv(input_file, chunksize=chunksize)):
        chunk = chunk.apply(lambda row: convert(row), axis=1)
        level2 = _postprocess_run(chunk, config)
        level2.to_csv(output_file, header=header, mode="a", index=False)
        header = False
