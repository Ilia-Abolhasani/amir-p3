import json

class DotDict(dict):
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


with open('./src/config/titles.json') as json_file:
    titles = DotDict(json.load(json_file))


def convert(df):
    df[titles.pre_mfei] = df[titles.pre_mfei].apply(lambda x: 10 * 10 if x == "-" else x)
    for col in [titles.dg, titles.n_term_struc, titles.num_of_lnk_res, titles.boi_gc, titles.pre_mfei]:
        df[col] = df[col].apply(lambda x: float(x))
    for col in ['mismatch type', 'mismatch size', 'mismatch start', 'mismatch end',
                 'bulge type', 'bulge size', 'bulge start', 'bulge end', "bulge strand",
                 'internal type', 'internal loop total size', 'internal start', 'internal end', 'internal loop HSBL',
                 'internal loop SSBL']:
        df[col] = df[col].apply(lambda x: eval(x))
    return df
