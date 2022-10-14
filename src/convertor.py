from read_configs import DotDict, read_titles, read_erros


def convert(df):
    titles = read_titles()
    df[titles.pre_mfei] = df[titles.pre_mfei].apply(
        lambda x: 10 * 10 if x == "-" else x
    )
    for col in [
        titles.dg,
        titles.n_term_struc,
        titles.num_of_lnk_res,
        titles.boi_gc,
        titles.pre_mfei,
    ]:
        df[col] = df[col].apply(lambda x: float(x))
    for col in [
        "mismatch type",
        "mismatch size",
        "mismatch start",
        "mismatch end",
        "bulge type",
        "bulge size",
        "bulge start",
        "bulge end",
        "bulge strand",
        "internal type",
        "internal loop total size",
        "internal start",
        "internal end",
        "internal loop HSBL",
        "internal loop SSBL",
    ]:
        df[col] = df[col].apply(lambda x: eval(x))
    return df
