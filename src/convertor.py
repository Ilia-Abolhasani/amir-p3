from read_configs import read_titles, read_erros


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
        titles.mismatch_type,
        titles.mismatch_size,
        titles.mismatch_start,
        titles.mismatch_end,
        titles.bulge_type,
        titles.bulge_size,
        titles.bulge_start,
        titles.bulge_end,
        titles.bulge_strand,
        titles.inter_type,
        titles.inter_size,
        titles.inter_start,
        titles.inter_end,
        titles.inter_HSBL,
        titles.inter_SSBL
    ]:
        df[col] = df[col].apply(lambda x: eval(x))
    return df
