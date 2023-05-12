import os
from utils import *
import pandas as pd


def start(temp_path):
    # Protein coding elimination [Diamond]
    os.system("./diamond makedb --in ./NR/nr -d ./Temp/diamond_output")

    os.system(
        "./diamond blastx -d ./Temp/diamond_output.dmnd                   -q ./Temp/extended_modified.txt                   -o ./Temp/diamond_matches.tsv                   -p 22"
    )

    dmn = pd.read_csv(f"{temp_path}/diamond_matches.tsv",
                      sep="\t", header=None)
    dmn.columns = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split(
        " "
    )
    coding_seq = dmn["qseqid"].unique()

    def clear(inp):
        if inp[:9] == "reverse::":
            return inp[9:]
        if inp[:9] == "forward::":
            return inp[9:]
        return inp

    coding_seq = pd.Series(coding_seq).apply(lambda x: clear(x))

    ext = fasta_to_df(f"{temp_path}/extended_modified.txt")
    print(f"total:      {ext.shape[0]}")
    non_coding = ext[~ext["tag"].isin(coding_seq)]
    print(f"non_coding: {non_coding.shape[0]}")
    df_to_fasta(non_coding, f"{temp_path}/extended_modified_non_coding.txt")
    coding = ext[ext["tag"].isin(coding_seq)]
    print(f"coding:     {coding.shape[0]}")
    df_to_fasta(coding, f"{temp_path}/extended_modified_coding.txt")
