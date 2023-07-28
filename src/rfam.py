#!/usr/bin/env python
import os
import pandas as pd
from read_configs import *
from utils import *


def start(
    extended_path,
    temp_path,
    evalue,
    num_cpus,
):
    # concat all .fa in rfam
    rfam_file = "./data/rfam/rfam.fa"
    if os.path.isfile(rfam_file):
        os.remove(rfam_file)
    os.system("cat ./data/rfam/*.fa > " + rfam_file)

    # create DB
    os.system(
        f"makeblastdb -in {rfam_file} -dbtype nucl -out {temp_path}/rfam_blastn_database"
    )

    header = "qseqid sseqid qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe sstrand qcovs qcovhsp qlen slen"

    os.system(
        f"blastn -query {extended_path} -out {temp_path}/rfam_blastn_result -num_threads {num_cpus} -db {temp_path}/rfam_blastn_database -word_size 28 -penalty -2 -reward 1 -gapopen 5 -gapextend 2 -outfmt '6 {header}' "
    )

    df_blastn = pd.read_csv(
        f"{temp_path}/rfam_blastn_result", sep="\t", header=None)
    df_blastn.columns = header.replace("  ", " ").split(" ")
    df_blastn = df_blastn[df_blastn['evalue'] < evalue]

    #
    ext = fasta_to_df(f"{extended_path}")
    print(f"total:      {ext.shape[0]}")
    filtered = ext[~ext["tag"].isin(df_blastn['qseqid'])]
    print(f"after remove rfam: {filtered.shape[0]}")
    df_to_fasta(filtered, f"{temp_path}/extended_modified_rfam.txt")
