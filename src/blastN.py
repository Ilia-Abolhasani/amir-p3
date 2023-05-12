#!/usr/bin/env python
import os
import pandas as pd
from read_configs import *
from utils import *


def start(
    input_genome_path,
    temp_path,
    nonconformity,
    flanking_value,
    num_cpus
):
    # BlastN
    os.system(
        "makeblastdb -in {input_genome_path} -dbtype nucl -out {temp_path}/blastn_database"
    )

    header = "qseqid sseqid qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe sstrand qcovs qcovhsp qlen slen"

    os.system(
        "blastn -query ./{temp_path}/BLASTn_queries.fasta         -out ./{temp_path}/BLASTn_result         -num_threads {num_cpus}         -db ./{temp_path}/blastn_database         -word_size 7         -penalty -3         -reward 2         -gapopen 5         -gapextend 2         -outfmt '6 {header}'       "
    )

    df_blastn = pd.read_csv(
        f"{temp_path}/BLASTn_result", sep="\t", header=None)
    df_blastn.columns = header.replace("  ", " ").split(" ")

    # alignment length adjustment
    def blastn_adjust(row):
        if row["sstrand"] == "plus":
            row["sstart"] = max(1, row["sstart"] - (row["qstart"] - 1))
            row["send"] = min(row["slen"], row["send"] +
                              (row["qlen"] - row["qend"]))
        if row["sstrand"] == "minus":
            row["send"] = max(1, row["send"] - (row["qstart"] - 1))
            row["sstart"] = min(row["slen"], row["sstart"] +
                                (row["qlen"] - row["qend"]))
        return row

    df_blastn = df_blastn.apply(lambda row: blastn_adjust(row), axis=1)

    df_blastn["Nonconformity"] = (
        df_blastn["qlen"]
        - (abs(df_blastn["qend"] - df_blastn["qstart"]) + 1)
        + df_blastn["gaps"]
        + df_blastn["mismatch"]
    )
    df_blastn = df_blastn[df_blastn["Nonconformity"] <= nonconformity]

    # remore redundancy and hold best one base of Nonconformity value
    df_blastn = df_blastn.sort_values(
        ["Nonconformity", "evalue"], ascending=(True, True))
    df_blastn = df_blastn.drop_duplicates(
        subset=["sseqid", "sstart", "qseqid", "send", "sstrand"], keep="first"
    )
    df_blastn.to_csv(f"{temp_path}/filtered_out_blastn.csv")

    # Result of the blastn to bed file
    df = df_blastn[["qseqid", "sseqid", "sstart", "send", "sstrand", "slen"]]
    df["ones"] = 1

    def switch(row):
        if row["sstart"] > row["send"]:
            temp = row["sstart"]
            row["sstart"] = row["send"]
            row["send"] = temp
        return row

    df = df.apply(lambda row: switch(row), axis=1)

    def convert(inp):
        if inp == "plus":
            return "forward"
        if inp == "minus":
            return "reverse"
        raise Exception(
            'Error, sstrand contains illegal word! only "plus" and "minus" are allowed'
        )

    df["strand"] = df["sstrand"].apply(lambda x: convert(x))

    def convert2sign(inp):
        if inp == "plus":
            return "+"
        if inp == "minus":
            return "-"
        raise Exception(
            'Error, sstrand contains illegal word! only "plus" and "minus" are allowed'
        )

    df["sign"] = df["sstrand"].apply(lambda x: convert2sign(x))

    df["hit_length"] = df.apply(lambda row: abs(
        row["send"] - row["sstart"]) + 1, axis=1)

    # convert sstart and send from location to index (range)
    df["sstart"] = df["sstart"].apply(lambda x: x - 1)

    df["downstream_flanking"] = df["sstart"].apply(
        lambda x: flanking_value if x > flanking_value else x
    )
    df["upstream_flanking"] = df.apply(
        lambda row: flanking_value
        if (row["send"] + flanking_value) <= row["slen"]
        else row["slen"] - row["send"],
        axis=1,
    )

    df["hit_start"] = df.apply(
        lambda row: row["downstream_flanking"]
        if row["sign"] == "+"
        else row["upstream_flanking"],
        axis=1,
    )
    df["hit_end"] = df.apply(
        lambda row: row["downstream_flanking"] + row["hit_length"]
        if row["sign"] == "+"
        else row["upstream_flanking"] + row["hit_length"],
        axis=1,
    )
    df["sstart"] = df["sstart"].apply(lambda x: max(x - flanking_value, 0))
    df["send"] = df.apply(
        lambda row: min(row["send"] + flanking_value, row["slen"]), axis=1
    )
    df["tag"] = df.apply(
        lambda row: f">{row['sseqid']}:{row['sstart']}-{row['send']}({row['sign']})", axis=1
    )
    df["reformated_tag"] = df["tag"].apply(lambda t: reformat(t))
    df[["tag", "reformated_tag", "hit_start", "hit_end"]].to_csv(
        f"./{temp_path}/hit_index_info.csv"
    )  # , index=False)
    df["location_tag"] = df.apply(
        lambda row: f">{row['sseqid']}|{row['sign']}|{row['sstart'] + 1}-{row['send']}|{row['hit_start']+1}-{row['hit_end']}",
        axis=1,
    )
    df[["location_tag", "qseqid"]].to_csv(
        f"{temp_path}/pipe_seprated_location_list.csv", index=False, sep="\t"
    )
    df[["sseqid", "sstart", "send", "strand", "ones", "sign"]].to_csv(
        f"{temp_path}/extension_index.bed", index=False, header=False, sep="\t"
    )
