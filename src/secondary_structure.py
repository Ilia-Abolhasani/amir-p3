# RNA 2d prediction
import os
import glob
from tqdm import tqdm
import shutil
import pandas as pd
import multiprocessing as mp
from utils import reformat, fasta_to_df, bracket_row, bracket_to_ct

## Mfold
def _mfold(result_path, num_cpus):
    base = f"{result_path}/secondary_structure/mfold/"
    os.system("rm -r {base}")
    os.system("mkdir -p {base}")
    df = fasta_to_df(f"{temp_path}/extended_modified_non_coding.txt")

    for index, row in df.iterrows():
        tag = reformat(row["tag"])
        if not os.path.exists(base + tag):
            os.makedirs(base + tag)
        with open(base + f"{tag}/SEQ.FASTA", "w") as file:
            file.write(f">{row['tag']}\n{row['data']}")

    os.run_cell_magic(
        "capture",
        "",
        'remove_lock = False\ndef run_mfold(tag):\n    tag = reformat(tag)\n    %cd {base + tag}\n    !mfold  SEQ="SEQ.FASTA" T=22   \n    if(not remove_lock):\n        !find . -name "SEQ*" -not -name "*.ct" -not -name "*.pdf" -not -name "*SEQ.FASTA" -not -type d -delete\n    %cd {current_path}\n\nif __name__ == \'__main__\':        \n    pool = mp.Pool(mp.cpu_count() - 3)      \n    pool.map(run_mfold, df[\'tag\'])  ',
    )


## Mxfold2
def _mxfold2(result_path, num_cpus):
    os.system(
        "mxfold2 predict ./extended.txt > Result/secondary_structure/mxfold2_result.txt"
    )

    df = fasta_to_df("./Result/secondary_structure/mxfold2_result.txt")
    df = df.apply(lambda row: bracket_row(row), axis=1)
    df.head(2)
    base = "./Result/secondary_structure/mxfold2/"
    os.system("rm -r {base}")
    os.system("mkdir -p {base}")
    for index, row in df.iterrows():
        if not os.path.exists(base + reformat(row["tag"])):
            os.makedirs(base + reformat(row["tag"]))
        tag = reformat(row["tag"])
        with open(base + f"{tag}/{tag}.ct", "w") as file:
            bracket = row["bracket"].split(" ")[0]
            deltaG = row["bracket"].split(" ")[1]
            ct = bracket_to_ct(row["tag"], row["data"], bracket, deltaG)
            file.write(ct)


## Vienna package
def _viennarna(result_path, num_cpus):
    base = f"{result_path}/secondary_structure/viennarna/"
    os.system("rm -r {base}")
    os.system("rm {result_path}/secondary_structure/viennarna_result.txt")
    os.system("mkdir -p {base}")

    os.system(
        "RNAfold --jobs=0 --infile {current_path}/{temp_path}/extended_modified.txt  --noPS -T 22 > {current_path}/{base}/viennarna_result.txt"
    )

    df = fasta_to_df(
        f"{result_path}/secondary_structure/viennarna/viennarna_result.txt"
    )
    df = df.apply(lambda row: bracket_row(row), axis=1)
    print(df.shape)
    df.head(2)

    for index, row in df.iterrows():
        tag = reformat(row["tag"])
        if not os.path.exists(base + tag):
            os.makedirs(base + tag)
        with open(base + f"{tag}/{tag}.ct", "w") as file:
            bracket = row["bracket"].split(" ")[0]
            deltaG = row["bracket"].split(" ")[1]
            ct = bracket_to_ct(row["tag"], row["data"], bracket, deltaG, False)
            file.write(ct)

    for file in glob.glob(f"{base}*.ps"):
        f = file[len(base) : -6]  # _ss.ps
        f = reformat(f)
        shutil.move(file, f"{base}{f}/{f}.ps")


## ContraFold
def _contrafold(result_path, num_cpus):
    counter = 0
    base = f"./{result_path}/secondary_structure/contrafold/"
    os.system("rm -r {base}")
    os.system("mkdir -p {base}")
    df = fasta_to_df(f"{temp_path}/extended.txt")

    for index, row in tqdm(df.iterrows()):
        tag = reformat(row["tag"])
        if not os.path.exists(base + tag):
            os.makedirs(base + tag)
        with open(base + f"{tag}/{tag}.FASTA", "w") as file:
            file.write(f">{row['tag']}\n{row['data']}")
        counter += 1

    def run_contrafold(tag):
        tag = reformat(tag)
        os.run_line_magic("cd", "Software/contrafold/src")
        os.system(
            "./contrafold predict ../..{base[1:]}{tag}/{tag}.FASTA > ../..{base[1:]}{tag}/{tag}.dot"
        )
        with open(f"../..{base[1:]}{tag}/{tag}.dot", "r") as file:
            text = file.read()
        text = [l for l in text.split("\n") if l[: len(">structure")] != ">structure"]
        header = text[0]
        with open(f"../..{base[1:]}{tag}/{tag}.dot", "w") as file:
            file.write("\n".join(text[1:]))
        os.system(
            "RNAeval  ../..{base[1:]}{tag}/{tag}.dot -T 20 > ../..{base[1:]}{tag}/{tag}.dotdg    "
        )
        with open(f"../..{base[1:]}{tag}/{tag}.dotdg", "r") as file:
            text = file.read()
        with open(f"../..{base[1:]}{tag}/{tag}.dot", "w") as file:
            file.write(header + "\n" + text)

        df = fasta_to_df(f"../..{base[1:]}{tag}/{tag}.dot")
        df = df.apply(lambda row: bracket_row(row), axis=1)
        tag = reformat(df["tag"][0])
        with open(f"../..{base[1:]}{tag}/{tag}.ct", "w") as file:
            bracket = df["bracket"][0].split(" ")[0]
            deltaG = df["bracket"][0].split(" ")[1]
            ct = bracket_to_ct(df["tag"][0], df["data"][0], bracket, deltaG, False)
            file.write(ct)
        #!rm ../..{base[1:]}{tag}/{tag}.dot
        #!rm ../..{base[1:]}{tag}/{tag}.dotdg
        os.system("rm ../..{base[1:]}{tag}/{tag}.FASTA")
        os.run_line_magic("cd", "{current_path}")

    if __name__ == "__main__":
        pool = mp.Pool(num_cpus)
        pool.map(run_contrafold, df["tag"].iloc[:10])


def secondary_structure(method, result_path, num_cpus):
    method = method.lower()
    if method == "mfold":
        _mfold(result_path, num_cpus)
    elif method == "viennarna":
        _viennarna(result_path, num_cpus)
    elif method == "contrafold":
        _contrafold(result_path, num_cpus)
    elif method == "mxfold2":
        _mxfold2(result_path, num_cpus)
    else:
        pass
