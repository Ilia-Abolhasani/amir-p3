# RNA 2d prediction
import os
import glob
from tqdm import tqdm
import shutil
import pandas as pd
import multiprocessing as mp
from utils import *
import functools


def _bracket_row(row):
    s = row['data']
    index = min(s.find('.'), s.find('('))
    data = row['data']
    row['data'] = data[0:index]
    row['bracket'] = data[index:]
    return row

# Mfold


def _run_mfold(tag, current_path, base, folding_temperature, remove_extra_files):
    tag = reformat(tag)
    os.chdir(f"{base + tag}")
    os.system(f'''mfold  SEQ="SEQ.FASTA" T={folding_temperature}''')
    if (remove_extra_files):
        os.system(
            '''find . -name "SEQ*" -not -name "*.ct" -not -name "*.pdf" -not -name "*SEQ.FASTA" -not -type d -delete'''
        )
    os.system("rm SEQ.ct")
    os.chdir(f"{current_path}")


def _mfold(extended_path, result_path, folding_temperature,  num_cpus, remove_extra_files=False):
    current_path = os.getcwd()
    base = f"{result_path}/secondary_structure/mfold/"
    os.system(f"rm -r {base}")
    os.system(f"mkdir -p {base}")
    df = fasta_to_df(extended_path)

    for index, row in df.iterrows():
        tag = reformat(row["tag"])
        if not os.path.exists(base + tag):
            os.makedirs(base + tag)
        with open(base + f"{tag}/SEQ.FASTA", "w") as file:
            file.write(f">{row['tag']}\n{row['data']}")
            print(len(row['data']))
    input("wait....")
    if __name__ == 'secondary_structure':
        pool = mp.Pool(num_cpus)
        pool.map(functools.partial(_run_mfold, current_path=current_path, base=base,
                 folding_temperature=folding_temperature, remove_extra_files=remove_extra_files), df['tag'])


# Mxfold2
def _mxfold2(extended_path, result_path, folding_temperature,  num_cpus):
    os.system(
        f"mxfold2 predict {extended_path} > {result_path}/secondary_structure/mxfold2_result.txt"
    )

    df = fasta_to_df(f"{result_path}/secondary_structure/mxfold2_result.txt")
    df = df.apply(lambda row: _bracket_row(row), axis=1)
    base = f"{result_path}/secondary_structure/mxfold2/"
    os.system(f"rm -r {base}")
    os.system(f"mkdir -p {base}")
    for index, row in df.iterrows():
        if not os.path.exists(base + reformat(row["tag"])):
            os.makedirs(base + reformat(row["tag"]))
        tag = reformat(row["tag"])
        with open(base + f"{tag}/{tag}.ct", "w") as file:
            bracket = row["bracket"].split(" ")[0]
            deltaG = row["bracket"].split(" ")[1]
            ct = bracket_to_ct(row["tag"], row["data"], bracket, deltaG)
            file.write(ct)


# Vienna package
def _viennarna(extended_path, result_path, folding_temperature,  num_cpus):
    current_path = os.getcwd()
    base = f"{result_path}/secondary_structure/viennarna/"
    os.system(f"rm -r {base}")
    os.system(f"rm {result_path}/secondary_structure/viennarna_result.txt")
    os.system(f"mkdir -p {base}")

    os.system(
        f"RNAfold -v --jobs={num_cpus} --infile {current_path}/{extended_path}  --noPS -T {folding_temperature} > {current_path}/{base}/viennarna_result.txt"
    )

    df = fasta_to_df(
        f"{result_path}/secondary_structure/viennarna/viennarna_result.txt"
    )
    df = df.apply(lambda row: _bracket_row(row), axis=1)

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
        f = file[len(base): -6]  # _ss.ps
        f = reformat(f)
        shutil.move(file, f"{base}{f}/{f}.ps")


# ContraFold
def _run_contrafold(tag, base, current_path, folding_temperature):
    tag = reformat(tag)
    os.system(
        f'''./software/contrafold/src/contrafold predict "{base}{tag}/{tag}.FASTA" > "{base}{tag}/{tag}.dot" '''
    )
    with open(f"{base}{tag}/{tag}.dot", "r") as file:
        text = file.read()
    text = [l for l in text.split(
        "\n") if l[:len(">structure")] != ">structure"]
    header = text[0]
    with open(f"{base}{tag}/{tag}.dot", "w") as file:
        file.write("\n".join(text[1:]))
    os.system(
        f'''RNAeval "{base}{tag}/{tag}.dot" -T {folding_temperature} > "{base}{tag}/{tag}.dotdg" '''
    )
    with open(f"{base}{tag}/{tag}.dotdg", "r") as file:
        text = file.read()
    with open(f"{base}{tag}/{tag}.dot", "w") as file:
        file.write(header + "\n" + text)

    df = fasta_to_df(f'{base}{tag}/{tag}.dot')
    df = df.apply(lambda row: _bracket_row(row), axis=1)
    tag = reformat(df["tag"][0])
    with open(f"{base}{tag}/{tag}.ct", "w") as file:
        bracket = df["bracket"][0].split(" ")[0]
        deltaG = df["bracket"][0].split(" ")[1]
        deltaG = "".join(deltaG).replace(
            " ", "").replace("(", "").replace(")", "")
        ct = bracket_to_ct(df["tag"][0], df["data"]
                           [0], bracket, deltaG, False)
        file.write(ct)


def _contrafold(extended_path, result_path, folding_temperature,  num_cpus):
    current_path = os.getcwd()
    counter = 0
    base = f"{result_path}/secondary_structure/contrafold/"
    os.system(f"rm -r {base}")
    os.system(f"mkdir -p {base}")
    df = fasta_to_df(extended_path)

    for index, row in tqdm(df.iterrows()):
        tag = reformat(row["tag"])
        if not os.path.exists(base + tag):
            os.makedirs(base + tag)
        with open(base + f"{tag}/{tag}.FASTA", "w") as file:
            file.write(f">{row['tag']}\n{row['data']}")
        counter += 1

    if __name__ == 'secondary_structure':
        pool = mp.Pool(num_cpus)
        pool.map(functools.partial(_run_contrafold, base=base,
                 current_path=current_path, folding_temperature=folding_temperature), df["tag"])


def start(method, extended_path, result_path, folding_temperature,  num_cpus):
    print("Starting secondary structure prediction...")
    print("Please wait, as this process may take a while...")
    method = method.lower()
    if method == "mfold":
        _mfold(extended_path, result_path, folding_temperature,  num_cpus)
    elif method == "viennarna":
        _viennarna(extended_path, result_path, folding_temperature,  num_cpus)
    elif method == "contrafold":
        _contrafold(extended_path, result_path, folding_temperature,  num_cpus)
    elif method == "mxfold2":
        _mxfold2(extended_path, result_path, folding_temperature,  num_cpus)
    else:
        raise Exception('Error, unrecognized secondary structure method.')
    print("Secondary structure prediction completed.")
