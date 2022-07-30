import numpy as np
import pandas as pd

if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO


def adjust(text, n=7):
    text = str(text)
    return " " * (n - len(text)) + text


def bracket_to_ct(tag, data, bracket, deltaG, negative_deltaG=True):
    deltaG = deltaG.replace("(", "").replace(")", "")
    deltaG = float(deltaG)
    if deltaG > 0 and negative_deltaG:  # negetive?!
        deltaG = -1 * deltaG
    stack = []
    index = np.zeros((len(bracket)), dtype=int)
    values = np.zeros((len(bracket)), dtype=int)
    for i in range(len(bracket)):
        index[i] = i + 1
        if bracket[i] == ".":
            values[i] = 0
        elif bracket[i] == "(":
            stack.append(i)
        elif bracket[i] == ")":
            if len(stack) == 0:
                print("structure error!")
            values[stack[-1]] = i + 1
            values[i] = stack[-1] + 1
            stack.pop()
        else:
            print("structure error!")
    if len(stack) != 0:
        print("structure error!")
    # body
    ct = f"{adjust(len(data),6)} dG ={adjust(deltaG,10)} {tag}\n"
    for i in range(len(bracket)):
        ct += f"{adjust(index[i],6)} {data[i]} {adjust(i,6)} {adjust((i+2)%(len(data)+1),6)} {adjust(values[i],6)} {adjust(index[i],7)}\n"
    return ct


def fasta_to_df(path):
    with open(path, "r") as file:
        text = file.read()
    lines = [line for line in text.split("\n") if len(line) > 0]
    s = ""
    tags = []
    data = []
    for l in lines:
        if l[0] == ">":
            tags.append(l)
            data.append(s)
            s = ""
        else:
            s += l
    data.append(s)
    df = pd.DataFrame({"tag": tags, "data": data[1:]})
    df["tag"] = df["tag"].apply(lambda x: x[1:])
    return df


def df_to_fasta(df, path):
    lines = []
    df.apply(lambda row: lines.append(f">{row['tag']}\n{row['data']}\n"), axis=1)
    with open(path, "w") as file:
        file.write("".join(lines))


def reformat(path):
    return path.replace("(", "_").replace(")", "_").replace(".", "").replace(":", "_")


def reformatCT(path):
    with open(path, "r") as file:
        text = file.read()
    text = [l for l in text.split("\n") if len(l) > 0]  # remove blank lines
    text = "\n".join(text)
    text = text.replace("\t", " ")
    while "  " in text:
        text = text.replace("  ", " ")
    lines = [l for l in text.split("\n")]
    for i in range(len(lines)):
        if lines[i][0] == " ":
            lines[i] = lines[i][1:]
        if lines[i][-1] == " ":
            lines[i] = lines[i][:-1]
    text = "\n".join(lines)
    return text


def get_ct_data(ct):
    ct = "\n".join(ct.split("\n")[1:])
    df = pd.read_csv(StringIO(ct), sep=" ", header=None)
    nucleotide = df.iloc[:, 1]
    index = df.iloc[:, 5]
    values = df.iloc[:, 4]
    return [nucleotide, index, values]


def ct2dot_bracket(path):
    [nucleotide, index, values] = get_ct_data(reformatCT(path))
    text = "".join(nucleotide) + "\n"
    watch = []
    for i, v in zip(index, values):
        if v == 0:
            text += "."
        else:
            if v not in watch:
                text += "("
                watch.append(i)
            if v in watch:
                text += ")"
    return text


def is_nested(index, values):
    max_value = max(index) + 10  # inf
    for i, v in zip(index, values):
        if v < max_value and v != 0:
            max_value = v
        if i >= max_value:
            max_value = max(index) + 10  # inf
        if v > max_value:
            return False
    return True
