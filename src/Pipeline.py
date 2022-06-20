#!/usr/bin/env python
# coding: utf-8

# # Config

# In[1]:


experiment = "A.indica"
input_genome_name = "GCA_022749755.1_ASM2274975v1_genomic.fna"


experiment_dir = "Experiment"


# In[2]:


input_genome_path = f'{experiment_dir}/{experiment}/{input_genome_name}'

temp_path = f"{experiment_dir}/{experiment}/Temp"
result_path = f"{experiment_dir}/{experiment}/Result"

temp_path_f = temp_path.replace(" ", "\ ")
result_path_f = result_path.replace(" ", "\ ")


# # Common

# In[3]:


#!pip install tqdm


# In[4]:


import json
import time
from subprocess import Popen, PIPE, STDOUT
import math
import numpy as np
import pandas as pd
import hashlib
import requests
import os, sys, subprocess
from tqdm.contrib.concurrent import process_map
from tqdm.notebook import tqdm
tqdm.pandas()
import multiprocessing as mp
import shutil
import urllib.parse
import glob
import os
import sys
import networkx
from networkx.algorithms.clique import find_cliques as maximal_cliques
from ast import literal_eval
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO
sys.path.append("./src/")
from ct_analizer import get_row
from filter1 import filter1_run
from filter2 import filter2


# In[5]:


if(not os.path.exists(temp_path)):
    os.mkdir(temp_path)
    
if(not os.path.exists(result_path)):
    os.mkdir(result_path)


# In[6]:


current_path = os.getcwd()


# In[7]:


def bracket_row(row):    
    s = row['data']
    index = min(s.find('.'), s.find('('))
    data = row['data']
    row['data'] = data[0:index]
    row['bracket'] = data[index:]
    return row


# In[8]:


def adjust(text,n=7):
    text = str(text)    
    return " " * (n - len(text)) + text


# In[9]:


def bracket_to_ct(tag, data, bracket, deltaG, negative_deltaG=True):    
    deltaG = deltaG.replace('(','').replace(')','')
    deltaG = float(deltaG)
    if(deltaG > 0 and negative_deltaG ): # negetive?!
        deltaG = -1 * deltaG
    stack = []
    index = np.zeros((len(bracket)), dtype = int)
    values = np.zeros((len(bracket)), dtype = int)
    for i in range(len(bracket)):
        index[i] = i + 1
        if(bracket[i] == '.'):
            values[i] = 0
        elif(bracket[i] == '('):
            stack.append(i)
        elif(bracket[i] == ')'):
            if(len(stack) == 0 ):
                print('structure error!')
            values[stack[-1]] = i + 1
            values[i]  = stack[-1] + 1
            stack.pop()
        else:
            print('structure error!')
    if(len(stack) != 0 ):
        print('structure error!')
    # body    
    ct = f"{adjust(len(data),6)} dG ={adjust(deltaG,10)} {tag}\n"   
    for i in range(len(bracket)):
        ct += f"{adjust(index[i],6)} {data[i]} {adjust(i,6)} {adjust((i+2)%(len(data)+1),6)} {adjust(values[i],6)} {adjust(index[i],7)}\n"
    return ct


# In[10]:


def fasta_to_df(path):
    with open(path, 'r') as file:
        text = file.read()
    lines = [line for line in text.split('\n') if len(line) > 0]
    s = ''
    tags = []
    data = []
    for l in lines:
        if(l[0]=='>'):
            tags.append(l)        
            data.append(s)
            s = ''
        else:
            s += l    
    data.append(s)
    df = pd.DataFrame(
            {
                'tag': tags,
                'data': data[1:]
            })
    df['tag'] = df['tag'].apply(lambda x: x[1:])    
    return df


# In[11]:


def df_to_fasta(df, path):
    lines = []
    df.apply(lambda row: lines.append(f">{row['tag']}\n{row['data']}\n"),axis=1)
    with open(path,'w') as file:
        file.write(''.join(lines))


# In[12]:


def reformat(path):
    return path.replace('(','_').replace(')','_').replace('.','').replace(':','_')


# In[13]:


def reformatCT(path):
    with open(path, 'r') as file:
        text = file.read()
    text = [l for l in text.split('\n') if len(l) > 0 ] # remove blank lines
    text = '\n'.join(text)
    text = text.replace("\t"," ")
    while("  " in text):
        text = text.replace("  ", " ")
    lines = [l for l in text.split('\n')]
    for i in range(len(lines)):
        if(lines[i][0] == " "):
            lines[i] = lines[i][1:]
        if(lines[i][-1] == " "):
            lines[i] = lines[i][:-1]
    text = '\n'.join(lines)
    return text


# In[14]:


def get_ct_data(ct):
    ct = "\n".join(ct.split('\n')[1:])
    df = pd.read_csv(StringIO(ct), sep=" ", header=None)               
    nucleotide = df.iloc[:,1]
    index = df.iloc[:,5]
    values = df.iloc[:,4]
    return [nucleotide, index, values]


# In[15]:


def ct2dot_bracket(path):
    [nucleotide, index, values] = get_ct_data(reformatCT(path))
    text = ''.join(nucleotide) + "\n"
    watch = []
    for i, v in zip(index,values):
        if(v == 0):
            text += '.'
        else:
            if( v not in watch):
                text += '('
                watch.append(i)
            if( v in watch):
                text += ')'
    return text


# In[16]:


def is_nested(index, values):
    max_value = max(index) + 10 # inf
    for i, v in zip(index, values):
        if(v < max_value and v != 0):
            max_value  = v
        if(i >= max_value):
            max_value = max(index) + 10 # inf
        if(v > max_value):
            return False               
    return True


# # Download data from Mirbase

# In[17]:


directory = './miRBase_driven_data'
base = "https://www.mirbase.org/ftp/CURRENT"        


# In[18]:


get_ipython().system('rm -r {directory}')
get_ipython().system('mkdir -p {directory}')

get_ipython().system('wget {base}/aliases.txt.gz -P ./{directory}/       ; gzip -d ./{directory}/aliases.txt.gz ')
get_ipython().system('wget {base}/hairpin.fa.gz -P ./{directory}/           ; gzip -d ./{directory}/hairpin.fa.gz ')
get_ipython().system('wget {base}/hairpin_high_conf.fa.gz -P ./{directory}/ ; gzip -d ./{directory}/hairpin_high_conf.fa.gz ')
get_ipython().system('wget {base}/mature.fa.gz -P ./{directory}/            ; gzip -d ./{directory}/mature.fa.gz ')
get_ipython().system('wget {base}/mature_high_conf.fa.gz -P ./{directory}/  ; gzip -d ./{directory}/mature_high_conf.fa.gz')
get_ipython().system('wget {base}/miRNA.str.gz -P ./{directory}/            ; gzip -d ./{directory}/miRNA.str.gz ')
get_ipython().system('wget {base}/miRNA.xls.gz -P ./{directory}/            ; gzip -d ./{directory}/miRNA.xls.gz ')
get_ipython().system('wget {base}/organisms.txt.gz -P ./{directory}/        ; gzip -d ./{directory}/organisms.txt.gz')


# In[18]:


mature = fasta_to_df(f'{directory}/mature.fa')
mature_high_conf = fasta_to_df(f'{directory}/mature_high_conf.fa')
mature['trim tag'] = mature['tag'].apply(lambda line: ' '.join(line.split(' ')[:2]))
mature['confidence'] = mature['trim tag'].isin(mature_high_conf['tag'])


# In[19]:


mature['organism'] = mature['tag'].apply(lambda x: x[:3])
print(mature.shape)
mature.head(2)


# In[20]:


organism = pd.read_csv(f'./{directory}/organisms.txt',sep='\t')
organism.columns = [c.replace('#','') for c in organism.columns] # remove sharp from columns
print(organism.shape)
organism.head(2)


# In[21]:


items = list(organism['tree'].unique())
items.sort(key=len)
items


# In[22]:


selectedTree = organism[organism['tree'].apply(lambda x: "Viridiplantae;" in x)]
print(selectedTree.shape)
selectedTree.head(5)


# In[23]:


#selectedTree = selectedTree[selectedTree['name'] == ""]


# In[24]:


selected = mature[mature['organism'].isin(selectedTree['organism'])]
print(selected.shape)
selected.head(1)


# In[25]:


df_to_fasta(selected,f'{temp_path}/mature_microRNA_queries.fasta')


# # Remove redundant

# ## cdhit-est

# In[26]:


get_ipython().system('./Software/cdhit/cd-hit-est -i ./{temp_path_f}/mature_microRNA_queries.fasta  -o ./{temp_path_f}/NR_mature_microRNA_queries.fasta     -c 1 -r 0 -G 1 -g 1 -b 30 -l 10 -aL 0 -AL 99999999 -aS 0     -AS 99999999 -s 0 -S 0')


# ## reformat

# In[27]:


with open(f'./{temp_path}/NR_mature_microRNA_queries.fasta.clstr','r') as file:
    text = file.read()
lines = [line for line in text.split('\n') if len(line) > 0]
cluster = []
seqid = []
last_cluster = ""
for l in lines:
    if(l[0]=='>'):        
        last_cluster = l.replace('>Cluster ',"C")
    else:        
        cluster.append(last_cluster)
        seqid.append(l.split(', >')[1].split('...')[0])                
seq2cluster = pd.DataFrame({'seqid': seqid,'cluster': cluster})
print(seq2cluster.shape)
seq2cluster.head(2)    


# In[28]:


df = fasta_to_df(f"./{temp_path}/mature_microRNA_queries.fasta")
df['accession'] = df['tag'].apply(lambda x : x.split(' ')[0])
seq2cluster = pd.merge(df,seq2cluster,how="inner",left_on='accession',right_on="seqid")
seq2cluster = pd.merge(seq2cluster, mature,how="inner",left_on='tag',right_on="tag")[['cluster','seqid','tag', 'confidence']]
print(seq2cluster.shape)
display(seq2cluster.head(2))
seq2cluster.to_csv(f'./{temp_path}/seq2cluster.csv',index=False)


# In[29]:


# todo: sorted first by cluster then by seqid
seq2cluster.sort_values("cluster").head(2)


# In[30]:


df = fasta_to_df(f"./{temp_path}/NR_mature_microRNA_queries.fasta")
df['tag'] = df['tag'].apply(lambda x : x.split(' ')[0])
df = pd.merge(df,seq2cluster,how="inner",left_on='tag',right_on="seqid")[['cluster','data']]

lines = []
df.apply(lambda row: lines.append(f">{row['cluster']}\n{row['data']}\n"),axis=1)
print(df.shape)
with open(f'./{temp_path}/BLASTn_queries.fasta','w') as file:
    file.write(''.join(lines))


# # BlastN

# !sudo apt-get install ncbi-blast+
# 

# In[31]:


path = f'./{experiment_dir}/{experiment}/{input_genome_name}'
get_ipython().system('makeblastdb -in {path} -dbtype nucl -out ./{temp_path_f}/blastn_database')


# In[32]:


header = 'qseqid sseqid qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe sstrand qcovs qcovhsp qlen slen'


# In[33]:


get_ipython().system("blastn -query ./{temp_path}/BLASTn_queries.fasta         -out ./{temp_path}/BLASTn_result         -num_threads {mp.cpu_count()}         -db ./{temp_path}/blastn_database         -word_size 7         -penalty -3         -reward 2         -gapopen 5         -gapextend 2         -outfmt '6 qseqid sseqid qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe sstrand qcovs qcovhsp qlen slen'       ")


# In[34]:


df_blastn = pd.read_csv(f'./{temp_path}/BLASTn_result', sep='\t',header=None)
df_blastn.columns = header.replace("  "," ").split(" ")
print(df_blastn.shape)
df_blastn.head(2)


# In[35]:


# alignment length adjustment
def blastn_adjust(row):
    if(row['sstrand'] == "plus"):        
        row['sstart'] = max(1, row['sstart'] - (row['qstart'] - 1))
        row['send'] = min(row['slen'], row['send'] + (row['qlen'] - row['qend']))
    if(row['sstrand'] == "minus"):        
        row['send'] = max(1, row['send'] - (row['qstart'] - 1 ))
        row['sstart'] = min(row['slen'], row['sstart'] + (row['qlen'] - row['qend'] ))
    return row
    
df_blastn = df_blastn.apply(lambda row: blastn_adjust(row), axis=1)


# In[36]:


threshold = 3
df_blastn['Nonconformity'] = df_blastn['qlen'] - (abs(df_blastn['qend'] - df_blastn['qstart']) + 1) + df_blastn['gaps'] + df_blastn['mismatch']
df_blastn = df_blastn[df_blastn['Nonconformity'] <= threshold]
print(df_blastn.shape)
df_blastn.head(2)


# In[37]:


# remore redundancy and hold best one base of Nonconformity value
df_blastn = df_blastn.sort_values(["Nonconformity", "evalue"], ascending = (True, True))
df_blastn = df_blastn.drop_duplicates(subset=['sseqid','sstart', 'qseqid', 'send','sstrand'], keep='first')
df_blastn.to_csv(f'./{temp_path}/filtered_out_blastn.csv')
print(df_blastn.shape)


# # Result of the blastn to bed file

# In[38]:


flanking_value = 200
df = df_blastn[['qseqid', 'sseqid', 'sstart', 'send', 'sstrand','slen']]
df['ones'] = 1


# In[39]:


def switch(row):
    if(row['sstart'] > row['send']):        
        temp = row['sstart']
        row['sstart'] = row['send']
        row['send'] = temp
    return row
df = df.apply(lambda row: switch(row), axis=1)


# In[40]:


def convert(inp):
    if(inp == "plus"):
        return "forward"
    if(inp == "minus"):
        return "reverse"
    raise Exception('Error, sstrand contains illegal word! only "plus" and "minus" are allowed')
df['strand'] = df['sstrand'].apply(lambda x: convert(x))


# In[41]:


def convert2sign(inp):
    if(inp == "plus"):
        return "+"
    if(inp == "minus"):
        return "-"
    raise Exception('Error, sstrand contains illegal word! only "plus" and "minus" are allowed')
df['sign'] = df['sstrand'].apply(lambda x: convert2sign(x))


# In[42]:


df['hit_length'] = df.apply(lambda row: abs(row['send'] - row['sstart']) + 1 ,axis=1)


# ## convert sstart and send from location to index (range)

# In[43]:


df['sstart'] = df['sstart'].apply(lambda x: x - 1)


# In[44]:


df['downstream_flanking'] = df['sstart'].apply(lambda x:  flanking_value if x > flanking_value else x)


# In[45]:


df['upstream_flanking'] = df.apply(lambda row:  flanking_value if (row['send']+flanking_value) <= row['slen'] else row['slen'] - row['send'],axis=1)


# In[46]:


df['hit_start'] = df.apply(lambda row: row['downstream_flanking'] if row['sign'] == "+" else row['upstream_flanking'],axis=1)


# In[47]:


df['hit_end'] = df.apply(lambda row: row['downstream_flanking'] + row['hit_length'] if row['sign'] == "+" else row['upstream_flanking'] + row['hit_length'],axis=1)


# In[48]:


df['sstart'] = df['sstart'].apply(lambda x: max(x - flanking_value, 0))
df['send'] = df.apply(lambda row: min(row['send'] + flanking_value , row['slen']),axis=1)


# In[49]:


df['tag'] = df.apply(lambda row: f">{row['sseqid']}:{row['sstart']}-{row['send']}({row['sign']})",axis=1)
df['reformated_tag'] = df['tag'].apply(lambda t: reformat(t))
df[['tag', 'reformated_tag', 'hit_start', 'hit_end']].to_csv(f'./{temp_path}/hit_index_info.csv')#, index=False)


# In[50]:


df['location_tag'] = df.apply(lambda row: f">{row['sseqid']}|{row['sign']}|{row['sstart'] + 1}-{row['send']}|{row['hit_start']+1}-{row['hit_end']}",axis=1)
df[['location_tag','qseqid']].to_csv(f'./{temp_path}/pipe_seprated_location_list.csv',index=False,sep='\t')


# In[51]:


df[['sseqid','sstart','send','strand','ones', 'sign']].to_csv(f'./{temp_path}/extension_index.bed', 
        index=False, header=False, sep="\t")


# # Extention
# 

# In[52]:


# !sudo apt-get install bedtools


# In[53]:


get_ipython().system('bedtools getfasta -fi {input_genome_path} -fo ./{temp_path}/extended_original.txt -s -bed ./{temp_path}/extension_index.bed')
get_ipython().system('rm input_genome.fna.fai')


# In[54]:


# todo: remove duplicated
'''
df = fasta_to_df("./Temp/extended.txt")
df = df.drop_duplicates(subset=['tag'], keep='first')
df_to_fasta(df,"./Temp/extended.txt")
len(df['tag'].unique())
''';


# # Convert hit region to upper case and other region to lower case

# In[55]:


ext = fasta_to_df(f'./{temp_path}/extended_original.txt')
info = pd.read_csv(f'./{temp_path}/hit_index_info.csv')
info['tag'] = info['tag'].apply(lambda x: x[1:])
print(info.shape)
info.head(2)


# In[56]:


ext = ext.sort_values(by=['tag']).reset_index()
ext['help_tag'] = ext.apply(lambda r: r['tag'] + str(r.name),axis=1)
del ext['tag']

info = info.sort_values(by=['tag']).reset_index()
info['help_tag'] = info.apply(lambda row: row['tag']+ str(row.name),axis=1)
def redefined_tag(row):
    tag = row['tag']
    [sstart, send] = tag.split(':')[-1].split('(')[0].split('-')
    sstart = int(sstart) + 1
    sign = tag.split('(')[-1].split(')')[0]    
    return f"{tag.split(':')[0]}|{sign}|{sstart}-{send}|{row['hit_start']+1}-{row['hit_end']}"
info['tag'] = info.apply(lambda row: redefined_tag(row),axis=1)
ext = pd.merge(ext,info,how='inner', on='help_tag')

def emphasis_hit(row):
    seq = list(row['data'].lower())            
    s = row['hit_start']
    e = row['hit_end']
    seq[s:e] = list(''.join(seq[s:e]).upper())    
    return ''.join(seq)
    
ext['data'] = ext.apply(lambda row: emphasis_hit(row),axis=1)
df_to_fasta(ext[['tag','data']],f"./{temp_path}/extended_modified.txt")


# # Protein coding elimination [Download nr]

# In[ ]:


get_ipython().system('wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz')


# # Protein coding elimination [Diamond]

# In[57]:


#!wget http://github.com/bbuchfink/diamond/releases/download/v2.0.13/diamond-linux64.tar.gz
#!tar xzf diamond-linux64.tar.gz


# In[ ]:


get_ipython().system('./diamond makedb --in ./NR/nr -d ./Temp/diamond_output')


# In[ ]:


get_ipython().system('./diamond blastx -d ./Temp/diamond_output.dmnd                   -q ./Temp/extended_modified.txt                   -o ./Temp/diamond_matches.tsv                   -p 22')


# In[58]:


dmn = pd.read_csv(f"./{temp_path}/diamond_matches.tsv", sep='\t', header=None)
dmn.columns = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'.split(' ')
coding_seq = dmn['qseqid'].unique()


# In[59]:


def clear(inp):
    if(inp[:9] == "reverse::"):
        return inp[9:]
    if(inp[:9] == "forward::"):
        return inp[9:]
    return inp
coding_seq = pd.Series(coding_seq).apply(lambda x : clear(x))

ext = fasta_to_df(f'./{temp_path}/extended_modified.txt')
print(f'total:      {ext.shape[0]}')
non_coding = ext[~ext['tag'].isin(coding_seq)]
print(f'non_coding: {non_coding.shape[0]}')
df_to_fasta(non_coding,f'./{temp_path}/extended_modified_non_coding.txt')
coding = ext[ext['tag'].isin(coding_seq)]
print(f'coding:     {coding.shape[0]}')
df_to_fasta(coding,f'./{temp_path}/extended_modified_coding.txt')


# # RNA 2d prediction

# ## Mfold

# In[129]:


'''
# installation
!wget http://www.unafold.org/download/mfold-3.6.tar.gz
!tar -xvf ./mfold-3.6.tar.gz; rm ./mfold-3.6.tar.gz
%cd ./mfold-3.6
!./configure
!make
!make install
%cd ..
!sudo apt install texlive-font-utils
''';


# In[62]:


counter = 0
base = f"./{result_path}/secondary_structure/mfold/"
get_ipython().system('rm -r {base}')
get_ipython().system('mkdir -p {base}')
df = fasta_to_df(f'./{temp_path}/extended_modified_non_coding.txt')

for index, row in df.iterrows():    
    tag = reformat(row['tag'])
    if(not os.path.exists(base + tag)):
        os.makedirs(base + tag)            
    with open(base + f"{tag}/SEQ.FASTA",'w') as file:
        file.write(f">{row['tag']}\n{row['data']}")
    counter += 1    
    #if(counter >= 100):
      #  break


# In[63]:


get_ipython().run_cell_magic('capture', '', 'remove_lock = False\ndef run_mfold(tag):\n    tag = reformat(tag)\n    %cd {base + tag}\n    !mfold  SEQ="SEQ.FASTA" T=22   \n    if(not remove_lock):\n        !find . -name "SEQ*" -not -name "*.ct" -not -name "*.pdf" -not -name "*SEQ.FASTA" -not -type d -delete\n    %cd {current_path}\n\nif __name__ == \'__main__\':        \n    pool = mp.Pool(mp.cpu_count() - 3)      \n    pool.map(run_mfold, df[\'tag\'])  ')


# ## Mxfold2

# In[ ]:


#!wget https://github.com/keio-bioinformatics/mxfold2/releases/download/v0.1.1/mxfold2-0.1.1.tar.gz
#!pip3 install mxfold2-0.1.1.tar.gz
#!rm mxfold2-0.1.1.tar.gz


# In[ ]:


get_ipython().system('mxfold2 predict ./extended.txt > Result/secondary_structure/mxfold2_result.txt')


# In[ ]:


df = fasta_to_df('./Result/secondary_structure/mxfold2_result.txt')
df = df.apply(lambda row: bracket_row(row) , axis=1)
df.head(2)


# In[ ]:


base = "./Result/secondary_structure/mxfold2/"
get_ipython().system('rm -r {base}')
get_ipython().system('mkdir -p {base}')
for index, row in df.iterrows():    
    if(not os.path.exists(base + reformat(row['tag']))):
        os.makedirs(base + reformat(row['tag']))        
    tag = reformat(row['tag'])
    with open(base + f"{tag}/{tag}.ct",'w') as file:
        bracket = row['bracket'].split(' ')[0]
        deltaG = row['bracket'].split(' ')[1]
        ct = bracket_to_ct(row['tag'], row['data'], bracket, deltaG)
        file.write(ct)    


# ## Vienna package

# In[ ]:


#!wget https://www.tbi.univie.ac.at/RNA/download/ubuntu/ubuntu_20_04/viennarna_2.4.18-1_amd64.deb -O viennarna.deb
#!sudo dpkg -i ./viennarna.deb
#!sudo apt-get -f install
#!rm viennarna.deb


# In[ ]:


base = "./Result/secondary_structure/viennarna/"
get_ipython().system('rm -r {base}')
get_ipython().system('rm ./Result/secondary_structure/viennarna_result.txt')
get_ipython().system('mkdir -p {base}')


# In[ ]:


get_ipython().run_line_magic('cd', '{base}')
get_ipython().system('RNAfold --jobs=0 --infile ../../Temp/extended_modified.txt  --noPS -T 22 > ../viennarna_result.txt')
get_ipython().run_line_magic('cd', '{current_path}')


# In[ ]:


df = fasta_to_df('./Result/secondary_structure/viennarna_result.txt')
df = df.apply(lambda row: bracket_row(row) , axis=1)
print(df.shape)
df.head(2)


# In[ ]:


for index, row in df.iterrows():    
    tag = reformat(row['tag'])
    if(not os.path.exists(base + tag)):
        os.makedirs(base + tag)      
    with open(base + f"{tag}/{tag}.ct",'w') as file:
        bracket = row['bracket'].split(' ')[0]
        deltaG = row['bracket'].split(' ')[1]
        ct = bracket_to_ct(row['tag'], row['data'], bracket, deltaG, False)
        file.write(ct)    


# In[ ]:


import glob
for file in glob.glob(f"{base}*.ps"):    
    f = file[len(base):-6] # _ss.ps 
    f = reformat(f)        
    shutil.move(file, f"{base}{f}/{f}.ps")    


# ## ContraFold

# In[ ]:


#!wget http://contra.stanford.edu/contrafold/contrafold_v2_02.tar.gz
#!tar -xvzf contrafold_v2_02.tar.gz && rm contrafold_v2_02.tar.gz
#%cd contrafold/src
#!make clean
#!make 
# to file must changed to be complieable # utility.hpp and optimization.c++ files


# In[ ]:


counter = 0
base = f"./{result_path}/secondary_structure/contrafold/"
get_ipython().system('rm -r {base}')
get_ipython().system('mkdir -p {base}')
df = fasta_to_df(f'./{temp_path}/extended.txt')

for index, row in tqdm(df.iterrows()):    
    tag = reformat(row['tag'])
    if(not os.path.exists(base + tag)):
        os.makedirs(base + tag)            
    with open(base + f"{tag}/{tag}.FASTA",'w') as file:
        file.write(f">{row['tag']}\n{row['data']}")
    counter += 1        


# In[ ]:


def run_contrafold(tag):
    tag = reformat(tag)    
    get_ipython().run_line_magic('cd', 'Software/contrafold/src')
    get_ipython().system('./contrafold predict ../..{base[1:]}{tag}/{tag}.FASTA > ../..{base[1:]}{tag}/{tag}.dot')
    with open(f"../..{base[1:]}{tag}/{tag}.dot", 'r') as file:
        text = file.read()
    text = [l for l in text.split("\n") if l[:len(">structure")] != ">structure"]    
    header = text[0]
    with open(f"../..{base[1:]}{tag}/{tag}.dot", 'w') as file:
        file.write('\n'.join(text[1:]))    
    get_ipython().system('RNAeval  ../..{base[1:]}{tag}/{tag}.dot -T 20 > ../..{base[1:]}{tag}/{tag}.dotdg    ')
    with open(f"../..{base[1:]}{tag}/{tag}.dotdg", 'r') as file:
        text = file.read()
    with open(f"../..{base[1:]}{tag}/{tag}.dot", 'w') as file:
        file.write(header + "\n" + text)    
    
    df = fasta_to_df(f'../..{base[1:]}{tag}/{tag}.dot')
    df = df.apply(lambda row: bracket_row(row) , axis=1)        
    tag = reformat(df['tag'][0])
    with open(f'../..{base[1:]}{tag}/{tag}.ct','w') as file:
        bracket = df['bracket'][0].split(' ')[0]        
        deltaG = df['bracket'][0].split(' ')[1]
        ct = bracket_to_ct(df['tag'][0], df['data'][0], bracket, deltaG, False)
        file.write(ct)    
    #!rm ../..{base[1:]}{tag}/{tag}.dot
    #!rm ../..{base[1:]}{tag}/{tag}.dotdg
    get_ipython().system('rm ../..{base[1:]}{tag}/{tag}.FASTA')
    get_ipython().run_line_magic('cd', '{current_path}')

if __name__ == '__main__':        
    pool = mp.Pool(mp.cpu_count() - 1)  
    pool.map(run_contrafold, df['tag'].iloc[:10])


# In[ ]:


'''path = 'secondary_structure/contrafold/AMWY020333941_469-893_-_/AMWY020333941_469-893_-_.dot'
!RNAeval  {path} -T 20 -v'''; 


# # CTAnalizer

# In[64]:


# only select those not ran before
base = f"./{result_path}/secondary_structure/mfold/"
df = fasta_to_df(f'./{temp_path}/extended_modified_non_coding.txt')
index_list =[]
for index, row in df.iterrows():    
    tag = reformat(row['tag'])    
    if(len(glob.glob(f'{base + tag}/*.ct')) != 0):
        index_list.append(index)
df = df.iloc[index_list,:]
print(df.shape)


# In[65]:


def run(tag, path, extra):                    
    try:
        return get_row(tag, path,extra)
    except Exception as e:
        print(str(e), tag)        
        return pd.Series()
        
def get_df_by_tag(tag , extra=0):           
    ct_files = glob.glob(f'{base}{reformat(tag)}/SEQ_*.ct')    
    return pd.Series(ct_files).apply(lambda path: run(tag, path,extra))    


# In[66]:


d = get_df_by_tag("NC_050105.1|-|1-225|201-224", extra=0)
#d.to_csv('Result/d.csv')#.iloc[0,:]['full seq visualization']


# # Apply on current data

# In[67]:


seq2cluster = pd.read_csv(f"{temp_path}/seq2cluster.csv")
seq2cluster['tag'] = seq2cluster.groupby(['cluster'])['tag'].transform(lambda x: ','.join(x))
seq2cluster['seqid'] = seq2cluster.groupby(['cluster'])['seqid'].transform(lambda x: ','.join(x))
seq2cluster = seq2cluster.drop_duplicates()
tag2cluster = pd.read_csv(f'./{temp_path}/pipe_seprated_location_list.csv',sep='\t')
tag2cluster['location_tag'] = tag2cluster['location_tag'].apply(lambda x : x[1:])
data = pd.merge(seq2cluster,tag2cluster,how='inner', left_on='cluster', right_on='qseqid')
data['Reference miRNA cluster'] = data['cluster']
data['Reference miRNA IDs'] = data['seqid']
data['Reference miRNA IDs and species'] = data['tag']
data = data[['location_tag','Reference miRNA cluster', 'Reference miRNA IDs', 'Reference miRNA IDs and species','confidence']]


# In[68]:


rcols_ref = ['Reference miRNA cluster',
             'Reference miRNA IDs',
             'Reference miRNA IDs and species']
rcols_boi = ['boi seq', 'boi name', 'boi dotbracket']

rcols = [*rcols_ref,
         *rcols_boi]

rcols_dg = [*rcols, 'delta G']

def selection(row):        
    global repeted
    tuple_row = tuple(row)            
    if(tuple_row not in repeted):
        repeted[tuple_row] = row.name
        return True
    return False

def boi_selection(row):    
    global repeted_boi        
    tuple_row_boi = tuple(row[rcols_boi])
    dg = row['delta G']
    if(tuple_row_boi not in repeted_boi):
        repeted_boi[tuple_row_boi] = { 
            "counter": 1,
            "dg": dg,
            "ref clusters": row['Reference miRNA cluster'],
            "ref ids": row['Reference miRNA IDs'],
            "ref species": row['Reference miRNA IDs and species'],
            "lock": False
        }
    else:            
        value =  repeted_boi[tuple_row_boi]
        value['counter'] += 1
        value['dg'] = min(value['dg'] , dg)
        value['ref clusters'] += "," + row['Reference miRNA cluster']
        value['ref ids'] += "," + row['Reference miRNA IDs']
        value['ref species'] += "," +  row['Reference miRNA IDs and species']                
        repeted_boi[tuple_row_boi] = value        


# In[69]:


get_ipython().system('rm ./{result_path}/ct_analizer.csv')
chunksize = 1 * (10 ** 4)
max_workers = mp.cpu_count() - 4
num_terminal = 5 # acceptable_terminal_structures

repeted = {}
repeted_boi = {}
header = True
orders = None
arr = np.array_split(df['tag'], max(df['tag'].shape[0]//chunksize , 1))
for chunk in tqdm(arr):
    dfs = []
    for row in process_map(get_df_by_tag , chunk, tqdm_class=tqdm, max_workers=max_workers, chunksize=5):        
        dfs.append(row)
    chunk = pd.concat(dfs,axis=0)
    chunk = pd.merge(data, chunk, how='right', left_on = 'location_tag', right_on ='seq name')
    del chunk['location_tag']
    if(header):
        orders = chunk.columns
    for col in orders:
        if(col not in chunk.columns):
            chunk[col] = np.nan            
        
    for col in chunk.columns:
        if(col not in orders ):
            print(f"Error in {col}")        
    chunk = chunk.reindex(columns=orders)
    chunk = chunk.replace(np.nan, '-').replace('', '-')                    
    # delete repeated
    selected = chunk[rcols].apply(lambda row: selection(row), axis=1)
    chunk = chunk[selected]
    
    # cluster refs
    chunk[rcols_dg].apply(lambda row: boi_selection(row), axis=1)        
    chunk.to_csv(f"./{result_path}/ct_analizer.csv", header=header, mode='a', index=False)    
    header = False


# In[70]:


get_ipython().system('rm ./{result_path}/ct_analizer_clustered.csv')
def isKeepCluster(row):
    global repeted_boi    
    dg = row['delta G']    
    if(row["boi name"] == "-"):
        return False
    tuple_row_boi = tuple(row[rcols_boi])    
    value = repeted_boi[tuple_row_boi]
    if(value['counter'] == 1):
        return True
    if(value['dg'] != dg):        
        return False
    if(value['lock']):
        return False
    value['lock'] = True
    repeted_boi[tuple_row_boi] = value
    return True


def makeCluster(row):            
    tuple_row_boi = tuple(row[rcols_boi])    
    value = repeted_boi[tuple_row_boi]
    if(value['counter'] != 1):
        for ref_c in ['ref clusters', 'ref ids', 'ref species']:
            value[ref_c] = value[ref_c].replace(' ,', ',').replace(', ', ',')                        
            value[ref_c] = value[ref_c].split(',')
            value[ref_c] = set(value[ref_c])
            value[ref_c] = ", ".join(value[ref_c])                    
        row['Reference miRNA cluster'] = value["ref clusters"]
        row['Reference miRNA IDs'] = value["ref ids"]
        row['Reference miRNA IDs and species'] = value["ref species"]
    return row


header = True
for chunk in tqdm(pd.read_csv(f"./{result_path}/ct_analizer.csv", chunksize=10 ** 5)):
    chunk = chunk[chunk[rcols_dg].apply(lambda row:  isKeepCluster(row), axis=1)]
    chunk = chunk.apply(lambda row : makeCluster(row), axis=1)
    chunk.to_csv(f"./{result_path}/ct_analizer_clustered.csv", mode='a', index=False)
    header = False


# # Filters

# In[71]:


get_ipython().system('rm ./{result_path}/result_level1_filter.csv')
filter1_run(input_file=  f"./{result_path}/ct_analizer_clustered.csv",
            output_file= f"./{result_path}/result_level1_filter.csv")


# In[107]:


config = {'delta_g_min': -999,
          'delta_g_max': 1,
          'hit_len_min': 21,
          'hit_len_max': 21,
          'hit_complementarity_percentage_min': 0.5,
          'hit_complementarity_percentage_max': 1.0,
          'number_of_terminal_structure_min': 0,
          'number_of_terminal_structure_max': 5,
          'boi_gc_content_min': 45,
          'boi_gc_content_max': 94, 
          'num_of_linking_residues_min': 5, 
          'num_of_linking_residues_max': 159,
          'hit_gc_content_percentage_min': 37, 
          'hit_gc_content_percentage_max': 86,
          'precursor_mfei_min': 0.87,
          'precursor_mfei_max': 1.3695556794836854, 
          'border_line_mismatch_max': 0,
          'border_line_bulge_max': 0,
          'border_line_internal_max': 0,
          'total_num_of_nonmatching_positions': 5,
          'total_num_of_mismached_positions': 5,
          'total_num_of_positions_in_bulges_and_loops': 2,
          'max_allowed_mismatch_size_in_hit_region': 2,
          'max_allowed_bulge_size_in_hit_region': 1,
          'max_allowed_internal_loop_size_in_hit_region': 3,
          'max_allowed_hsbl_ssbl_size': 2,
          'minimum_required_clear_region': 0,
          'acceptable_num_for_hit_locations_in_bulges_or_loops': 2,
          'acceptable_num_for_unmatched_locations_in_hit_region': 5,
          'delete_if_mature_duplex_involvement_in_apical_loop': 'YES',
          'border_line_structure_allowance': 'NOT ACCEPTED'}

filter2(input_file  = f"./{result_path}/result_level1_filter.csv",
            output_file = f"./{result_path}/result_level2_filter.csv",
            config=config)


# # Cluster JSC

# In[20]:


result = pd.read_csv(f"./{result_path}/result_level2_filter.csv")
print(result.shape)
result.head(2)


# In[21]:


def jaccard(A, B):
    a1 = int(A.split('-')[0])
    a2 = int(A.split('-')[1])
    b1 = int(B.split('-')[0])
    b2 = int(B.split('-')[1])
    s1 = set([i for i in range(min(a1,a2), max(a1,a2) + 1)])
    s2 = set([i for i in range(min(b1,b2), max(b1,b2) + 1)])    
    intersection = len(s1.intersection(s2))
    union = (len(s1) + len(s2)) - intersection
    return float(intersection) / union


# In[22]:


def same2cluster(same_dict, threshold=0.8):        
    counter = 1
    for key in same_dict:        
        item2cluster = {}
        SET = list(set(same_dict[key]))
        G = networkx.Graph()
        # add nodes
        for s in SET:
            G.add_node(s)
        # add edges
        for i in range(0, len(SET)):
            for j in range(i+1, len(SET)):
                if(jaccard(SET[i],SET[j]) >= threshold):        
                    G.add_edge(SET[i], SET[j], weight=jaccard(SET[i],SET[j]))
                            
        # get maximal
    
        for clique in maximal_cliques(G):    
            unique = str(counter).zfill(4)
            counter += 1
            for item in clique:
                if(item in item2cluster):
                    item2cluster[item].append(unique)
                else:
                    item2cluster[item] = [unique]    
        same_dict[key] = item2cluster
    return same_dict 


# ## hit jaccard similarity

# In[23]:


hit_threshold = 0.8
same_strand_hit = {}

def same_strand(row):    
    global same_strand_hit
    chrom = row['chromosome']
    sign = row['sign']
    hit = row['hit position on chromosome']      
    key = f'{chrom}{sign}'
    if(key in same_strand_hit):        
        same_strand_hit[key].append(hit)
    else:
        same_strand_hit[key] = [hit]
        
rcols_hit = ['chromosome', 'sign', 'hit position on chromosome']
result[rcols_hit].apply(lambda row: same_strand(row), axis=1)
same_strand_hit = same2cluster(same_strand_hit, threshold=hit_threshold)

def _f(row):        
    key = f"{row['chromosome']}{row['sign']}"
    return same_strand_hit[f'{key}'][row['hit position on chromosome']]
result['hit cluster number'] = result[rcols_hit].apply(lambda row: _f(row),axis=1)


# In[24]:


rcols_boi = ['boi seq', 'boi name', 'boi dotbracket']
boi_threshold = 0.8
same_strand_boi = {}

def same_strand(row):    
    global same_strand_boi    
    seq = row['boi seq'].lower()
    name = row['boi name']
    dotbracket = row['boi dotbracket'].lower()
    if(pd.isna(name) or name == "-"):
        return 
    name = name.split('|')    
    loction = name[0] + name[1] + name[2]
    hit = name[3]
    key = f'{seq}{dotbracket}{loction}'
    if(key in same_strand_boi):        
        same_strand_boi[key].append(hit)
    else:
        same_strand_boi[key] = [hit]
        
        
result[rcols_boi].apply(lambda row: same_strand(row), axis=1)

same_strand_boi = same2cluster(same_strand_boi, threshold=boi_threshold)

def _f(row):        
    seq = row['boi seq'].lower()
    name = row['boi name']
    dotbracket = row['boi dotbracket'].lower()
    if(pd.isna(name) or name == "-"):
        return 
    name = name.split('|')    
    loction = name[0] + name[1] + name[2]
    hit = name[3]
    key = f'{seq}{dotbracket}{loction}'        
    return same_strand_boi[f'{key}'][hit]
result['boi cluster number'] = result[rcols_boi].apply(lambda row: _f(row),axis=1)


# In[25]:


precursor_threshold = 0.8
same_strand_precursor = {}
rcols_pre = ['precursor seq', 'precursor name', 'precursor dotbracket']

def same_strand(row):    
    global same_strand_precursor    
    seq = row['precursor seq'].lower()
    name = row['precursor name']
    dotbracket = row['precursor dotbracket'].lower()
    if(pd.isna(name) or name == "-"):
        return 
    name = name.split('|')    
    loction = name[0] + name[1] + name[2]
    hit = name[3]
    key = f'{seq}{dotbracket}{loction}'
    if(key in same_strand_precursor):        
        same_strand_precursor[key].append(hit)
    else:
        same_strand_precursor[key] = [hit]
        
        
result[rcols_pre].apply(lambda row: same_strand(row), axis=1)
same_strand_precursor = same2cluster(same_strand_precursor, threshold=precursor_threshold)

def _f(row):        
    seq = row['precursor seq'].lower()
    name = row['precursor name']
    dotbracket = row['precursor dotbracket'].lower()
    if(pd.isna(name) or name == "-"):
        return 
    name = name.split('|')    
    loction = name[0] + name[1] + name[2]
    hit = name[3]
    key = f'{seq}{dotbracket}{loction}'
    return same_strand_precursor[f'{key}'][hit]
result['precursor cluster number'] = result[rcols_pre].apply(lambda row: _f(row),axis=1)


# In[26]:


hit2cluster = {}
hit_unique = result['hit seq'].unique()
for i in range(0, hit_unique.shape[0]):
    hit2cluster[hit_unique[i]] = str(i + 1).zfill(4)
result['identical hit cluster'] = result['hit seq'].apply(lambda hit: hit2cluster[hit])


# In[27]:


seed_start = 2
seed_end = 13
result['seed region'] = result['hit seq'].apply(lambda hit: hit[seed_start-1:seed_end])


# In[28]:


result.to_csv(f"./{result_path}/result_level2_filter_clustered.csv",index=False)
get_ipython().system('zip -r ./{result_path}/result_level2_filter_clustered.zip ./{result_path}/result_level2_filter_clustered.csv')


# # BlastX

# In[ ]:


get_ipython().system('makeblastdb -in ./NR/nr -dbtype prot -out ./NR/nr_database')

#!head -n 100 ./Temp/extended_modified.txt > ./input_blastx.txt

get_ipython().system('blastx -query ./input_blastx.txt         -db ./NR/nr_database         -out ./Temp/BlastX/blastx         -num_threads 20         -evalue 1e-3         -outfmt "6 qseqid sseqid qstart qend evalue bitscore score length frames qframe qcovs qcovhsp staxids"')

blx = pd.read_csv('./Temp/BlastX/blastx', sep='\t', header=None)
blx.columns = 'qseqid sseqid qstart qend evalue bitscore score length frames qframe qcovs qcovhsp staxids'.split(' ')
coding_seq = blx['qseqid'].unique()


# # psRNATarget

# In[153]:


url = "https://www.zhaolab.org/psRNATarget/analysis"

payload='allowbulge=yes&curschema=s2&cutpos1=10&cutpos2=11&expect=5&function=3&gapextp=0.5&gapstartp=2&gup=0.5&hspsize=19&maxnummismatchinseed=2&misp=1&seedfactor=1.5&seedpos1=2&seedpos2=13&srna=&srna_content=%3Eath-miR156a%0D%0AUGACAGAAGAGAGUGAGCAC%0D%0A%3Eath-miR157a%0D%0AUUGACAGAAGAUAGAGAGCAC%0D%0A%3Eath-miR158a%0D%0AUCCCAAAUGUAGACAAAGCA%0D%0A%3Eath-miR398a%0D%0AUGUGUUCUCAGGUCACCCCUU%0D%0A%3Eath-miR398b%0D%0AUGUGUUCUCAGGUCACCCCUG%0D%0A%3Eath-miR398c%0D%0AUGUGUUCUCAGGUCACCCCUG%0D%0A%3Eath-miR834%0D%0AUGGUAGCAGUAGCGGUGGUAA%0D%0A%3Eath-miR390a%0D%0AAAGCUCAGGAGGGAUAGCGCC%0D%0A%3Eath-miR390b%0D%0AAAGCUCAGGAGGGAUAGCGCC&srna_uploaded=&target=&target_content=%3EAT1G27360.1%0D%0AAAGGTATCTATTTGCCTAGCCAGAGTTATATATAGGATTGATTGTCTAGTCTTTTCTTATATGATTTTTGTTCTCATTTACTAATCAAAGTTCTGCAAACTTGTAGTTGTTGTAGGATTTGTTGCTCTGGCTCTGGTGGTAGGTCTATGAAATCAACCCATATCGTGAATGGACTGCAACATGGTATCTTCGTCCCAGTGGGATTGGGAGCATTTGATCATGTCCAATCCGTCAAGGACTGAAGATGACAGCAAACAGCTACCTACTGAGTGGGAAATTGAAAAAGGTGAAGGAATTGAATCTATAGTTCCACATTTCTCAGGCCTTGAGAGAGTCAGTAGTGGCTCTGCCACCAGCTTCTGGCACACTGCTGTATCGAAAAGCTCACAGTCGACCTCTATCAACTCATCATCTCCCGAAGCCAAACGATGCAAGCTTGCATCAGAAAGTTCCCCTGGAGATTCTTGCAGCAACATAGACTTTGTCCAGGTGAAGGCTCCCACAGCTCTCGAGGTATCCGTTGCCTCAGCTGAATCAGATCTTTGTTTAAAACTAGGAAAGCGGACATACTCTGAAGAATACTGGGGTAGAAACAATAATGAAATTTCAGCGGTTTCTATGAAGTTGTTAACTCCATCTGTTGTCGCTGGGAAATCCAAATTGTGTGGTCAGAGCATGCCAGTCCCGCGTTGCCAAATTGATGGCTGTGAACTGGATCTCTCATCTGCTAAGGGTTATCATCGTAAGCACAAAGTCTGCGAAAAGCATTCAAAGTGCCCAAAAGTTAGCGTGAGTGGCCTGGAACGTCGGTTCTGCCAACAGTGTAGCAGGTTCCATGCTGTCTCTGAATTTGATGAGAAGAAACGAAGCTGCCGAAAACGTCTTTCTCATCATAATGCGAGGCGTCGTAAGCCACAAGGAGTATTTTCAATGAATCCCGAGAGGGTGTATGATCGAAGACAGCATACAAATATGTTGTGGAATGGGGTGTCCCTTAACGCGAGATCTGAAGAAATGTATGAATGGGGTAATAACACTTATGATACAAAGCCTAGACAAACGGAAAAAAGCTTTACTCTGAGCTTCCAGAGAGGTAATGGCTCTGAGGACCAGCTGGTTGCTAGTAGCAGCCGTATGTTCTCTACATCTCAAACCTCAGGTGGGTTCCCAGCAGGAAAGTCCAAGTTTCAACTTCATGGCGAAGATGTGGGAGAATACTCAGGAGTCCTCCATGAATCTCAAGATATCCACCGTGCTCTCTCTCTTCTGTCAACCTCTTCGGATCCCCTGGCCCAACCACATGTGCAGCCATTTTCTCTACTCTGTTCATATGATGTTGTACCAAAATAGATGAGTAAGTAATGTGTAATTTGTAAACCTGTTACTCAGTTGGTGGATACTTTTCCAAACCTATGATAAAAACCTCGTCCTAGATCCCGTTAAATGCCAAACTTTCGGCTACTATAACTATGTTATCGTTATCATTATCATTGTTTAACACCCT%0D%0A%3EAT1G27360.4%20%7C%20Symbols%3A%20%20%7C%20squamosa%20promoter-binding%20protein-like%2011%20(SPL11)%20%7C%20chr1%3A9501808-9503856%20FORWARD%0D%0ACTGGGTGAAACATAGAAAAGTTTCTCTTGCTCAAGTTAATGATAAAAGGGTGAGAGCAATAAACGCTGATAAGCCTTGTCTGGTCCTTGGAATTTTGAATTTTCTTTTTCTATCTTACTTATAGTATTGGTAGTTGAGGGTGTCGTCGATAAGTTGTTGTAGGATTTGTTGCTCTGGCTCTGGTGGTAGGTCTATGAAATCAACCCATATCGTGAATGGACTGCAACATGGTATCTTCGTCCCAGTGGGATTGGGAGCATTTGATCATGTCCAATCCGTCAAGGACTGAAGATGACAGCAAACAGCTACCTACTGAGTGGGAAATTGAAAAAGGTGAAGGAATTGAATCTATAGTTCCACATTTCTCAGGCCTTGAGAGAGTCAGTAGTGGCTCTGCCACCAGCTTCTGGCACACTGCTGTATCGAAAAGCTCACAGTCGACCTCTATCAACTCATCATCTCCCGAAGCCAAACGATGCAAGCTTGCATCAGAAAGTTCCCCTGGAGATTCTTGCAGCAACATAGACTTTGTCCAGGTGAAGGCTCCCACAGCTCTCGAGGTATCCGTTGCCTCAGCTGAATCAGATCTTTGTTTAAAACTAGGAAAGCGGACATACTCTGAAGAATACTGGGGTAGAAACAATAATGAAATTTCAGCGGTTTCTATGAAGTTGTTAACTCCATCTGTTGTCGCTGGGAAATCCAAATTGTGTGGTCAGAGCATGCCAGTCCCGCGTTGCCAAATTGATGGCTGTGAACTGGATCTCTCATCTGCTAAGGGTTATCATCGTAAGCACAAAGTCTGCGAAAAGCATTCAAAGTGCCCAAAAGTTAGCGTGAGTGGCCTGGAACGTCGGTTCTGCCAACAGTGTAGCAGGTTCCATGCTGTCTCTGAATTTGATGAGAAGAAACGAAGCTGCCGAAAACGTCTTTCTCATCATAATGCGAGGCGTCGTAAGCCACAAGGAGTATTTTCAATGAATCCCGAGAGGGTGTATGATCGAAGACAGCATACAAATATGTTGTGGAATGGGGTGTCCCTTAACGCGAGATCTGAAGAAATGTATGAATGGGGTAATAACACTTATGATACAAAGCCTAGACAAACGGAAAAAAGCTTTACTCTGAGCTTCCAGAGAGGTAATGGCTCTGAGGACCAGCTGGTTGCTAGTAGCAGCCGTATGTTCTCTACATCTCAAACCTCAGGTGGGTTCCCAGCAGGAAAGTCCAAGTTTCAACTTCATGGCGAAGATGTGGGAGAATACTCAGGAGTCCTCCATGAATCTCAAGATATCCACCGTGCTCTCTCTCTTCTGTCAACCTCTTCGGATCCCCTGGCCCAACCACATGTGCAGCCATTTTCTCTACTCTGTTCATATGATGTTGTACCAAAATAGATGAGTAAGTAATGTGTAATTTGTAAACCTGTTACTCAGTTGGTGGATACTTTTCCAAACCTATGATAAAAACCTCGTCCTAGATCCCGTTAAATGCCAAACTTTCGGCTACTATAACTATGTTATCGTTATCATTATCATTGTTTAACACCCT%0D%0A%3EAT1G32140.1%20%7C%20Symbols%3A%20%20%7C%20F-box%20family%20protein%20%7C%20chr1%3A11562722-11564813%20REVERSE%0D%0AATGACGATGATGTCCGACCTTTCACTTGATTTAGTCGAAGAGATATTGTGTAGGGTTCCGATAACTTCTCTTAAAGCAGTGAGATCTAGTTGCAAACTATGGAACGTTCTTTCCAAGAACCGGATTTTATGTAAAACAGAAGCTAGAAATCAGTTTTTAGGGTTCACGATAATGAATCATAGGCTTTATTCCATGAGATTCAATCTCCATGGAATCGGCCTCAATGAAAACAGTGAAGAGTTCATTGATCCATCTATAAAGCCAATAGGTAATTTACTTAATCAAGTCGAGATATCTAAAGTGTTTTATTGCGAAGGTTTATTGTTATGCGTCACAAGGAACCACTCAAGCAAGCTCGTAGTTTGGAACCCGTATTTGGGAGAAATTCGTTGGATCAAAACTAGGAATGATTACCACATAGGCGTTACATATGCTCTCGGGTACGACAACAACAAGAACCACATGATCTTGAGGTTTTTTTCTGAACAAGGCTACTACGAGATTTACGACATGAACTCTTCTGACTCATGGGATTGTTTTTATGGCATTCCCAACAAGGGGTTAAAATGTTATCAGCCCGGCGCGTCGTTAAATGGAAATGCTTATTTTTTGACTGAGGGAAGAGAAGTAATGGAAGGGTATGATTGCTTACTCGGTTTTGATTTTACAACAAAGAAATTTGGACCACTTCTTTCTTTGTCGTTTTCGCATGATTTTATAGAGACTGGGAGACTATCTTGTGTTAAAGGAGAGAAACTTGCGGTCTTATATCAGCGCTGCTATACCTATGAGATGGAGGTATGGGTGACAACTAAGATAGAGCCGAATGCGGTGTCATGGAGCAAATTCTTAGCAGTTGAAATGGAACCACTCACTAGTCTAAAGTTTAACGATGATTCTGGCAGTTTCTTCATTGACGAAGAGAAGAAAATCGTCGTGGTTTTTGATATAGACGAATCTGAACGCAACAATACGGCTTACATCATTGGAGATTATGGATGCTTGAAAGAAGTTGATCTTGATGAAGTTGTGAACCCACAAGAATCTGTGGAGGTCGGAGACCGCATTTATTCTTTTTCACCATTTGTGTGCTCTTGCTCTTATGTTCCAAGTTTAGTGAAATTTAAAGAAGATGCAGAACATGAAAGGAAAGATAAGAAGAGGAAGAGTAAGAGGAAGCGAACCAACAAGGATGGATATGATTTTATACTGTGTTTTGATTTTACAACCGAGAGATTTGGACAGATTCTTCCTCTGCCGTTTAAACATTCTTTTAGGGATACTTGGACTCTATCTTCTGTTAAAGAAGAGAAACTTGCAGTCGCAGTGTTATACTGGAAAAATACATGTGTGATGATAGAGATATGGATGACAATTAAGATTGATCCTAATGTCGAGTCGTGGAGCAAATTCTTAAGAGTTGATAGGAAACCATGCATTGATCTTCGCTTTGATGATCGTAATGACAGTTTTTTCATTGACGAAGAGAAGAAAGTTGTCGTGTTTTTTAGTTCAGACAAAGTTAAAACCTCTACGGCTTACGTCATTGGAGATAATAGATATTTGAGAACAGTGGATCTTGAAAAAGCTGCAAACTCCCAAGAATCTGTGGAGGTCGGAGAACGCGTGTATTGTTTTTCGCCGCTTGTGTGCTCTTGCTCTTATTATGTTCCAAGTTTGGTGAAAATCAACCACAATGCAGGACGCAAAAGGAAAGAGAAGAAGACGAAGCGCAAAAGTAAAGACAAGCAGATGAAACTAAGCAACAAGGTGTAA%0D%0A%3EAT2G42200.1%20%7C%20Symbols%3A%20%20%7C%20squamosa%20promoter-binding%20protein-like%209%20(SPL9)%20%7C%20chr2%3A17594485-17596708%20FORWARD%0D%0AACCACTCTCGTCTCTTTCTTTTTTCCTTCTGTTCTGTTTCTCTCTCTAAACCCAAAACAGTCAAAATCAGGGAAGCCGAAATTTTCTTTGCTTTCTTCTCCTTTGGTCCTTTCTTTAAACCCGAGACAGTTAGGTTTGTGTGAGAGAGAGAATGATGAGTAAAACCCTTTCTGTCTGAGTAAGAGGAAACCAACATGGAGATGGGTTCCAACTCGGGTCCGGGTCATGGTCCGGGTCAGGCAGAGTCGGGTGGTTCCTCCACTGAGTCATCCTCTTTCAGTGGAGGGCTCATGTTTGGCCAGAAGATCTACTTCGAGGACGGTGGTGGTGGATCCGGGTCTTCTTCCTCAGGTGGTCGTTCAAACAGACGTGTCCGTGGAGGCGGGTCGGGTCAGTCGGGTCAGATACCAAGGTGCCAAGTGGAAGGTTGTGGGATGGATCTAACCAATGCAAAAGGTTATTACTCGAGACACCGAGTTTGTGGAGTGCACTCTAAAACACCTAAAGTCACTGTGGCTGGTATCGAACAGAGGTTTTGTCAACAGTGCAGCAGGTTTCATCAGCTTCCGGAATTTGACCTAGAGAAAAGGAGTTGCCGCAGGAGACTCGCTGGTCATAATGAGCGACGAAGGAAGCCACAGCCTGCGTCTCTCTCTGTGTTAGCTTCTCGTTACGGGAGGATCGCACCTTCGCTTTACGAAAATGGTGATGCTGGAATGAATGGAAGCTTTCTTGGGAACCAAGAGATAGGATGGCCAAGTTCAAGAACATTGGATACAAGAGTGATGAGGCGGCCAGTGTCGTCACCGTCATGGCAGATCAATCCAATGAATGTATTTAGTCAAGGTTCAGTTGGTGGAGGAGGGACAAGCTTCTCATCTCCAGAGATTATGGACACTAAACTAGAGAGCTACAAGGGAATTGGCGACTCAAACTGTGCTCTCTCTCTTCTGTCAAATCCACATCAACCACATGACAACAACAACAACAACAACAACAACAACAACAACAACAACAATACATGGCGAGCTTCTTCAGGTTTTGGCCCGATGACGGTTACAATGGCTCAACCACCACCTGCACCTAGCCAGCATCAGTATCTGAACCCGCCTTGGGTATTCAAGGACAATGATAATGATATGTCTCCTGTTTTGAATTTAGGTCGATACACCGAGCCAGATAATTGTCAGATAAGTAGTGGCACGGCAATGGGTGAGTTCGAGTTATCTGATCACCATCATCAAAGTAGGAGACAGTACATGGAAGATGAGAACACAAGGGCTTATGACTCTTCTTCTCACCATACCAACTGGTCTCTCTGACTTGTCTTTGCATCAGAGAATCTTCTTACAATGAACGATTCTGCAATATCTTATCTTTTTGCTTCTTTGTTTATTCTGTTATCTGCTATCAATAAACCAGACAATTGTTGCCAGATAATGGCTTTTGATTTTGATTTGTTGTTTTATCTCCATGAAAATCCAAGTTATGAGATCAGATT%0D%0A%3EAT1G49910.1%20%7C%20Symbols%3A%20%20%7C%20WD-40%20repeat%20family%20protein%20%2F%20mitotic%20checkpoint%20protein%2C%20putative%20%7C%20chr1%3A18482693-18485143%20FORWARD%0D%0AATGACTTTGGTGCCGGCCATTGGTCGCGAGCTCTCGAATCCACCGTCCGATGGGATTTCTAATCTTCGATTTTCTAATAACAGCGATCATTTACTAGTCTCTTCATGGGATAAGAGTGTAAGATTGTATGATGCGAACGGCGATTTGATGAGAGGGGAGTTTAAACATGGTGGAGCGGTACTCGATTGCTGCTTCCATGATGATTCTTCTGGATTCAGTGTTTGCGCCGATACTAAAGTTAGAAGAATTGACTTCAATGCTGGCAAAGAAGACGTTTTGGGTACGCATGAGAAGCCAGTTCGATGTGTTGAGTATTCTTATGCTGCAGGGCAAGTGATCACTGGAAGTTGGGATAAAACGATTAAATGTTGGGATCCAAGAGGTGCAAGTGGGACTGAGCGCACACAGATTGGAACTTATATGCAACCTGAGCGTGTTAACTCTCTTTCTCTTGTTGGAAATCGTTTGGTAGTGGCAACAGCAGGAAGGCATGTCAACATTTATGATCTTAGAAATATGTCCCAGCCTGAGCAAAGAAGAGAGTCCTCACTTAAATACCAGACAAGATGTGTACGTTGTTATCCCAACGGAACAGGATATGCCCTTAGCTCTGTTGAAGGGAGGGTTTCAATGGAGTTTTTTGATCTATCAGAAGCTGCTCAGGCTAAAAAATATGCTTTCAAATGTCACCGGAAATCAGAGGATGGAAGGGACATTGTCTACCCTGTAAATGCAATTGCTTTCCATCCGATTTATGGCACTTTTGCTTCCGGAGGCTGTGATGGTTTTGTCAACATTTGGGACGGTAACAATAAGAAGAGGCTTTATCAGTACTCTAAGTATCCAACAAGTATTGCGGCGCTGTCATTCAGCCGAGATGGTGGATTACTGGCTGTTGCTTCTAGTTACACGTTTGAAGAGGGAGACAAACCGCATGAACCGGACGCCATCTTTGTTAGAAGTGTTAATGAAATTGAAGTGAAACCGAAACCCAAAGTATACCCAAATCCCCCGGTATAGTCAAGAAATAATGGAATGAGCAGAGTCAAATTCGACTTGTGTGTTGTTGTATTGTAGCACTTGAAAGTGAGTTATAAAATCTTATTTTGGCTGTAAAGTGAAATGTGAACGTTATAATGGCTTTCGAATCTGAGATGGTGTTCCATTTACTCTCTGGGTTGCCTCCAATTTTTCTTTTAGGACCAACTATCTTATTTTTACCTT%0D%0A%3EAT1G48460.1%20%7C%20Symbols%3A%20%20%7C%20similar%20to%20unknown%20protein%20%5BArabidopsis%20thaliana%5D%20(TAIR%3AAT5G63040.2)%3B%20similar%20to%20Os01g0704200%20%5BOryza%20sativa%20(japonica%20cultivar-group)%5D%20(GB%3ANP_001044004.1)%3B%20similar%20to%20hypothetical%20protein%20MtrDRAFT_AC124952g33v1%20%5BMedicago%20truncatula%5D%20(GB%3AABE93586.1)%3B%20contains%20domain%20Multidrug%20resistance%20ABC%20transporter%20MsbA%2C%20N-terminal%20domain%20(SSF90123)%20%7C%20chr1%3A17915014-17916968%20FORWARD%0D%0AAGAAAACAATTCCAAAAAATAAAATGTCAGAAAAGAATTTTCTTTTAGAATAAAGACAGTGAAGAGATTTATTTCAAAGCCTGGGTTTAAGCTGCTGAGAGAACACAAAAAACCCTAACAAAAATGGAATCGAAAGCAATTTGCTTAGGGTTTCTTCCTCCAAGACTTCGATTTTCATCTCCACGTTTACTCTCTCTTCCTCCTTCTCCTCCTGCTTCTTCCACATTTGCGACGCGTCACAAACTTGATTCCAGACAAACCCTCCTTTGGAACAAACCGCAATTGAGCCGAGTTCGTGTAGCGTGTTCTTCTTCTCAATCTGACTCAAGACCTGAGAAGAAGCAATCGGATAAGAGTAACTATGCTCGAGCTGAGCTGTTCCGTGGGAAATCAGGTTCTGTTTCTTTCAATGGTCTGACTCATCAGCTGGTTGAAGAAAGTAAACTGGTTTCAGCTCCGTTTCAAGAAGAGAAAGGTTCTTTCTTGTGGGTTTTGGCTCCTGTTGTTTTGATTTCTTCGTTGATTCTTCCTCAGTTCTTTCTAAGTGGTATCATTGAAGCTACCTTCAAAAACGACACTGTTGCTGAAATTGTTACTTCTTTTTGCTTTGAGACGGTGTTTTATGCTGGTCTTGCGATATTCCTGTCTGTGACTGACCGAGTGCAGAGGCCGTACTTAGACTTCAGCTCCAAGAGATGGGGTCTGATCACTGGACTGAGGGGATACCTTACGTCTGCATTCCTCACGATGGGTTTAAAAGTTGTAGTTCCCGTATTTGCTGTTTACATGACTTGGCCAGCTCTTGGAATAGATGCTTTGATTGCAGTGCTTCCTTTCTTGGTTGGCTGTGCAGTTCAAAGAGTTTTCGAGGCTCGGCTTGAAAGACGTGGCTCATCCTGTTGGCCCATTGTTCCAATAGTCTTTGAGGTGTATAGGCTGTATCAGGTGACAAGAGCAGCGACTTTTGTTCAGAGGCTGATGTTTATGATGAAAGATGCGGCAACGACTGCTGAAATAACAGAGCGAGGAGTTGCACTAGTTGGTTTGGTTGTGACTTTGCAGTTTCTAGCTGTTATGTGTCTCTGGTCGTTTATCACTTTTCTTATGCGCCTCTTTCCTTCTAGACCTGTAGGTGAAAACTACTAGATCTCAGTGTTTAGTGATTGTTAGATGTAGCCAAATCCCATCGGTTTTGTTTTGTTTCTGTGTTCATTTCAGTAGTAATGAATTGTATTAAGTCACTTTAAGAATTGGTTGATCATGTGAAATGAGAATTGGCTGGAAATGTTATAGAACG%0D%0A%3EAT5G50570.2%20%7C%20Symbols%3A%20%20%7C%20squamosa%20promoter-binding%20protein%2C%20putative%20%7C%20chr5%3A20599309-20601106%20REVERSE%0D%0AAAAAAGGACAAATCTTGATATTGCTTTGATTGCTGTTGTGTATGTATGTGTTTTTATAGTGAGAGAAGAAAAAAAAGCACAATCTTTGAATGGACTGGAATTTCAAACTTAGCTCTGGTTATTTATCTGGATTCGATCAAGAACCAGATTTATCACCAATGGATGGTTCGATCTCGTTTGGTGGGTCGTCACAGTCAAAAGCGGATTTTTCATTTGATCTAAAACTTGGAAGAAACATTGGAAACTCTTCCTCTGTTTTTGGTGATACAGAGCAAGTGATTAGTCTTAGTAAGTGGAAAGATAGTGCTTTAGCTAAACCAGAAGGTTCAAGAAGCTCGAGTTCAAAGAGAACAAGAGGGAATGGTGTTGGAACCAACCAGATGCCGATTTGTCTTGTTGATGGATGTGATTCTGATTTTAGTAATTGTAGAGAGTATCATAAGAGACATAAAGTTTGTGATGTTCATTCAAAAACTCCTGTGGTTACTATTAATGGTCATAAACAGAGGTTTTGTCAACAATGCAGCAGGTTTCATGCTTTGGAGGAGTTTGATGAAGGGAAGAGAAGTTGTAGGAAACGTCTTGATGGACATAATCGAAGACGACGGAAGCCGCAGCCTGAACATATCGGTCGTCCTGCCAACTTCTTTACGGGTTTTCAAGGTAGCAAATTGCTAGAGTTTTCTGGTGGTTCACATGTGTTTCCAACTACATCTGTGTTGAACCCGAGCTGGGGAAATAGTCTTGTAAGCGTTGCTGTAGCCGCCAATGGTTCGAGTTATGGGCAGAGCCAGAGCTATGTTGTTGGTTCTTCTCCTGCAAAGACAGGGATAATGTTTCCAATCTCTTCTTCTCCAAACAGTACCAGAAGCATAGCAAAACAATTCCCTTTCTTGCAAGAAGAAGAAAGCTCGAGAACCGCATCGTTGTGTGAGAGAATGACGAGTTGCATCCATGACTCTGATTGTGCTCTCTCTCTTCTGTCATCCTCCTCGTCGTCAGTCCCTCATTTGCTTCAACCACCACTTTCTTTGTCCCAAGAAGCAGTTGAGACAGTTTTTTACGGGTCGGGATTGTTTGAGAATGCGAGTGCAGTCTCTGATGGATCGGTTATATCCGGTAACGAGGCTGTCCGTCTTCCGCAGACATTCCCGTTTCATTGGGAGTAGTAGAAGAAGAAGTAGGTAGATAGATAGAATCAGAAAGATCTATTTGTGTCTCTTCTCTTCTCCCTCATTTTTCAATGTTCTTTATCATCATCATTGTTCTTGTTAACACTACAAGAAATATGGACATTCTTAACACACCGAAAACGCTATAATAACGTTTACATAGCGGATTCATAAACGCTGTGTTTGCCGGAGCTATCTTAGAGTGGTCACATACAATAGCGTTTCTTGCTATGCTATTAAATGTTCACATATTATGGCGTAAAAGAAATGCTTTGATTTCCTTTGTTATTGCACAATTTTGATGTTATACTTTTGTAACTCTTTTTAAGGGCTATAAACTATTATTTTGTAGCTATATTTTATAGGCTATGATCTGATATGTTGTAGCTTTATTTTTTTGGCTATAAATCTATAAACAGCCTATAAGTTTGGGACATTTTGTTACACATTTGCAAGAAACCTCTTTTTG%0D%0A&target_uploaded=&top=200'
headers = {
  'sec-ch-ua': '" Not A;Brand";v="99", "Chromium";v="98", "Google Chrome";v="98"',
  'sec-ch-ua-mobile': '?0',
  'sec-ch-ua-platform': '"Windows"',
  'Upgrade-Insecure-Requests': '1',
  'Content-Type': 'application/x-www-form-urlencoded',
  'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/98.0.4758.102 Safari/537.36',
  'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9',
  'Sec-Fetch-Site': 'same-origin',
  'Sec-Fetch-Mode': 'navigate',
  'Sec-Fetch-User': '?1',
  'Sec-Fetch-Dest': 'document',
  'Cookie': 'session=eyJzaW1wbGVwYWdlIjpmYWxzZX0.YhEQeA.N4sPemlWsAA1q7EwbfdkVXJbqs8'
}

response = requests.request("POST", url, headers=headers, data=payload)
type(response)


# In[ ]:




