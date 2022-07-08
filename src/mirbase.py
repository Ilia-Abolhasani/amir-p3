mirbase_url = "https://www.mirbase.org/ftp/CURRENT"

def download(mirbase_dir):
    os.system(f'rm -r {mirbase_dir}')
    os.system(f'mkdir -p {mirbase_dir}')
    items = [
        'aliases.txt.gz',
        'hairpin.fa.gz', 
        'hairpin_high_conf.fa.gz',
        'mature.fa.gz',
        'mature_high_conf.fa.gz',
        'miRNA.str.gz',
        'miRNA.xls.gz',
        'organisms.txt.gz'
        ]
    for item in items:
        os.system(f'wget {mirbase_url}/organisms.txt.gz -P {mirbase_dir}/')
        os.system(f'gzip -d {mirbase_dir}/organisms.txt.gz')


def select_mirs(mirbase_dir):
    mature = fasta_to_df(f'{mirbase_dir}/mature.fa')
    mature_high_conf = fasta_to_df(f'{mirbase_dir}/mature_high_conf.fa')
    mature['trim tag'] = mature['tag'].apply(lambda line: ' '.join(line.split(' ')[:2]))
    mature['confidence'] = mature['trim tag'].isin(mature_high_conf['tag'])

    mature['organism'] = mature['tag'].apply(lambda x: x[:3])

    organism = pd.read_csv(f'{mirbase_dir}/organisms.txt',sep='\t')
    organism.columns = [c.replace('#','') for c in organism.columns] # remove sharp from columns
    
    items = list(organism['tree'].unique())
    items.sort(key=len)

    selectedTree = organism[organism['tree'].apply(lambda x: "Viridiplantae;" in x)]
    #selectedTree = selectedTree[selectedTree['name'] == ""]

    selected = mature[mature['organism'].isin(selectedTree['organism'])]
    return selected