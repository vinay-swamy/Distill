#%%
import pandas as pd 
import numpy as np 
import os 
import pyfasta
import pybedtools
import concurrent.futures as cf
import argparse
from pyarrow.feather import read_feather
os.chdir('/data/swamyvs/Distill')
#%%

def encode_and_pad_seq(seq):
    res = [None]*len(seq)
    '''one hot encod nucleotide array. Any nonstandard nucleotide will be set to n
        
    '''
    for i, nt in enumerate(seq.lower()):
        if nt == 'a':
            res[i]=np.array([1,0,0,0,0])
        elif nt == 'c':
            res[i] = np.array([0,1,0,0,0])
        elif nt == 'g':
            res[i] = np.array([0,0,1,0,0])# leave space for n so array is in quasi alphabetical order 
        elif nt == 't':
            res[i] = np.array([0,0,0,0,1])
        else:# nt is notstandard 
            res[i] = np.array([0,0,0,1,0])
    res = np.array(res)
    return res
    # pad_size = max_seq_size- len(res) 
    # pad = np.zeros((pad_size, res.shape[1]))
    # return np.concatenate([res, pad], axis = 0)
def vectorize_variant(loc, pad):
    window_size = pad *2 +1
    chrom = loc[0]
    start = loc[1]
    end = loc[2]
    sub_type = loc[3]
    ref = loc[4]
    alt = loc[5]
    status = loc[6]
    if len(alt) > window_size:
        print(f'Insertion is longer than window_size {window_size} ; skipping\n')
        return None
    try:

        if sub_type == 'snp':
            #simple just swap out the middle nt with the alt
            alt_seq=genome[chrom][(start-pad):start] + alt + genome[chrom][end:(end+pad)]
        elif sub_type == 'ins':
            #
            new_pad = window_size - len(alt)
            lhs = int(new_pad/2)
            rhs = new_pad-lhs
            alt_seq=genome[chrom][(start-lhs):start] + alt + genome[chrom][end:(end+rhs)]
        elif sub_type == 'del':
            new_end =  start + len(alt)# technically this is the same as end, but i like to make it a little more explicit
            alt_seq = genome[chrom][(start-pad):(start+1)] + genome[chrom][new_end:(new_end+pad)]

        seq_vctr = encode_and_pad_seq(alt_seq)
        assert seq_vctr.shape == (window_size, 5)
        return((seq_vctr,status))
    except:
        print(f'{chrom}:{start}-{end} {sub_type} failed')
        return None


# %%
parser = argparse.ArgumentParser()
parser.add_argument('--pad', action='store', type = int)
parser.add_argument('--ncores', action = 'store', type = int)
parser.add_argument('--outpfx', action = 'store', type = str)
args = parser.parse_args()
pad = args.pad
NCORES = args.ncores
outpfx= args.outpfx

# %%
# pad = 5
# NCORES = 16
# outpfx = 'testing/all_vars'
#%%
all_variants = read_feather('feather_data/clean/all_variants.ftr')
genome = pyfasta.Fasta('ref/gencode_genome_grch37_clean.fa')
# convert numeric chroms to chr chroms to match gencode genome 
all_variants['chrom'] = all_variants.chrom.astype(str)
all_variants['start'] = all_variants.start.astype('int64')
all_variants['end'] = all_variants.end.astype('int64')
chr_chroms = [str(i) for i in range(1,24 )] + ['X' , 'Y' ]
is_numeric_chr = all_variants.chrom.isin( chr_chroms ) 
repl_chr =['chr'+i for i in  all_variants['chrom'][is_numeric_chr]]
all_variants['chrom'][is_numeric_chr] = repl_chr
#replace MT with M
all_variants['chrom'][all_variants['chrom'] == 'MT'] = 'chrM'
# convert ts/tp types t snp
is_snp= ~all_variants['sub_type'].isin(['ins', 'del']) 
all_variants['sub_type'][is_snp] = 'snp'
variant_list = all_variants[['chrom', 'start', 'end', 'sub_type', 'ref', 'alt', 'Status']].to_numpy()
#%%
with cf.ProcessPoolExecutor(max_workers=NCORES) as proc_exec:
    futures = {proc_exec.submit(vectorize_variant, var, pad) for var in variant_list }
    variant_vectors_raw = [ task.result() for task in cf.as_completed(futures)]

variant_vectors = np.asarray([tup[0] for tup in variant_vectors_raw if tup is not None ])
target = np.asarray([tup[1] for tup in variant_vectors_raw if tup is not None ])
#%%
np.save(f'{outpfx}_pad_{pad}_vectors.npy', variant_vectors)
np.save(f'{outpfx}_pad_{pad}_labels.npy', target)


# %%
