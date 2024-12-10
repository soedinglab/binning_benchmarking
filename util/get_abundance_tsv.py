#!/usr/bin/env python

import pandas as pd
import argparse
import glob
import os
import sys
import re

parser = argparse.ArgumentParser(prog='get_abundance_tsv.py', description="Process abundance files\n")
parser.add_argument("-i", "--input_dir", type=str, help="directory of abundance files", required=True)
parser.add_argument("-l", "--contig_length", type=str, help="Contig length", required=True)
parser.add_argument("-m", "--minlength", type=int, help="Minimum length of contigs", default=1000)
args = parser.parse_args()

directory = args.input_dir
length_file = args.contig_length
minlength = args.minlength

if not os.path.isdir(directory):
    print(f"Error: The directory '{directory}' does not exist or is not a directory.")
    sys.exit(1)

if not os.path.isfile(length_file):
    print(f"Error: The contig length file '{args.contig_length}' does not exist or is not a file.")
    sys.exit(1)

abund_files = glob.glob(os.path.join(directory, 'abundances*.tsv'))
data_list = []
for file_path in abund_files:
    df = pd.read_csv(file_path, sep='\t', header=None)
    
    file_name = os.path.basename(file_path)
    file_id = os.path.splitext(file_name)[0].replace('abundances_','')
    
    df['2'] = file_id
    data_list.append(df)

fractional_counts = pd.concat(data_list, ignore_index=True)
contig_length = pd.read_csv(length_file,header=None,sep='\t')
length_dict = dict(zip(contig_length[0], contig_length[1]))
fractional_counts.columns = ['contig_id', 'abundance', 'sample_id']
read_counts = fractional_counts.pivot_table(index = 'contig_id', columns = 'sample_id', values = 'abundance').sort_index(axis=1)
sample_count = len(read_counts.columns)
read_counts['contigName'] = read_counts.index
read_counts['contigLen'] = read_counts['contigName'].map(length_dict)
read_counts = read_counts[read_counts['contigLen']>=minlength]
read_counts.loc[:,'totalAvgDepth'] = read_counts.iloc[:,0:sample_count].sum(axis=1)
read_counts = read_counts[list(read_counts.columns[-3:])+list(read_counts.columns[:-3])]
num_columns = len(read_counts.columns)
var_columns = pd.DataFrame(0.0, index=read_counts.index, columns=[f'{i}_var' for i in range(num_columns - 3)])
cols = []
for i in range(3, num_columns):
    cols.append(read_counts.columns[i])
    cols.append(f'{i-3}_var')

final_columns = list(read_counts.columns[:3]) + cols
read_counts = pd.concat([read_counts, var_columns], axis=1)
read_counts = read_counts[final_columns]

def natural_sort_key(s):
    return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', s)]

# sort rows by contig ids
read_counts = read_counts.sort_values(by='contigName', key=lambda col: col.map(natural_sort_key))
read_counts.to_csv(os.path.join(directory, 'abundances_gf.tsv'), sep='\t', index=False)


