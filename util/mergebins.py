#!/usr/bin/env python
import pandas as pd
import subprocess
import sys, os
import re
from collections import Counter
from itertools import combinations
from Bio import SeqIO
import numpy as np
import shutil
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster
from concurrent.futures import ThreadPoolExecutor
import warnings
import argparse

warnings.filterwarnings("ignore", message="UserWarning: The value of the smallest subnormal for <class 'numpy.float32'> type is zero.")

current_path = os.environ.get('PATH')
pyrodigal_path = os.path.dirname(os.path.realpath(__file__))
def check_arguments(args):

    # check bin directory
    if not os.path.isdir(args.bindir):
        raise FileNotFoundError(f"The directory '{args.bindir}' does not exist.")

    # check input file format
    input_binsfastalist = [(os.path.join(args.bindir, f), f) for f in os.listdir(args.bindir) if f.endswith('.'+args.format)]
    if not input_binsfastalist:
        raise FileNotFoundError(f"No files with extension '{args.format}' found in the directory '{args.bindir}'.")

    # check if plass is accessible
    plass_executable = os.path.join(args.plassdir, 'plass')
    if not os.path.isfile(plass_executable) or not os.access(plass_executable, os.X_OK):
        raise FileNotFoundError(f"The executable '{plass_executable}' is not accessible or not executable.")
    result = subprocess.run([plass_executable], shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Failed to execute '{plass_executable}'. {result.stderr}")

def select_contigwithsmgs(input_binfasta, idslist, output_fasta):
    ids_to_extract_set = set(idslist)

    with open(input_binfasta, "r") as infile, open(output_fasta, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id in ids_to_extract_set:
                SeqIO.write(record, outfile, "fasta")

def split_binbysample(input_binfasta, subbindir):
    with open(input_binfasta, 'r') as file:
        fasta_content = file.read()
    # extract sample id
    header_pattern = re.compile(r'(>.*S(\d+)C[^\n]*\n)')
    headers = [(m.start(), m.group(1), m.group(2)) for m in re.finditer(header_pattern, fasta_content)]

    sequences_by_sample_id = {}

    for i in range(len(headers)):
        start_pos = headers[i][0]
        sample_id = headers[i][2]
        end_pos = headers[i + 1][0] if i + 1 < len(headers) else len(fasta_content)
        sequence_block = fasta_content[start_pos:end_pos]
        if sample_id not in sequences_by_sample_id:
            sequences_by_sample_id[sample_id] = []
        sequences_by_sample_id[sample_id].append(sequence_block)
    
    for sample_id, sequences in sequences_by_sample_id.items():
        output_fasta = f'{sample_id}_sample.fasta'
        with open(os.path.join(subbindir, output_fasta), 'w') as out_file:
            out_file.writelines(sequences)

    return np.array(list(sequences_by_sample_id.keys())).astype(int)

def clearpath(path):
    try:
        assembly_files = os.listdir(path)

        for f in assembly_files:
            if f.endswith("assembly.fasta") or f.endswith("sample.fasta") or f.endswith("marker.fasta"):
                os.remove(os.path.join(path, f))
    except:
        pass

def get_markersstat(query_sequence, subbindir):
    subprocess.run([f'python {pyrodigal_path}/pyrodigal_prediction.py --seq {query_sequence} --outdir {subbindir}'],
        shell=True,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE)

    markerfile = os.path.join(subbindir,'marker_hits')
    if not os.stat(markerfile).st_size == 0:
        markerhits = pd.read_csv(subbindir+'marker_hits', header=None, sep='\t')
        cogs_counts = Counter(markerhits[1])
        total_cogs = len(cogs_counts)
        unique_cogs = sum(1 for count in cogs_counts.values() if count == 1)
        duplicate_cogs = sum(1 for count in cogs_counts.values() if count > 1)
        cont = duplicate_cogs / total_cogs
        os.remove(os.path.join(subbindir,'marker_hits'))
        return unique_cogs, duplicate_cogs, cont, list(markerhits[0])
    else:
        os.remove(os.path.join(subbindir,'marker_hits'))
        return 0, 0, 0, []

def segregate_subbins(sampleids, bindir, subbindir, outputname, foldertype):
    
    if foldertype == 'merge':
        filelist = [f'{subbindir}/{sid}_sample.fasta' for sid in sampleids]
        file_list_str = ' '.join(filelist)
        outputdir = os.path.join(*[bindir,'merge/', outputname])
        print(outputdir, 'outputdir merge')
        subprocess.run([f"cat {file_list_str} > {outputdir}"], shell=True)
    else:
        for sid in sampleids:
            inputfile = f'{subbindir}/{sid}_sample.fasta'
            outputfile = os.path.join(*[bindir,'nonmerge/',f'{outputname}_{sid}.fasta'])
            print(outputfile, 'outputfile')
            subprocess.run([f'mv {inputfile} {outputfile}'], shell=True)

def clustersubbins(samplewithscmgs, subbin_markerstat, subbindir):
    samplepairs = list(combinations(samplewithscmgs,2))
    sample_indexdict = {key : index for index, key in enumerate(subbin_markerstat.keys())}
    distance_matrix = np.full((len(subbin_markerstat.keys()), len(subbin_markerstat.keys())), 1e6)
    np.fill_diagonal(distance_matrix, 0)
    print(subbindir, len(samplepairs),'subbin directory and length of pairs', flush=True)
    for pair in samplepairs:
        pair1, pair2 = pair
        marker_file1 =  os.path.join(subbindir, f'{pair1}_marker.fasta')
        marker_file2 =  os.path.join(subbindir, f'{pair2}_marker.fasta')
        input_assembly = os.path.join(subbindir, f'{pair1}_{pair2}_marker.fasta')
        output_assembly = os.path.join(subbindir, f'{pair1}_{pair2}_assembly.fasta')
        subprocess.run(['cat ' + marker_file1 + ' ' + marker_file2 + ' > ' + input_assembly], shell=True)
        subprocess.run(['penguin nuclassemble ' + input_assembly + ' ' + output_assembly + ' tmp --num-iterations 1 --min-seq-id 0.995 --contig-output-mode 0 --max-seq-len 4000000'], \
            shell=True)#,
            # stdout = subprocess.PIPE,
            # stderr = subprocess.PIPE)
        try:
            shutil.rmtree(os.path.join(subbindir, 'tmp'))
        except:
            pass

        if not os.stat(output_assembly).st_size == 0:
            unique_cogs_ij, duplicate_cogs_ij, cont_ij, _ = get_markersstat(output_assembly, subbindir)
            merged_cogs = subbin_markerstat[pair1][0] + subbin_markerstat[pair2][0] - unique_cogs_ij
            delta_duplicate_cogs = duplicate_cogs_ij - subbin_markerstat[pair1][1] - subbin_markerstat[pair2][1]
            
            print(delta_duplicate_cogs, merged_cogs,
                unique_cogs_ij, duplicate_cogs_ij, cont_ij,
                subbin_markerstat[pair1][0], subbin_markerstat[pair1][1], subbin_markerstat[pair1][2],
                subbin_markerstat[pair2][0], subbin_markerstat[pair2][0], subbin_markerstat[pair2][2], flush=True)
            # NOTE: condition is strict.
            if delta_duplicate_cogs < 1 and merged_cogs >= 1: # cont_ij < 0.05 and
                distance_matrix[sample_indexdict[pair1], sample_indexdict[pair2]] = 1
                distance_matrix[sample_indexdict[pair2], sample_indexdict[pair1]] = 1
    print(distance_matrix, 'distance_matrix', flush=True)
    condensed_distance_matrix = squareform(distance_matrix)
    # Linkage clustering
    Z = linkage(condensed_distance_matrix, method='single')
    clusters = fcluster(Z, t=1, criterion='distance')
    clustersets = [samplewithscmgs[np.where(clusters==f)[0]] for f in set(clusters)]
    return clustersets
import threading
def parallelprocess_bins(binfasta, filename, bindir, format):
    binname = filename.replace('.'+format,'')
    print(binname, 'binname in thread:', threading.current_thread().name, flush=True)
    subbindir = os.path.join(bindir, filename.replace('.'+format,'/'))
    os.makedirs(subbindir, exist_ok=True)
    clearpath(subbindir)
    samplesinbin = split_binbysample(binfasta, subbindir)
    subbin_markerstat = {}
    for sample_id in samplesinbin:
        query_sequence = os.path.join(subbindir, f'{sample_id}_sample.fasta')
        unique_cogs_i, duplicate_cogs_i, cont, idslist = get_markersstat(query_sequence, subbindir)
        if idslist:
            subbin_markerstat[sample_id] = [unique_cogs_i, duplicate_cogs_i, cont]
            samplemarker_sequence = os.path.join(subbindir, f'{sample_id}_marker.fasta')
            select_contigwithsmgs(query_sequence, idslist, samplemarker_sequence)
    print(binname, subbin_markerstat, 'bin name and subbin markerstats', flush=True)
    if len(subbin_markerstat) > 1:
        samplewithscmgs = np.array(list(subbin_markerstat.keys())).astype(int)
        non_mergelist = np.setxor1d(samplesinbin, samplewithscmgs)        
        
        clustersets = clustersubbins(samplewithscmgs, subbin_markerstat, subbindir)
        
        for ids in clustersets:
            if len(ids) == 1:
                non_mergelist = np.append(non_mergelist,ids[0])
        clustersets = [cluster for cluster in clustersets if len(cluster)>1]

        if clustersets:
            mergecounter = 0
            for cluster in clustersets:
                segregate_subbins(cluster, bindir, subbindir, f'merge_{binname}_{mergecounter}.fasta', 'merge')
                mergecounter += 1
        if non_mergelist.size != 0:
            segregate_subbins(non_mergelist, bindir, subbindir, binname, 'nonmerge')
            clearpath(subbindir)
    else:
        segregate_subbins(samplesinbin, bindir, subbindir, binname, 'nonmerge')
        clearpath(subbindir)


def createsecondarypath(bindir):
    # create folder for merge bins
    paths = ['merge', 'nonmerge']
    # create folder for merge/non-merge bins
    for p in paths:
        try:
            secondarypath = os.path.join(bindir, p)
            if os.path.exists(secondarypath):
                shutil.rmtree(secondarypath)
            os.makedirs(secondarypath)
        except:
            pass

def mergebins(bindir,
    plassdir: str,
    threads:int,
    format: str = 'fasta',
    ):
    new_path = plassdir + os.pathsep + current_path
    os.environ['PATH'] = new_path

    createsecondarypath(bindir)

    # get extension from user
    input_binsfastalist = [(os.path.join(bindir, f), f) for f in os.listdir(bindir) if f.endswith('.'+format)]
    print(input_binsfastalist, 'bin fasta list')
    max_workers = min(threads, len(input_binsfastalist))
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        
        futures = [executor.submit(parallelprocess_bins, binfasta, filename, bindir, format) 
            for binfasta, filename in input_binsfastalist]
        for future in futures:
            future.result()

    try:
        shutil.rmtree(os.path.join(bindir, 'tmp'))
        # remove folders for each bin
        for _, filename in input_binsfastalist:
            binfoldername = filename.replace('.'+format,'/')
            subbindir = os.path.join(bindir, binfoldername)
            shutil.rmtree(subbindir)
    except:
        pass

# def fetchreads(bindir):
#     binlist = [f.replace('.'+format,'/') for f in os.listdir(bindir) if f.endswith('.'+format)]
#     # path = os.path.join(bindir, filename.replace('.'+format,'/'))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="mergebins",
        description="Merge redundant bins",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage="%(prog)s --bindir --format [optional] --plassdir [optional] --thread [optional]",
        add_help=True,
    )

    parser.add_argument("--bindir", type=str, \
        help="directory in which fasta files of bins located", required=True)
    parser.add_argument("--format", type=str, \
        help="bin file extension (fa/fas/fna/fasta)", default='fasta')
    parser.add_argument("--plassdir", type=str, \
        help="directory for plass executable", default="")
    parser.add_argument("--thread", type=int,\
        help="number of threads to be used", default=os.cpu_count())

    args = parser.parse_args()

    check_arguments(args)
    args.bindir = os.path.abspath(args.bindir)
    mergebins(args.bindir, args.plassdir, args.thread, args.format)