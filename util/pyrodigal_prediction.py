#!/usr/bin/env python
import os
import sys
import time
import pyhmmer
import pathlib
import argparse
import pyrodigal
import collections
import numpy as np
from Bio import SeqIO
import multiprocessing
import pandas as pd

parent_path = os.path.dirname(__file__)
sys.path.insert(0, parent_path)
# This code is referred from markersgf script of genomeface and https://pyhmmer.readthedocs.io/en/stable/examples/fetchmgs.html#Getting-the-cutoffs


# reference SCMG:
# 1) 40 genes: Creevey et al (2011), https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0022099#pone-0022099-t001
# 2) 120 bacterial and 122 archea genes: Parks et al (2017), https://www.nature.com/articles/s41564-017-0012-7


# BitScore cutoff obtained from https://github.com/motu-tool/fetchMGs/tree/master/lib
cutoff = {'COG0012': 210.0, 'COG0016': 240.0, 'COG0018': 340.0, \
    'COG0048': 100.0, 'COG0049': 120.0, 'COG0052': 140.0, \
    'COG0080': 90.0, 'COG0081': 130.0, 'COG0085': 1020.0, \
    'COG0086': 60.0, 'COG0087': 120.0, 'COG0088': 110.0, \
    'COG0090': 180.0, 'COG0091': 80.0, 'COG0092': 120.0, \
    'COG0093': 80.0, 'COG0094': 110.0, 'COG0096': 80.0, \
    'COG0097': 100.0, 'COG0098': 140.0, 'COG0099': 120.0, \
    'COG0100': 80.0, 'COG0102': 100.0, 'COG0103': 80.0, \
    'COG0124': 320.0, 'COG0172': 170.0, 'COG0184': 60.0, \
    'COG0185': 70.0, 'COG0186': 80.0, 'COG0197': 70.0, \
    'COG0200': 60.0, 'COG0201': 210.0, 'COG0202': 80.0, \
    'COG0215': 400.0, 'COG0256': 70.0, 'COG0495': 450.0, \
    'COG0522': 80.0, 'COG0525': 740.0, 'COG0533': 300.0, \
    'COG0541': 450.0, 'COG0552': 220.0}

gene_finder = pyrodigal.GeneFinder(meta=True)
alphabet = pyhmmer.easel.Alphabet.amino()
hmms = []
Result = collections.namedtuple("Result", ["name","query", "cog", "bitscore"])
def predict_genes(entries):
    """ predict genes """
    annotations = []
    proteins = []
    for entry in entries:
        genes = gene_finder.find_genes(str(entry.seq).encode())
        for idx, gene in enumerate(genes):
            protein = pyhmmer.easel.TextSequence(name=(str(idx) + "_" + str(entry.id)).encode(),
                sequence=gene.translate()).digitize(alphabet)
            proteins.append(protein)

    hits = pyhmmer.hmmsearch(hmms, proteins, cpus=1)
    for hit in hits:
        for domain_hit in hit:
            if domain_hit.included and not domain_hit.duplicate:
                p_nm = domain_hit.name.decode()
                record_id = p_nm[p_nm.index('_')+1:]
                if domain_hit.score >= cutoff[hit.query_name.decode()]:
                    annotations.append(\
                    Result(p_nm, record_id, hit.query_name.decode(), domain_hit.score))

    return annotations

def load_hmms(hmms_path) -> list:
    """ Load the HMMs """
    files = list(pathlib.Path(hmms_path).glob('*.hmm'))
    if len(files) == 0:
        raise RuntimeError(f'hmm files are missing')

    for hmm_file in files:
        with open(hmm_file, "rb") as hmm_f:
            with pyhmmer.plan7.HMMFile(hmm_f) as hmm_file:
                hmm = hmm_file.read()
                hmms.append(hmm)
    return hmms



def filter_results(all_results):
    """ filter results by bit score """
    best_results = {}
    keep_query = set()
    # TODO: when multiple SCMGs are mapped to contig, we remove it.
    # But can we make use of this information in some way?
    for result in all_results:
        query = result.query
        if query in best_results:
            previous_bitscore = best_results[query].bitscore
            if result.bitscore > previous_bitscore:
                best_results[query] = result
                keep_query.add(query)
            elif result.bitscore == previous_bitscore:
                if best_results[query].cog != result.cog:
                    if query in keep_query:
                        keep_query.remove(query)
        else:
            best_results[query] = result
            keep_query.add(query)

    # for k in best_results:
    #     print(k, 'k', flush=True)

    results = [best_results[k] for k in best_results if k in keep_query]

    return results

def gene_prediction(seq, outdir):
    s=time.time()
    record = SeqIO.parse(seq, "fasta")
    # contig_names = np.load(names, allow_pickle=True)['arr_0']
    # contig_names = pd.read_csv(names, header=None).to_numpy().flatten()

    records = []
    for f in record:
        records.append(f)
    del record

    num_seqs = len(records)
    CHUNK_SIZE_FACTOR = 1
    pool_size = 1
    chunk_size = num_seqs // (pool_size * CHUNK_SIZE_FACTOR)
    record_chunks = [records[i:i+chunk_size] for i in range(0, len(records), chunk_size)]

    hmms = load_hmms(parent_path + '/markerhmms')

    with multiprocessing.Pool(pool_size) as pool:
        predictions = pool.map(predict_genes, record_chunks)
    del record_chunks

    # remove empty predictions
    all_results = [result for results in predictions for result in results]

    filtered_results = filter_results(all_results)

    # with open(outdir + '/marker_hits', 'w', encoding='utf-8') as file:
    #     for hit in filtered_results:
    #         index_array = np.where(contig_names == hit.query)[0]
    #         if index_array.size > 0:
    #             index = index_array[0]
    #             file.write(str(index) +'\t' + hit.cog +'\n')
    # file.close()

    with open(outdir + '/marker_hits', 'w', encoding='utf-8') as file:
        for hit in filtered_results:
            file.write(str(hit.query) +'\t' + hit.cog +'\n')
    file.close()
    print(time.time()-s, 'seconds for SCMGs prediction')

    # s2 = time.time()
    # data = [[int(np.where(contig_names==hit.query)[0]), hit.cog] for hit in filtered_results]
    # with open(outdir + "marker_hits_withids.csv", 'w', newline='', encoding='utf-8') as file:
    #     writer = csv.writer(file, delimiter='\t')
    #     writer.writerows(data)



if __name__ == "__main__":

    start = time.time()
    parser = argparse.ArgumentParser(
        prog="pyrodigal_prediction",
        description="Predict scmgs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage="%(prog)s --seq --outdir",
        add_help=True,
    )

    parser.add_argument("--seq", type=str, \
        help="contig sequence in fasta", required=True)
    # parser.add_argument("--names", type=str, \
    #     help="ids of contigs", required=True)
    parser.add_argument("--outdir", type=str, \
        help="output directory", required=True)

    args = parser.parse_args()

    # if not args.seq.endswith(".fasta") or not args.seq.endswith(".fa"):
    #     raise IOError(f"{args.seq} is not fasta file")
    gene_prediction(args.seq, args.outdir)
    
