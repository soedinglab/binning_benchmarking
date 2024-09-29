import os
import argparse
from collections import defaultdict
# Function to extract sample id from a read_id
def extract_read_sample_id(read_id):
    return read_id.split('R')[0].replace('S','')  # Sample ID is between 'S' and 'R'

# Function to extract sample id from a contig id
def extract_contig_sample_id(contig_id):
    return contig_id.split('C')[0].replace('S','')  # Sample ID is between 'S' and 'C'

# Read file A (read_id -> otuid)
def parse_readmapping(read_mapping):
    sample_otu = defaultdict(lambda: defaultdict(list))  # {sample_id: {otuid: [read_ids]}}
    
    with open(read_mapping, 'r') as rmap:
        for line in rmap:
            read_id, otuid, tax_id = line.strip().split()
            sample_id = extract_read_sample_id(read_id)
            sample_otu[sample_id][otuid].append((read_id, tax_id))
    
    return sample_otu

# Read file B (read_id -> contig_id)
def parse_contigmapping(contig_mapping):
    read_to_contig = defaultdict(list)
    
    with open(contig_mapping, 'r') as cmap:
        for line in cmap:
            read_id, contig_id = line.strip().split()
            read_to_contig[read_id].append(contig_id)
    
    return read_to_contig

# Main function to calculate the fraction and assign OTUs
def calculate_contig_otu_mapping(sample_otu, read_to_contig):
    # {sample_id: {contig_id: {otuid: count}}}
    sample_contig_otu_fraction = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    # {sample_id: {otuid: total_read_count}}
    otu_read_counts = defaultdict(lambda: defaultdict(int))
    # {sample_id: {otuid: tax_id}}
    otu_taxid_mapping = defaultdict(lambda: defaultdict(str))


    # Step 1: Count reads per OTU and contig, only if samples match
    for sample_id, otus in sample_otu.items():
        for otuid, reads in otus.items():
            # Total reads for this OTU in this sample
            otu_read_counts[sample_id][otuid] = len(reads)
            otu_taxid_mapping[sample_id][otuid] = reads[0][1]

            for read_id, _ in reads:
                if read_id in read_to_contig:
                    contig_ids = read_to_contig[read_id]
                    for contig_id in contig_ids:
                        contig_sample_id = extract_contig_sample_id(contig_id)
                        # Ensure both read and contig come from the same sample.
                        # If not found due to poor alignment. Use 1-to-1 sample-wise reads and contigs for mapping
                        if sample_id == contig_sample_id:  
                            sample_contig_otu_fraction[sample_id][contig_id][otuid] += 1
    
    return sample_contig_otu_fraction, otu_read_counts, otu_taxid_mapping

# Select the OTU with the highest fraction for each contig in each sample
def assign_otu_to_contig(sample_contig_otu_fraction, otu_read_counts, otu_taxid_mapping):
    contig_otu_assignments = defaultdict(dict)  # {sample_id: {contig_id: otuid}}
    
    for sample_id, contigs in sample_contig_otu_fraction.items():
        for contig_id, otus in contigs.items():
            best_otu = None
            best_fraction = 0.0
            best_tax_id = None

            for otuid, count in otus.items():
                total_reads_for_otu = otu_read_counts[sample_id][otuid]
                fraction = count / total_reads_for_otu  # Fraction of reads from OTU assigned to this contig

                if fraction > best_fraction:
                    best_fraction = fraction
                    best_otu = otuid
                    best_tax_id = otu_taxid_mapping[sample_id][otuid]

            contig_otu_assignments[sample_id][contig_id] = (best_otu, best_tax_id)
    
    return contig_otu_assignments

def parse_contig_lengths(contig_length):
    contig_lengths = {}
    
    with open(contig_length, 'r') as file_c:
        for line in file_c:
            contig_id, length = line.strip().split()
            contig_lengths[contig_id] = int(length)  # Store lengths as integers
    
    return contig_lengths

# Main function to run the process
def main(read_mapping, contig_mapping, contig_length, outputdir):
    # Step 1: Parse the files
    sample_otu = parse_readmapping(read_mapping)
    read_to_contig = parse_contigmapping(contig_mapping)
    contig_lengths = parse_contig_lengths(contig_length)
    sample_contig_otu_fraction, otu_read_counts, otu_taxid_mapping = calculate_contig_otu_mapping(sample_otu, read_to_contig)
    contig_otu_assignments = assign_otu_to_contig(sample_contig_otu_fraction, otu_read_counts, otu_taxid_mapping)
    
    sample_id = list(contig_otu_assignments.keys())[0]
    with open(os.path.join(outputdir, f'{sample_id}_contigs_otuidsmapping.tsv'), 'w+') as f:
        for sample_id, contigs in contig_otu_assignments.items():
            for contig_id, (otuid, tax_id) in contigs.items():
                length = contig_lengths.get(contig_id)
                if length == None:
                    raise RuntimeError(f'{contig_id} length is {length}. Check input!')
                f.write(f"{contig_id}\t{otuid}\t{tax_id}\t{length}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="gsmapping",
        description="get contigs to otu pairs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage="%(prog)s --readmapping --contigmapping --outdir",
        add_help=True,
    )

    parser.add_argument("--readmapping", type=str, \
        help="a tsv file with read and otuid mapping", required=True)
    parser.add_argument("--contigmapping", type=str, \
        help="a tsv file with contig and read id mapping", required=True)
    parser.add_argument("--length", type=str, \
        help="a tsv file with contigs and length", required=True)
    parser.add_argument("--outdir", type=str, \
        help="output directory", required=True)

    args = parser.parse_args()

    args.outdir = os.path.abspath(args.outdir)
    main(args.readmapping, args.contigmapping, args.length, args.outdir)