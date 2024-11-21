import os
import glob
import argparse
from collections import defaultdict

def split_fasta_by_sample(input_dir, output_dir, fasta_format):
    # Get a list of all FASTA files in the input directory with the specified extension
    fasta_files = glob.glob(os.path.join(input_dir, f"*.{fasta_format}"))

    # Create the output directory if it does not exist
    os.makedirs(output_dir, exist_ok=True)

    for input_fasta in fasta_files:
        # Dictionary to hold sequences grouped by sample ID
        sample_dict = defaultdict(list)
        
        # Get the input file name without the extension
        input_name = os.path.splitext(os.path.basename(input_fasta))[0]

        # Read the input FASTA file
        with open(input_fasta, 'r') as file:
            current_header = None
            current_sequence = []

            for line in file:
                line = line.strip()

                if line.startswith(">"):
                    if current_header:
                        sample_id = current_header.split('S')[1].split('C')[0]
                        sample_dict[sample_id].append((current_header, ''.join(current_sequence)))

                    # Reset for the new header
                    current_header = line
                    current_sequence = []
                else:
                    # Sequence line
                    current_sequence.append(line)

            # Handle the last sequence in the file
            if current_header:
                sample_id = current_header.split('S')[1].split('C')[0]
                sample_dict[sample_id].append((current_header, ''.join(current_sequence)))

        # Write output FASTA files for each sample ID
        for sample_id, sequences in sample_dict.items():
            output_file = os.path.join(output_dir, f"S{sample_id}_{input_name}.fasta")
            with open(output_file, 'w') as out_file:
                for header, sequence in sequences:
                    out_file.write(f"{header}\n{sequence}\n")

        print(f"FASTA files split by sample ID from {input_fasta} successfully!")

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Split FASTA files by sample ID.')
    parser.add_argument('--input_dir', required=True, help='Input directory containing FASTA files')
    parser.add_argument('--output_dir', required=True, help='Output directory to save split FASTA files')
    parser.add_argument('--format', default='fasta', help='Extension of the input FASTA files (default: fasta)')

    args = parser.parse_args()

    # Run the function with the provided arguments
    split_fasta_by_sample(args.input_dir, args.output_dir, args.format)
