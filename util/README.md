# aligner2counts
This repository contains script to process read mapping file (.bam) on the fly with mapping process. It produces the following outputs,
 - read abundance of contigs (fractional count of total reads mapped to contig, `<sample_id>_count`)
 - fraction of reads shared between contigs if found (split-read mapping, `<sample_id>_countlinks`)
 - read abundance by uniquely mapped reads (`<sample_id>_uniqcount`)
 - read abundance by reads that mapped to multiple contigs (`<sample_id>_crosscount`)
 - read coverage ((read abundance * read length) / contig length, `<sample_id>_coverage`)
 - alignment file in `.sam` format (`<sample_id>.sam`)
 - read-contig pair from alignment (`<sample_id>_mapids`)
As of now, it is specialized to process only bowtie mapping to compute abundance and coverage. If you have alignment from different tool and want only to get read-contig pair, use `--only-mapids` option

# Installation
	g++ -o aligner2counts aligner2counts.cpp -O3
	export PATH=${PATH}/$(pwd) (save this line in .bashrc or .bash_profile for easy access)
Require gcc version `>=9.4.0`


# Usage
	Usage: aligner command | aligner2counts outdir outputname [--minlength N] [--no-coverage] [--single] [--only-mapids] [--qcov X] [--seq-id Y]
	Options:
  	<outdir>            Directory where output files will be stored.
  	<outputname>        Base name for the output files (without extensions).

 	 --minlength N      Minimum length of contigs to be considered. (Default: N=1000)
  	                    Alignments for contigs shorter than this length will be ignored.

  	--no-coverage       Do not output coverage. (Optional)
  	                    This flag disables the output of contig coverage.
	
  	--single            Process input as single-end reads. (Optional)
  	                    By default, paired-end reads are expected unless this flag is set.

  	--only-mapids       Output only the mapping identifiers. (Optional)
  	                    This flag restricts the output to just the IDs of reads and mapped contigs.

  	--qcov X            Minimum query/read coverage threshold X. (Optional, default=99.0%)
  	                    Specifies the minimum percentage of query sequence that must be aligned.

  	--seq-id Y          Minimum sequence identity percentage Y. (Optional, default=97.0%)
  	                    Filters alignments by requiring at least Y% identity.

  	-h, --help          Display this help message and exit.

# Example
One go mapping and processing,

 	bowtie2 -q -x <ref_index> <reads.fq> | aligner2counts <output_directory> <sample_id>

Process only `.bam` or `.sam` (require samtools preinstalled)

 	samtools view -h <input_alignment> | aligner2counts <output_directory> <sample_id>

If you have `.sam` file, you can also use `cat <input>.sam | aligner2counts <output_directory> <sample_id>`
