## Build CPP executables
	g++ -o aligner2counts aligner2counts.cpp -O3
	export PATH=${PATH}/$(pwd) (save this line in .bashrc or .bash_profile for easy access)
 	g++ -o extractreads extractreads.cpp -O3
	export PATH=${PATH}/$(pwd)
 	g++ -o convertfasta_multi2single convertfasta_multi2single.cpp -O3
	export PATH=${PATH}/$(pwd)
 	g++ -o splitreadsbysample splitreadsbysample.cpp -O3
	export PATH=${PATH}/$(pwd)
Require gcc version `>=9.4.0`

❗️If you compile from source under macOS we recommend installing and using `gcc` instead of `clang` as a compiler. `gcc` can be installed with Homebrew. Force `cmake` to use `gcc` as a compiler by running:

    CC="$(brew --prefix)/bin/gcc-14"
    CCX="$(brew --prefix)/bin/g++-14"

## aligner2counts
This repository contains script to process read mapping file (.bam) on the fly with mapping process. It produces the following outputs,
 - read abundance of contigs (fractional count of total reads mapped to contig, `<sample_id>_count`)
 - fraction of reads shared between contigs if found (split-read mapping, `<sample_id>_countlinks`)
 - read abundance by uniquely mapped reads (`<sample_id>_uniqcount`)
 - read abundance by reads that mapped to multiple contigs (`<sample_id>_crosscount`)
 - read coverage ((read abundance * read length) / contig length, `<sample_id>_coverage`)
 - alignment file in `.sam` format (`<sample_id>.sam`)
 - read-contig pair from alignment (`<sample_id>_mapids`)
As of now, it is specialized to process only bowtie mapping to compute abundance and coverage. If you have alignment from different tool and want only to get read-contig pair, use `--only-mapids` option

### Usage
	Usage: aligner command | ./aligner2counts outdir outputname [--minlength N] [--no-coverage] [--single] [--only-mapids] [--qcov X] [--seq-id Y] [--strobealign]
	Options:
	<outdir>            Directory where output files will be stored.
	<outputname>        Base name for the output files (without extensions).

	--minlength N       Minimum length of contigs to be considered. (Default: N=1000)
			    Alignments for contigs shorter than this length will be ignored.

	--no-coverage       Do not output coverage. (Optional)
			    This flag disables the output of contig coverage.

	--single            Process input as single-end reads. (Optional)
			    By default, paired-end reads are expected unless this flag is set.

	--strobealign       Alignment is generated from stobealign with `--eqx` flag. (Optional, default=bowtie2)

	--only-mapids       Output only the mapped pairs of read and contig identifiers. (Optional)
			    This flag restricts the output to just the IDs of reads and mapped contigs.

	--qcov X            Minimum query/read coverage threshold X. (Optional, default=99.0%)
			    Specifies the minimum percentage of query sequence that must be aligned.

	--seq-id Y          Minimum sequence identity percentage Y. (Optional, default=97.0%)
			    Filters alignments by requiring at least Y% identity.

	-h, --help          Display this help message and exit.

### Example
On the fly mapping and processing,

 	bowtie2 -q -x <ref_index> <reads.fq> | aligner2counts <output_directory> <sample_id>

Process only `.bam` or `.sam` (require samtools preinstalled)

 	samtools view -h <input_alignment> | aligner2counts <output_directory> <sample_id>

If you have `.sam` file, you can also use `cat <input>.sam | aligner2counts <output_directory> <sample_id>`

## Generate abundance matrix
`python get_abundance_tsv.py -i <inputdir> -l <contig_length> -m <minlength|1000>`

`inputdir` - a directory where sample-wise abundance text files are located

`contig_length` - a `tab` separated text file with contig_id and length

## Accuracy prediction of embedding space
`python linearclassifier_embeddingspace.py --latent <abundance.npy> --otuids <contig_otuids> --names <contig_names> --outdir <outputdir|pwd>`

`abundance.npy` - abundance matrix in `npy` format, generated from previous command using `get_abundance_tsv.py`.

`contig_otuids` - a text file with otuids of contigs with the same order as in abundance matrix

`contig_names` - a text file with contig's names with the same order as in abundance matrix

## Predict single-copy marker genes
`python pyrodigal_prediction.py --seq <contigseq.fasta> --outdir <outputdir>`

Required `pyhmmer` and `pyrodigal` installed.

## Extract reads from raw reads file to perform reassembly
`extractreads <fullpath/binfastafolder> <allsample_mapfile> <all_reads.corr.reads.fq> -f (binformat|fasta)`

## Filter scaffolds by length
`covertfasta_multi2single <contigseq.fasta> <outputdir> <minlength|1000> --length (print length, optional)`

To filter scaffolds obtained from reassembly, we used `minlength`=500 for this study. This script can also be used to get `contig_length` file.

## Mapping OTU IDs to assembled contigs based on the highest fraction of reads mapped from OTUs to contigs
This step was used to assess bins obtained from MEGAHIT assembled contigs using AMBER
`python gsmapping --readmapping read_otuid.tsv --contigmapping read_contig_mapped.tsv --length contig_length --outdir output_dir`

`read_otuids.tsv` - a `space` separated text file with read_id, otuid and tax_id

`read_contig_mapped.tsv` - a `space` separated text file with read_id and contig_id

`contig_length` - a `tab` separated text file with contig_id and length

## Split concatenated read file by sample id
This step was required to get sample-wise reads back from concatenated corrected reads output from CoCo.
`splitreadsbysample <sample_ids> <concatenatedreads.fastq> <outdir>`

`sample_ids` - a text file containing a list of sample ids
