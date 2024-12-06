# binning_benchmarking

Binning benchmarking involves the following steps

coassembly: `read correction -> assembly -> mapping -> generate_abundance_matrix -> binning -> assessment`

single-sample: `read correction -> single-sample_assembly -> single-sampleread_mapping -> generate_abundance_matrix -> binning -> assessment`
 
multi-sample: `read correction -> sample-wise_assembly -> pool_allsampleassembly_contigs -> mapping -> generate_abundance_matrix -> binning -> split_bins -> remove_redundantbins -> assessment`

## Read correction

Reads correction is done by CoCo (https://github.com/soedinglab/CoCo). The correction is efficitive if k-mer counts are computed using reads pooled from all samples.

**Concatenate reads for CoCo correction**

`cat *_reads.fq > all_reads.fq`

**Compute k-mer counts using dsk tool** (https://github.com/GATB/dsk)

`dsk -file all_reads.fq -kmer-size 41` (output: all_reads.h5)

**Correct reads using CoCo**

`coco correction --reads all_reads.fq --counts all_reads.h5 --outdir allreads_coco_corrected` (output: all_reads.corr.reads.fq)

**Split reads by sample origin**

`splitreadsbysample <sampleid_file> <fastq_file> <outdir>` (output: <sample_id>.fastq)

## Assembly
MEGAHIT must be installed (https://github.com/voutcn/megahit.git)

### pooled assembly
`megahit --12 *_reads.fastq -t 64 --presets meta-sensitive -o megahit_out`

### sample-wise assembly
`megahit --12 <sample_id>_reads.fastq -t 64 --presets meta-sensitive -o <sample_id>_megahit_out`

`cat *_megahit_out/final.contigs.fa > allcontigs_concatenatedallsamples.fa`  (concatenate all sample assemblies into one master file, make sure that headers of sample-wise assemblies are prefixed with sample id spearated by 'C')

## Mapping
Strobealign is the fast and accurate aligner. We used to obtain the abundance matrix. (https://github.com/ksahlin/strobealign.git)

`mkdir samfiles`

### pooled assembly

`strobealign -t 64 --aemb megahit_out/final.contigs.fa --eqx --interleaved <sample_id>.fastq > samfiles/abundances_<sample_id>.tsv`

`strobealign -t 64 megahit_out/final.contigs.fa --eqx --interleaved <sample_id>.fastq | samtools view -h -o samfiles/<sample_id>_strobealign.sam`

### concatenated sample-wise assembly
`strobealign -t 64 --aemb allcontigs_concatenatedallsamples.fa --eqx --interleaved ${sample_id}.fastq > samfiles/abundances_<sample_id>.tsv`

`strobealign -t 64 allcontigs_concatenatedallsamples.fa --eqx --interleaved ${sample_id}.fastq |samtools view -h -o samfiles/<sample_id>_strobealign.sam`

If you have used other aligners (eg. bowtie2, bwa-mem), use our in-house script

`samtools view samfiles/<sample_id>_strobealign.sam | aligner2counts samfiles <sample_id> --only-mapids`

## Generate abundance matrix
`python util/get_abundance_tsv.py -i <samfiles> -l <contig_length> -m <minlength|1000>`

contig_length is a tab separated `.txt` file that should contain contig ids and length (contig_id\tlength). This file can be generated using `convertfasta_multi2single` executable (see README.md in `util/`).

inputdir is the directory of sample-wise abundance.tsv file. `abundances_<sample_id>.tsv`

### Sort alignment files
`samtools sort samfiles/<sample_id>_strobealign.sam -o samfiles/<sample_id>_strobealign_sorted.bam`

## Binning
Refer to benchmarking_scripts.ipynb. Ensure the order of contigs in the abundance matrix matches the assembly FASTA file.

## Split bins (multi-sample binning)
By default, most deep learning methods can split bins by sample id in multi-sample binning mode (McDevol, VAMB and GenomeFace). But tools such as COMEBin and MetaBAT2 don't have an option for it. To perform splitting, use our script in `util/`.

`python splitfasta_bysampleids.py --input_dir <bindir> --output_dir <outputdir> --format <binformat|fasta>`

This script assumes that sample id is located in-between `S` and `C` characters. For example, from a contig id `S1C141_284`, it will detect `1` as the sample id.

## Remove redundancy (multi-sample binning)
For this benchmarking, we mapped bins to source genomes for AMBER assessment as described in README.md in `util/`. However, it can be performed with de-replication approach `dRep` (https://github.com/MrOlm/drep). We leave the choice to the user.

## Assessment
### CheckM2
CheckM2 is a neural network-based method that estimates bin completeness and purity. (https://github.com/chklovski/CheckM2.git)

`checkm2 predict --input <binning_tool>_results -o <binning_tool>_results/checkm2_results --thread 24 -x fasta`

### AMBER
For the binning of contigs from gold-standard sets, we used AMBER assessment. (https://github.com/CAMI-challenge/AMBER.git)

`amber.py <binning_tool>_cluster.tsv -g gsa_pooled_mapping_short.binning -o amber_results`

where gsa_pooled\_mapping\_short.binning files for marine, strain-madness and plant-associated datasets were obtained from the CAMI2 assessment study.

### CheckM
CheckM is used to validate MetaBAT2 and MetaWRAP bin_refinement results. (https://github.com/Ecogenomics/CheckM.git)

`checkm lineage_wf <binning_tool>_results <binning_tool>_results/checkm_results -x fasta -t 24`

## Reassembly (post-binning refinement)

**Extract reads for each bin**

For combined read fastq and mapfiles

`extractreads <fullpath/binfastafolder> <allsample_mapids> <all_reads.corr.reads.fq> -f (binformat|fasta)`

For sample-wise processing

    for sample in samplelist;
    do
        extractreads fullpath/binfastafolder ${sample}_mapids ${sample}.fastq -f (binformat|fasta);
    done

`allsample_mapids` is a text file containing mapped read_id and contig_id separated by `tab`. This file can be generated using `aligner2counts` executable (see README.md in `util/`) for each sample as `<sample_id>_mapids`. Concatenate these sample-wise mapids to obtain `allsample_mapids`. From extracreads run, you will get `<bin_id>.fastq`. The extracted <bin_id>.fastq file will have reads from all samples.

`spades.py --12 <bin_id>.fastq --trusted-contigs <bin_id>.fasta --only-assembler --careful -o <bin_id>_assembly/ -t 12 -m 128`
`SPAdes` must be installed.

Refer to `workflow_reassemble` workflow to run the entire steps in a single run.

## Plotting
Refer to plots.ipynb
