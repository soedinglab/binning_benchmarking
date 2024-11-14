# binning_benchmarking

## Read correction

Reads are corrected by CoCo (https://github.com/soedinglab/CoCo). The correction is efficitive if k-mer counts are computed from pooled reads from all samples.

**Example with cami2 strain-madness dataset**

`cat *_reads.fq > all_reads.fq`

**Compute k-mer counts using dsk tool** (https://github.com/GATB/dsk)

`dsk -file all_reads.fq -kmer-size 41` (output: all_reads.h5)

**Correct reads using CoCo**

`coco correction --reads all_reads.fq --counts all_reads.h5 --outdir allreads_coco_corrected` (output: all_reads.corr.reads.fq)

**Split reads by sample origin**

`splitreadsbysample <sampleid_file> <fastq_file> <outdir>` (output: <sample_id>.fastq)

## Assemble reads
Requires MEGAHIT installed (https://github.com/voutcn/megahit.git)
### pooled assembly
`megahit --12 *_reads.fastq -t 64 --presets meta-sensitive -o megahit_out`

### sample-wise assembly

`megahit --12 <sample_id>_reads.fastq -t 64 --presets meta-sensitive -o <sample_id>_megahit_out`

`cat *_megahit_out/final.contigs.fa > allcontigs_concatenatedallsamples.fa`  (concatenate all sample assemblies into one master contig file)

## Prepare abundance file from aligner output
Strobealign is the fast and accurate aligner and we used to obtain abundance matrix. (https://github.com/ksahlin/strobealign.git)

`mkdir samfiles`

### pooled assembly

`strobealign -t 64 --aemb megahit_out/final.contigs.fa --eqx --interleaved <sample_id>.fastq > samfiles/abundances_<sample_id>.tsv`

`strobealign -t 64 megahit_out/final.contigs.fa --eqx --interleaved <sample_id>.fastq | samtools view -h -o samfiles/<sample_id>_strobealign.sam`

### concatenated sample-wise assembly
`strobealign -t 64 --aemb allcontigs_concatenatedallsamples.fa --eqx --interleaved ${sample_id}.fastq > samfiles/abundances_<sample_id>.tsv`

`strobealign -t 64 allcontigs_concatenatedallsamples.fa --eqx --interleaved ${sample_id}.fastq |samtools view -h -o samfiles/<sample_id>_strobealign.sam`

If you have used other aligners (eg. Bowtie2, bwa-mem), use our in-house script

`samtools view samfiles/<sample_id>_strobealign.sam | aligner2counts samfiles <sample_id> --only-mapids`

### Generate abundance matrix
`python util/get_abundance_tsv.py -i <inputdir> -l <contigslength> -m <minlength|1000>`

contigslength is a tab separated textfile that should contain contig ids and length (contig_id\tlength)

inputdir is the directly of sample-wise abundance.tsv file.`<samfiles/>` 

## Binning
Refer to benchmarking_scripts.ipynb. Make sure order of contigs data in abundance matrix and assembly file are the same as GenomeFace assumes so by default.

# Post binning steps

## Reassembly

**Extract mapped reads for each bin**

For combined read fastq and mapfile

`extractreads <fullpath/binfastafolder> <allsample_mapfile> <all_reads.corr.reads.fq>`

For sample-wise processing, end extracted fastq file will have reads from all samples

    for sample in samplelist;
    do
        extractreads fullpath/binfastafolder ${sample}_mapfile ${sample}.fastq;
    done

Refer to README of workflow to run in the entire steps of reassembly with a single run.

## Assessment
### CheckM2
CheckM2 is a neural network-based method that estimates bin completeness and purity reliably. (https://github.com/chklovski/CheckM2.git) \
`checkm2 predict --input <binning_tool>_results -o <binning_tool>_results/checkm2_results --thread 24 -x fasta`

### AMBER
For the binning of contigs from gold-standard sets, we used AMBER assessment tool. (https://github.com/CAMI-challenge/AMBER.git) \
`amber.py <binning_tool>_cluster.tsv -g gsa_pooled_mapping_short.binning -o amber_results` where gsa_pooled\_mapping\_short.binning files for marine, strain-madness and plant-associated datasets are provided from CAMI2 assessment study.

### CheckM
CheckM is used to validate MetaBAT2 and MetaWRAP bin_refinement results. (https://github.com/Ecogenomics/CheckM.git) \
`checkm lineage_wf <binning_tool>_results <binning_tool>_results/checkm_results -x fasta -t 24`

## Plotting
Refer to plots.ipynb
