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
### pooled assembly
`megahit --12 *_reads.fastq -t 64 --presets meta-sensitive -o megahit_out`

### sample-wise assembly

`megahit --12 <sample_id>_reads.fastq -t 64 --presets meta-sensitive -o <sample_id>_megahit_out`
`cat *_megahit_out/final.contigs.fa > allcontigs_concatenatedallsamples.fa`  (concatenate all sample assemblies into one master contig file)

## Prepare abundance file from aligner output
`mkdir samfiles`
### pooled assembly
`strobealign -t 64 --aemb megahit_out/final.contigs.fa --eqx --interleaved <sample_id>.fastq > samfiles/abundances_<sample_id>.tsv`
`strobealign -t 64 megahit_out/final.contigs.fa --eqx --interleaved <sample_id>.fastq | samtools view -h -o samfiles/<sample_id>_strobealign.sam`

### concatenated sample-wise assembly
`strobealign -t 64 --aemb allcontigs_concatenatedallsamples.fa --eqx --interleaved ${sample_id}.fastq > samfiles/abundances_<sample_id>.tsv`
`strobealign -t 64 allcontigs_concatenatedallsamples.fa --eqx --interleaved ${sample_id}.fastq |samtools view -h -o samfiles/<sample_id>_strobealign.sam`

## Binning
Refer to benchmarking_scripts.ipynb

# Post binning steps

## Reassembly

**Extract mapped reads for each bin**

For combined read fastq and mapfile

`extractreads <fullpath/binfastafolder> <allsample_mapfile> <all_reads.corr.fq>`

For sample-wise processing, end extracted fastq file will have reads from all samples

    for sample in samplelist;
    do
        extractreads fullpath/binfastafolder ${sample}_mapfile ${sample}_reads.corr.fq; 
    done

## Plotting
Refer to plots.ipynb
