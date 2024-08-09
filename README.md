# binning_benchmarking

# Read correction

Reads are corrected by CoCo (https://github.com/soedinglab/CoCo). The correction is efficitive if k-mer counts are computed from pooled reads from all samples.

**Example with cami2 strain-madness dataset**

`cat *_reads.fq > all_reads.fq`

**Compute k-mer counts using dsk tool** (https://github.com/GATB/dsk)

`dsk -file all_reads.fq -kmer-size 41` (output: all_reads.h5)

**Correct reads using CoCo**

`coco correction --reads all_reads.fq --counts all_reads.h5 --outdir allreads_coco_corrected` (output: all_reads.corr.reads.fq)

**Split reads by sample origin**

`splitreadsbysample sample_id reads.fq outdir`




# Post binning steps

## Reassembly

**Extract mapped reads for each bin**

For combined read fastq and mapfile

`extractreads fullpath/genomeface_results/ allsample_mapfile all_reads.corr.fq`

For sample-wise processing, end extracted fastq file will have reads from all samples

    for sample in samplelist;
    do
        extractreads fullpath/genomeface_results/ ${sample}_mapfile ${sample}_reads.corr.fq; 
    done


