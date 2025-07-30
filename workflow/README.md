### Snakemake Workflow for Benchmarking Binning Tools

Create a Snakemake environment (preferably using conda, [https://anaconda.org/bioconda/snakemake](https://anaconda.org/bioconda/snakemake)) to run the Snakemake pipeline. Singularity container platform should be preinstalled. To reproduce the results, install the binning tools and assessment tool used here, namely:

- **VAMB** (v4.1.3, [https://github.com/RasmussenLab/vamb.git](https://github.com/RasmussenLab/vamb.git))
- **SemiBin2** (v2.2.0, [https://github.com/BigDataBiology/SemiBin.git](https://github.com/BigDataBiology/SemiBin.git))
- **COMEBin** (v1.0.4, [https://github.com/ziyewang/COMEBin.git](https://github.com/ziyewang/COMEBin.git))
- **GenomeFace** (prerelease, [https://gist.github.com/richardlett/f56156631ea8167983f6e670b119a767](https://gist.github.com/richardlett/f56156631ea8167983f6e670b119a767))
- **TaxVAMB** (v5.0.4, [https://github.com/RasmussenLab/vamb.git](https://github.com/RasmussenLab/vamb.git))
- **MetaBAT2** (v2.17, [https://anaconda.org/bioconda/metabat2](https://anaconda.org/bioconda/metabat2))
- **CheckM2** (v1.0.3, [https://github.com/chklovski/CheckM2.git](https://github.com/chklovski/CheckM2.git))

Install these tools with version specified here as instructed in their respective repositories (recommended). Otherwise, use the `.yml` files in the `environments` folder for the conda environments used in this study. Reassembly module has been tested and can work on both Linux and MacOS systesm whereas workflows for coassembly multi-sample and multi-sample binning work only on Linux OS due to non-availability of some tools for MacOS.

####  Configuring paths for workflow
Set correct paths in `config.yaml` of respective workflows before start running.

`STROBEALIGNPATH`: Specifies the path to Strobealign executable file used for sequence alignment. (eg. <parentpath>/strobealign/build)

`UTILPATH`: Specifies the path to the utility directory and contains helper scripts for binning benchmarking tasks. It is located in download directory of binning_benchmarking. (eg. <downloadpath>/binning_benchmarking/util)

`ENVIRONMENT`: Specifies environment path to find environment file for binning tools (`.yml`), can be found at <downloadpath>/binning_benchmarking/workflow/environments.

`SPADESPATH`: Specifies the path to executable file of SPAdes assembler.

**! Note:** argument inputs for files and folders should contain `absolute path` for all workflows. If the conda environment for any tool fails to create during a Snakemake run, you can manually create it using the following command `conda env create --file=<tool>_env.yml`. After creating the environment manually, update the corresponding Snakemake rule to reference the environment by its name (`<tool_env>`) instead of the environment YAML path (`f{ENVIRONMENT}/comebin_env.yml`).  

#### To Run the coassembly binning pipeline
Use the command below:

```
snakemake --config dataset=<datasetname> mode=pooled threads=24 readpath=<fastaqfilespath> outputpath=<outputpath> minlength=1000 --cores=24 --use-conda
```
![workflow_coassembly](https://github.com/user-attachments/assets/971796e3-8bb5-4a64-87b6-22eb2f8befc3)
#### To run the multi-sample binning pipeline
Navigate to \`workflow_multisample\` and use the command below:

```
snakemake --config dataset=<datasetname> mode=multisample threads=24 readpath=<fastaqfilespath> outputpath=<outputpath> minlength=1000 --cores=24 --use-conda
```
Bins have to be de-replicated or selected the best bin per genome based on gold-standard mapping before subjecting them to assessment.
![workflow_multisample](https://github.com/user-attachments/assets/9f3a312c-07f4-40bf-8b69-6de72ec43099)

#### Post-binning reassembly bins
Reassembly is performed using contigs in the bin and reads mapped to those contigs from all samples using SPAdes assembler (https://github.com/ablab/spades/releases/tag/v4.0.0).

```
snakemake --config threads=24 binpath=<binpath> binformat=<fasta|faa|fa> sampath=<alginmentpath> readpath=<fastaqfilespath> readformat=<fastq|fq> outputpath=<outputpath> --cores 24 --use-conda
```
<img src="https://github.com/user-attachments/assets/4f311ced-602f-450b-9dd2-8d3cbcbd0d38" alt="workflow_reassemble" width="800"/>

#### Analysing assessment results
This study focused on couting number of near-complete (90% completeness and 5% contamination), higher quality (70% completeness and 10% contamination) and medium quality bins (50% completeness and 10% contamination) bins. Use the bash command below on CheckM2 output `quality_report.tsv`

```
awk '{if($2>=90&&$3<5) print $0}' quality_report.tsv | wc -l && awk '{if($2>=70&&$3<10) print $0}' quality_report.tsv | wc -l && awk '{if($2>=50&&$3<10) print $0}' quality_report.tsv | wc -l

```
