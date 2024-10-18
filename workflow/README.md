### Snakemake Workflow for Benchmarking Binning Tools

Create a Snakemake environment (preferably using conda, [https://anaconda.org/bioconda/snakemake](https://anaconda.org/bioconda/snakemake)) to run the Snakemake pipeline. To reproduce the results, install the binning tools and assessment tool used here, namely:

- **McDevol** (v0.1.0, [git@github.com:soedinglab/McDevol.git](git@github.com:soedinglab/McDevol.git) *not yet public*)
- **VAMB** (v4.1.3, [https://github.com/RasmussenLab/vamb.git](https://github.com/RasmussenLab/vamb.git))
- **COMEBin** (v1.0.4, [https://github.com/ziyewang/COMEBin.git](https://github.com/ziyewang/COMEBin.git))
- **GenomeFace** (prerelease, [https://gist.github.com/richardlett/f56156631ea8167983f6e670b119a767](https://gist.github.com/richardlett/f56156631ea8167983f6e670b119a767))
- **MetaBAT2** (v2.17, [https://anaconda.org/bioconda/metabat2](https://anaconda.org/bioconda/metabat2))
- **CheckM2** (v1.0.3, [https://github.com/chklovski/CheckM2.git](https://github.com/chklovski/CheckM2.git))

Install these tools as instructed in their respective GitHub repositories or use the `.yml` files in the `environments` folder for the conda environments used in this study. For McDevol, set the McDevol download path in `config.yaml`.

#### To Run the Pooled Assembly Binning Pipeline
Use the command below:

```
snakemake --config dataset=<datasetname> mode=pooled threads=24 readpath=<fastaqfilespath> outputpath=<outputpath> minlength=1000 --cores=24 --use-conda
```
#### To Run the Multisample Assembly Binning Pipeline
Navigate to \`workflow_multisample\` and use the command below:

```
snakemake --config dataset=<datasetname> mode=multisample threads=24 readpath=<fastaqfilespath> outputpath=<outputpath> minlength=1000 --cores=24 --use-conda
```