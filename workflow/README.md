### Snakemake workflow for benchmarking binning of pooled assembly contigs

Create a snakemake environment (preferrably using conda, https://anaconda.org/bioconda/snakemake) to run snakemake pipeline. To reproduce the results, install binning tools and asssment tool used here namely \
        - McDevol (v0.1.0, git@github.com:soedinglab/McDevol.git not yet public) \
        - VAMB (v4.1.3, https://github.com/RasmussenLab/vamb.git) \
        - COMEBin (v1.0.4, https://github.com/ziyewang/COMEBin.git) \
        - GenomeFace (prerelease, https://gist.github.com/richardlett/f56156631ea8167983f6e670b119a767) \
        - MetaBAT2 (v2.17, https://anaconda.org/bioconda/metabat2) \
        - CheckM2 (v1.0.3, https://github.com/chklovski/CheckM2.git) \
in your system as instructed in their respective github repositories or use `.yml` files in `environments` folder for the conda environments used in this study. For McDevol, set the mcdevol download path in `config.yaml`.


To run pooled assembly binning pipeline, use the command below
`snakemake --config dataset=<datasetname> mode=pooled threads=24 readpath=<fastaqfilespath> outputpath=<outputpath> minlength=1000 --cores=24 --use-conda`

To run multisample assembly binning pipeline, navigate to workflow_multisample use the command below
`snakemake --config dataset=<datasetname> mode=multisample threads=24 readpath=<fastaqfilespath> outputpath=<outputpath> minlength=1000 --cores=24 --use-conda`