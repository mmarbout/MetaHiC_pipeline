# Phage Host association in metaHiC samples

The goal of this pipeline is to bin contigs in MAGs abd viral bins, characterize
the viral genomes obtained and associate them with their host using several 
metaHiC smaples. The analysis done in this repository have been described in:  

Phages with a broad host range are common across ecosystems
 

## Usage

To download the pipeline.

```sh
git clone https://github.com/ABignaud/meta_analysis.git
```

To run the analysis, update the base_dir in the config file and download fastq 
files in a fastq_dir folder, and the corresponding assembly in a fasta_dir 
folder. Once it is done, run the following command:  

```sh
cd meta_analysis
snakemake --conda-create-envs-only --use-conda -j 1
snakemake -j 32 --use-conda
```



