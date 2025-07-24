# Phage Host association in metaHiC samples

The goal of this pipeline is to bin contigs in MAGs and vMAGs, characterize
the viral genomes obtained and associate them with their host. 

The analysis have been described in:  

Bignaud et al.
Phages with a broad host range are common across ecosystems.
 
The different data mandatory to reproduce this analysis are all publicly available and described in the publication.


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



