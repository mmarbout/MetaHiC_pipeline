#!/bin/env snakemake -s

# snakemake --rulegraph | dot -Tsvg > images/rulegraph.svg
# snakemake --dag | dot -Tsvg > images/dag.svg

# This file can be run using snakemake. It was tested on snakemake 5.10.0.
# It orchestrates the analysis of several metaHiC samples.

import numpy as np
import pandas as pd
from os.path import join
from snakemake.utils import validate

# Set parameters.
shell.prefix("set -euo pipefail;")

# LOAD CONFIG FILES
configfile: 'config/config.yaml'

samples = pd.read_csv(
    config['samples'], 
    sep=';', 
    dtype=str,
    comment='#',
).set_index(['library'], drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

BASE_DIR = config['base_dir']
OUT_DIR = join(BASE_DIR, config['out_dir'])
TMP = join(BASE_DIR, config['tmp_dir'])
REF_DIR = config['ref_dir']
FASTQ_DIR = config['fastq_dir']
DB_DIR = join(BASE_DIR, config['db_dir'])

hic_library = np.unique(samples.library[samples.type == 'hic'])
sg_library = np.unique(samples.library[samples.type == 'sg'])
samples_name = np.unique(samples['ref'])

sample_ref = {}
for lib in hic_library:
    ref = samples.ref[lib]
    sample_ref[lib] = ref
sample_ref = pd.DataFrame.from_dict(sample_ref, orient='index', columns=['ref'])
sample_ref['hic_library'] = sample_ref.index

# remove sample with no MAGs to avoid error message while doing a meaningless
# QC.
qc_exception = [
    "mat2a_proxi",
    "Hum4f_HiC",
    "mat3a_proxi",
    "Hum1c_HiC",
    "Hum1b_HiC",
]
sample_ref_qc = sample_ref.copy()
for exception in qc_exception:
    sample_ref_qc = sample_ref_qc.drop(exception)

# Split fastq
N_SPLITS = config['n_splits']
split_names = [f'part_{s:03}' for s in range(1, N_SPLITS + 1)] #base names of split files

wildcard_constraints:
    hic_library = "|".join(np.unique(samples.library[samples.type == 'hic'])),
    undigest_library = "|".join(np.unique(samples.library[samples.digestion == 'NO'])),
    digest_library = "|".join(np.unique(samples.library[samples.digestion == 'YES'])),
    sg_library = "|".join(np.unique(samples.library[samples.type == 'sg'])),
    ref = "|".join(samples_name),
    reseq = "|".join(['_nvq', '-a_nvq', '-b_nvq', '-c_nvq', '_nxq', '_NVQ']),

# Pipeline sub-workflows
include: 'rules/00_common.smk'
include: 'rules/01_hic_processing.smk'
include: 'rules/02_metator_binning.smk'
include: 'rules/03_metavir.smk'
include: 'rules/04_reassemble_virus.smk'

rule all:
  input:
      # 02 - Metator binning.
      expand(
        join(OUT_DIR, '{ref}', 'metator', 'bin_summary.txt'),
        ref=samples_name,
      ),
      expand(
        join(OUT_DIR, '{ref}', 'metator', 'qc', '{ref}_camembert_plot.pdf'),
        zip, **sample_ref_qc,
      ),
      expand(
        join(OUT_DIR, '{ref}', 'metator', 'qc_{hic_library}', '{hic_library}_camembert_plot.pdf'),
        zip, **sample_ref_qc,
      ),
      expand(
        join(OUT_DIR, '{ref}', 'checkM', 'output.tsv'),
        ref=samples_name,
      ),
      expand(
          join(OUT_DIR, '{ref}', 'gtdbtk', '{ref}.bac120.summary.tsv'),
          ref=samples_name,
      ),
      # expand(
      #     join(OUT_DIR, '{ref}', 'metator', 'scaffolding.done'),
      #     ref=samples_name,
      # ),
      # 03 - Phage annotation and association
      # expand(
      #   join(OUT_DIR, '{ref}', "annotation", "phage_list_{algo}.txt"),
      #   ref=samples_name,
      #   algo=['virsorter', 'genomad'],
      # ),
      expand(
        join(OUT_DIR, '{ref}', 'metavir_MGE', 'phages_bin_summary.tsv'),
        ref=samples_name,
      )
      # expand(
      #   join(OUT_DIR, '{ref}', 'metavir', 'contact_map', 'metavir_00001', '1.fa'),
      #   ref=samples_name,
      # ),
      # 04 - Phage genome reconstruction
      # expand(
      #   join(TMP, '{ref}', 'metavir', 'contact_map', 'phage_reassembly.done'),
      #   ref=samples.ref[samples.type == "sg"],
      # ),
