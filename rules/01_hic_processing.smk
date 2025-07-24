#!/bin/env snakemake -s

# Rules to generates the HiC matrices.


# Alignment of a single fastq split from a Hi-C library with hicstuff iteralign.
rule cutsite_hic:
    input:
        R1 = join(
            TMP,
            'split_reads', 
            '{undigest_library}_R1', 
            '{undigest_library}_R1.{split}.fq.gz'
        ),
        R2 = join(
            TMP,
            'split_reads', 
            '{undigest_library}_R2', 
            '{undigest_library}_R2.{split}.fq.gz'
        ),
    output: 
        R1 = join(
            TMP,
            'cutsite', 
            '{undigest_library}', 
            '{undigest_library}_{split}_R1.fq.gz'
        ),
        R2 = join(
            TMP,
            'cutsite', 
            '{undigest_library}', 
            '{undigest_library}_{split}_R2.fq.gz'
        ),
    params:
        out_dir = join(TMP, 'cutsite', '{undigest_library}'),
        enzyme = lambda w: np.unique(samples.enzyme[w.undigest_library])[0],
        prefix = join(
            TMP, 'cutsite', '{undigest_library}', '{undigest_library}_{split}'
        ),
    threads: config['threads_small']
    conda: "../envs/metator_binning.yaml"
    shell:
        """
        mkdir -p {params.out_dir}
        hicstuff cutsite \
            --threads {threads} \
            --enzyme {params.enzyme} \
            --mode all \
            --prefix {params.prefix} \
            --forward {input.R1} \
            --reverse {input.R2}
        rm {input.R1} {input.R2}
        """

# Generate HiC pairs even for shotgun reads.
rule generates_pairs:
    input:
        index = lambda w: join(TMP, 'ref', f'{samples.loc[w.hic_library, "ref"]}.bt2_index.done'),
        R1 = join(TMP, 'cutsite', '{hic_library}', '{hic_library}_{split}_R1.fq.gz'),
        R2 = join(TMP, 'cutsite', '{hic_library}', '{hic_library}_{split}_R2.fq.gz'),
    params:
        fasta = lambda w: join(TMP, 'ref', f'{samples.loc[w.hic_library, "ref"]}.genome'),
        enzyme = lambda w: samples.enzyme[w.hic_library],
        prefix = '{hic_library}_{split}',
        qmin = config['metator']['min_qual'],
        outdir = join(TMP, 'pairs', '{hic_library}_{split}'),
        out_pairs = join(TMP, 'pairs'),
        tmpdir = join(TMP, 'pairs', '{hic_library}_{split}', 'tmp'),
        pairs = join(TMP, 'pairs', '{hic_library}_{split}', 'tmp', '{hic_library}_{split}.valid_idx_pcrfree.pairs'),
    output: 
        pairs = join(TMP, 'pairs', '{hic_library}_{split}.valid_idx_pcrfree.pairs'),
    threads: config['threads']
    conda: '../envs/metator_binning.yaml'
    shell:
        """
        mkdir -p {params.out_pairs}
        hicstuff pipeline \
            --duplicates \
            --enzyme {params.enzyme} \
            --force \
            --genome {params.fasta} \
            --no-cleanup \
            --outdir {params.outdir} \
            --prefix {params.prefix} \
            --quality-min {params.qmin} \
            --threads {threads} \
            --tmpdir {params.tmpdir} \
            {input.R1} \
            {input.R2}
        mv {params.pairs} {output.pairs}
        rm -rf {input.R1} {input.R2} {params.tmpdir} {params.outdir}
        """
