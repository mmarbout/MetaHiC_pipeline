#!/bin/env snakemake -s

# Rules to bin the scaffolded genomes.

# Rules to build louvain functions:
rule build_louvain:
    input: join(BASE_DIR, 'meta_analysis', 'external', 'louvain-generic.tar.gz')
    output: touch(join(BASE_DIR, 'meta_analysis', 'external', 'louvain.done'))
    params: 
        louvain_dir = join(BASE_DIR, 'meta_analysis', 'external')
    threads: 1
    conda: "../envs/metator_binning.yaml"
    shell:
        """
        cd {params.louvain_dir}
        tar -xzf louvain-generic.tar.gz
        cd gen-louvain
        make
        """

# Build pyfastx index from the fasta file.
rule build_pyfastx_index:
    input: join(TMP, 'ref', '{ref}.fa'),
    output: touch(join(REF_DIR, '{ref}.pyfastx_index.done')),
    threads: 1
    conda: '../envs/metator_binning.yaml'
    shell: 'pyfastx index {input}'

# Sort the pair using pairtools and pypairix using metator wrapper.
rule sort_metator_pairs:
    input: 
        join(TMP, 'pairs', '{hic_library}_{split}.valid_idx_pcrfree.pairs'),
    output: 
        join(TMP, 'pairs', '{hic_library}_{split}.valid_idx_pcrfree_sorted.pairs.gz'),
    threads: config['threads']
    conda: '../envs/metator_binning.yaml'
    shell: 'metator pairs -t {threads} -rF {input}'

# Run MetaTOR binning.
rule metator_binning:
    input:
        pairs = lambda w: expand(
            join(TMP, 'pairs', '{hic_library}_{split}.valid_idx_pcrfree_sorted.pairs.gz'),
            hic_library=samples.library[np.logical_and(
                list(samples.type == 'hic'),
                list(samples.ref == w.ref),
            )],
            split=split_names,
        ),
        assembly = join(TMP, 'ref', '{ref}.fa'),
        louvain_check = join(BASE_DIR, 'meta_analysis', 'external', 'louvain.done'),
    params:
        algorithm = config['metator']["algorithm"],
        aligner = config['metator']["aligner"],
        mode = config['metator']["mode"],
        louvain_path = join(BASE_DIR, 'meta_analysis', 'external', 'gen-louvain'),
        enzyme = lambda w: samples.enzyme[np.logical_and(
            list(samples.type == 'hic'),
            list(samples.ref == w.ref),
        )][0],
        iterations = config['metator']["iter"],
        rec_iter = config['metator']["rec_iter"],
        normalization = config['metator']["normalization"],
        out_dir = join(OUT_DIR, '{ref}', 'metator'),
        overlap = config['metator']["overlap"],
        rec_overlap = config['metator']["rec_overlap"],
        size = config['metator']["size"],
        tmp_dir = join(TMP, 'metator', '{ref}'),
        prefix = '{ref}',
        pairs = lambda w: ','.join([','.join([
            join(
                TMP, 'pairs', f'{hic_library}_{split}.valid_idx_pcrfree_sorted.pairs.gz'
            ) for hic_library in samples.library[np.logical_and(
            list(samples.type == 'hic'),
            list(samples.ref == w.ref),
        )] for split in split_names])]),
    output:
        contig_data = join(OUT_DIR, '{ref}', 'metator', 'contig_data_final.txt'),
        network = join(OUT_DIR, '{ref}', 'metator', 'network.txt'),
        bin_summary = join(OUT_DIR, '{ref}', 'metator', 'bin_summary.txt'),
        binning = join(OUT_DIR, '{ref}', 'metator', 'binning.txt'),
    conda: '../envs/metator_binning.yaml'
    threads: config["threads"]
    benchmark: "benchmarks/metator/{ref}.benchmark.txt"
    shell:
        """
        rm -rf {params.out_dir} {params.tmp_dir}
        LOUVAIN_PATH={params.louvain_path}
        metator pipeline \
            --assembly {input.assembly} \
            --algorithm {params.algorithm} \
            --aligner {params.aligner} \
            --aligner-mode {params.mode} \
            --enzyme {params.enzyme} \
            --force \
            --iterations {params.iterations} \
            --rec-iter {params.rec_iter} \
            --normalization {params.normalization} \
            --outdir {params.out_dir} \
            --overlap {params.overlap} \
            --prefix {params.prefix} \
            --rec-overlap {params.rec_overlap} \
            --size {params.size} \
            --start pair \
            --threads {threads} \
            --tmpdir {params.tmp_dir} \
            -1 {params.pairs} 
        """

# Run quality check from MetaTOR.
rule qc_metator:
    input:
        assembly = join(TMP, 'ref', '{ref}.fa'),
        pairs = lambda w: expand(
            join(TMP, 'pairs', '{hic_library}_{split}.valid_idx_pcrfree_sorted.pairs.gz'),
            hic_library=samples.library[np.logical_and(
                list(samples.type == 'hic'),
                list(samples.ref == w.ref),
            )],
            split=split_names,
        ),
        contig_data = join(OUT_DIR, '{ref}', 'metator', 'contig_data_final.txt'),
        bin_summary = join(OUT_DIR, '{ref}', 'metator', 'bin_summary.txt'),
    params:
        enzyme = lambda w: np.unique(samples.enzyme[np.logical_and(
            list(samples.type == 'hic'),
            list(samples.ref == w.ref),
        )])[0],
        prefix = '{ref}',
        outdir = join(OUT_DIR, '{ref}', 'metator', 'qc'),
        tmpdir = join(TMP, 'metator', '{ref}'),
    output:
        join(OUT_DIR, '{ref}', 'metator', 'qc', '{ref}_camembert_plot.pdf'),
        join(OUT_DIR, '{ref}', 'metator', 'qc', '{ref}_event.pdf'),
    threads: 1
    conda: '../envs/metator_binning.yaml'
    shell:
        """
        metator qc \
            --assembly {input.assembly} \
            --bin-summary {input.bin_summary} \
            --contig-data {input.contig_data} \
            --enzyme {params.enzyme} \
            --outdir {params.outdir} \
            --prefix {params.prefix} \
            --plot \
            --tmpdir {params.tmpdir} \
            {input.pairs}
        """

# Run quality check from MetaTOR one per sample.
rule qc_metator_lib:
    input:
        assembly = join(TMP, 'ref', '{ref}.fa'),
        pairs = lambda w: expand(
            join(TMP, 'pairs', '{{hic_library}}_{split}.valid_idx_pcrfree_sorted.pairs.gz'),
            split=split_names,
        ),
        contig_data = join(OUT_DIR, '{ref}', 'metator', 'contig_data_final.txt'),
        bin_summary = join(OUT_DIR, '{ref}', 'metator', 'bin_summary.txt'),
    params:
        enzyme = lambda w: np.unique(samples.enzyme[np.logical_and(
            list(samples.type == 'hic'),
            list(samples.ref == w.ref),
        )])[0],
        prefix = '{hic_library}',
        outdir = join(OUT_DIR, '{ref}', 'metator', 'qc_{hic_library}'),
        tmpdir = join(TMP, 'metator', '{ref}_{hic_library}'),
    output:
        join(OUT_DIR, '{ref}', 'metator', 'qc_{hic_library}', '{hic_library}_camembert_plot.pdf'),
        join(OUT_DIR, '{ref}', 'metator', 'qc_{hic_library}', '{hic_library}_event.pdf'),
    threads: 1
    conda: '../envs/metator_binning.yaml'
    shell:
        """
        rm -rf {params.tmpdir}
        metator qc \
            --assembly {input.assembly} \
            --bin-summary {input.bin_summary} \
            --contig-data {input.contig_data} \
            --enzyme {params.enzyme} \
            --outdir {params.outdir} \
            --prefix {params.prefix} \
            --plot \
            --no-clean-up \
            --tmpdir {params.tmpdir} \
            {input.pairs}
        reads=$( \
            grep -v '#' {params.tmpdir}/*/*pairs | \
            cut -f1 | \
            sed 's/:..$//' | \
            sed 's/:..$//' | \
            sort -u | \
            wc -l \
        )
        echo "Number of reads before digestion: $reads"
        """

# CheckM annotation and completness evaluation.
rule run_checkM:
    input:
        bin_summary = join(OUT_DIR, '{ref}', 'metator', 'bin_summary.txt'),
    params:
        fasta_dir = join(OUT_DIR, '{ref}', 'metator', 'final_bin_unscaffold'),
        out_dir = join(OUT_DIR, '{ref}', 'checkM'),
        tmp_dir = '{ref}_metator_checkM_tmp',
        end = 'fa',
    output:
        checkm_out = join(OUT_DIR, '{ref}', 'checkM', 'output.tsv'),
    conda: '../envs/metator_binning.yaml'
    threads: config['threads']
    shell:
        """
        module add hmmer/3.3 pplacer/v1.1-alpha19 prodigal/2.6.3 CheckM/1.1.3
        mkdir -p {params.out_dir} {params.tmp_dir}

        rm -rf {params.tmp_dir}
        echo "data in progress"
        checkm tree \
            --threads {threads} \
            --extension {params.end} \
            {params.fasta_dir} \
            {params.tmp_dir}
        checkm tree_qa {params.tmp_dir} -o 1 -f {params.out_dir}/taxonomy.txt

        echo "checkM lineage marker set"
        checkm lineage_set \
            {params.tmp_dir} \
            {params.tmp_dir}/markers.txt

        echo "checkM analyse and qa"
        checkm analyze \
            --threads {threads} \
            --tmpdir {params.tmp_dir} \
            --extension {params.end} \
            {params.tmp_dir}/markers.txt \
            {params.fasta_dir} \
            {params.tmp_dir}
        checkm qa \
            --threads {threads} \
            --tmpdir {params.tmp_dir} \
            {params.tmp_dir}/markers.txt \
            {params.tmp_dir} \
            -o 2 > {params.out_dir}/checkM_results_complete.txt
        sed 's/ \+/ /g' {params.out_dir}/checkM_results_complete.txt \
            > {output.checkm_out}
        rm -r {params.tmp_dir}
        """

# GTDBtk annotation and completness evaluation.
rule run_gtdbtk:
    input:
        bin_summary = join(OUT_DIR, '{ref}', 'metator', 'bin_summary.txt'),
    params:
        outdir = join(OUT_DIR, '{ref}', 'gtdbtk'),
        prefix = '{ref}',
        fasta_dir = join(OUT_DIR, '{ref}', 'metator', 'final_bin_unscaffold'),
        tmpdir = join(TMP, '{ref}', 'gtdbtk'),
    output: join(OUT_DIR, '{ref}', 'gtdbtk', '{ref}.bac120.summary.tsv'),
    threads: config["threads_large"]
    conda: '../envs/metator_binning.yaml'
    shell:
        """
        module add \
            prodigal \
            R/3.6.2 \
            hmmer/3.3 \
            pplacer/v1.1-alpha19 \
            FastANI/1.33 \
            FastTree/2.1.11 \
            Mash/2.2 \
            GTDBTk/
        mkdir -p {params.tmpdir} 
        gtdbtk classify_wf \
            --genome_dir {params.fasta_dir} \
            --out_dir {params.outdir} \
            --extension fa \
            --prefix {params.prefix} \
            --cpus {threads} \
            --force \
            --tmpdir {params.tmpdir}
        """

# Generates scaffold for binning.
rule scaffold:
    input:
        assembly = join(TMP, 'ref', '{ref}.fa'),
        pairs = lambda w: expand(
            join(TMP, 'pairs', '{hic_library}_{split}.valid_idx_pcrfree_sorted.pairs.gz'),
            hic_library=samples.library[np.logical_and(
                list(samples.type == 'hic'),
                list(samples.ref == w.ref),
            )],
            split=split_names,
        ),
        contig_data = join(OUT_DIR, '{ref}', 'metator', 'contig_data_final.txt'),
        bin_summary = join(OUT_DIR, '{ref}', 'metator', 'bin_summary.txt'),
    params:
        junctions = config['scaffold']['junctions'],
        threshold = config['scaffold']['threshold'],
        window_size = config['scaffold']['window_size'],
        outdir_fasta = join(OUT_DIR, '{ref}', 'metator', 'scaffold'),
        outdir_tab = join(OUT_DIR, '{ref}', 'metator', 'scaffold_info'),
        input_fasta = join(OUT_DIR, '{ref}', 'metator', 'final_bin_unscaffold'),
    output:
        touch(join(OUT_DIR, '{ref}', 'metator', 'scaffolding.done')),
    threads: 1
    conda: '../envs/metator_binning.yaml'
    shell:
        """
        mkdir -p {params.outdir_fasta} {params.outdir_tab}
        for bin_full in $(ls {params.input_fasta}/*.fa | cut -f1 -d '.') ;
        do  
            bin=$(basename $bin_full)
            metator scaffold \
                -b $bin \
                -i {params.input_fasta}/"$bin".fa \
                -j {params.junctions} \
                -o {params.outdir_fasta}/"$bin".fa \
                -O {params.outdir_tab}/"$bin".txt \
                -T {params.threshold} \
                -w {params.window_size} \
                {input.pairs}
        done
        """
