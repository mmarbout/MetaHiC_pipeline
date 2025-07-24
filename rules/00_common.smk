#!/bin/env snakemake -s

# Rules to generates the build bowtie2 index and split fastq for alignement.

# Generate fasta without small contigs
rule remove_small_contigs:
    input:
        fasta = join(REF_DIR, '{ref}.fa'),
    output: 
        fasta = join(TMP, 'ref', '{ref}.fa'),
    params:
        min_size = config['metator']['min_size'],
    threads: 1
    conda: "../envs/metator_binning.yaml"
    shell: 'seqkit seq -m {params.min_size} {input.fasta} > {output.fasta}'

# Build bowtie2 index of the reference genome.
rule bt2_index_metator:
    input: join(TMP, 'ref', '{ref}.fa'),
    output: touch(join(TMP, 'ref', '{ref}.bt2_index.done')),
    params: idx = join(TMP, 'ref', '{ref}.genome'),
    threads: config['threads']
    conda: "../envs/metator_binning.yaml"
    shell: "bowtie2-build --threads {threads} {input} {params.idx}"

# Make splits from Hi-C fastq files to speed up mapping. [between 1 and 999].
rule split_fastq:
    input: join(FASTQ_DIR, '{library}_R{end}.fq.gz'),
    output: 
        expand(
            join(
                TMP,
                'split_reads',
                '{{library}}_R{{end}}',
                '{{library}}_R{{end}}.{split}.fq.gz'
            ),
            split=split_names
        ),
    params:
        n_splits = N_SPLITS,
        split_dir = join(TMP, 'split_reads', "{library}_R{end}"),
        tmp_dir = join(TMP, 'split_reads', "tmp", "{library}_R{end}"),
        fq_tmp = join(
            TMP,
            'split_reads',
            "tmp",
            "{library}_R{end}",
            '{library}_R{end}.fq'
        ),
        end = '{end}',
    message: "Splitting {wildcards.library}_R{wildcards.end} into {params.n_splits} split fastq."
    conda: "../envs/metator_binning.yaml"
    threads: config['threads_small']
    shell:
        """
        mkdir -p {params.split_dir}
        # Correct /1 -/2 at the end of the read ID
        mkdir -p {params.tmp_dir}
        set +euo pipefail
        ID=$(zcat {input} | head -n 1)
        ID2=$(zcat {input} | head -n 1 | cut -f1 -d '/')
        if [[ "$ID" == "$ID2" ]] ; then
            INPUT={input}
        else
            echo "Changing read IDs."
            zcat {input} | sed 's/\/{params.end}//' > {params.fq_tmp}
            INPUT={params.fq_tmp}
        fi
        # 100 split fastqs will be created with name pattern 001.fq - 100.fq
        seqkit split2 -p {params.n_splits} \
            -w 0 \
            -f \
            -e .gz \
            -j {threads} \
            -1 "$INPUT" \
            -O {params.split_dir}
        rm -r {params.tmp_dir}
        """

# Make splits from Hi-C fastq files to speed up mapping. [between 1 and 999].
rule split_fastq_digest:
    input: join(FASTQ_DIR, '{digest_library}_digest_R{end}.fq.gz'),
    output: 
        expand(
            join(
                TMP,
                'cutsite',
                '{{digest_library}}',
                '{{digest_library}}_{split}_R{{end}}.fq.gz',
            ),
            split=split_names
        ),
    params:
        n_splits = N_SPLITS,
        split_dir = join(TMP, 'cutsite', "{digest_library}"),
        tmp_dir = join(TMP, 'cutsite', "tmp", "{digest_library}"),
        fq_tmp = join(
            TMP,
            'cutsite',
            "tmp",
            "{digest_library}",
            '{digest_library}_R{end}.fq.gz'
        ),
        end = '{end}',
    message: "Splitting {wildcards.digest_library}_R{wildcards.end} into {params.n_splits} split fastq."
    conda: "../envs/metator_binning.yaml"
    threads: config['threads_small']
    shell:
        """
        mkdir -p {params.split_dir} {params.split_dir}_R{wildcards.end}

        # Correct /1 -/2 at the end of the read ID
        mkdir -p {params.tmp_dir}
        if [[ "{input}" == *nvq* ]] ; then
            echo "Changing read IDs."
            zcat {input} | sed 's/\/{params.end}//' > {params.fq_tmp}
            INPUT={params.fq_tmp}
        else
            INPUT={input}
        fi

        # 100 split fastqs will be created with name pattern 001.fq - 100.fq
        seqkit split2 -p {params.n_splits} \
            -w 0 \
            -f \
            -j {threads} \
            -1 "$INPUT" \
            -O {params.split_dir}_R{wildcards.end}
            
        for i in $(seq -f "%03g" 1 {params.n_splits})
            do
            mv \
                {params.split_dir}_R{wildcards.end}/{wildcards.digest_library}_digest_R{wildcards.end}.part_"$i".fq.gz \
                {params.split_dir}/{wildcards.digest_library}_part_"$i"_R{wildcards.end}.fq.gz
        done
        rm -rf {params.tmp_dir} {params.split_dir}_R{wildcards.end}
        """
