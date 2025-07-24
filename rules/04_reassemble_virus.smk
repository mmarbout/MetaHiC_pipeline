#!/bin/env snakemake -s

# Rules to reassemble potential phage contigs to have clean phage from the 
# metagenomics samples. In order to do it we select shotgun reads mapping on the
# identified phage contigs. After checking the coverage and doing a subsample to
# not go upper than 100x coverage, we do a new with only these reads using
# spades. Then the new fasta is kept.

import itertools

os.makedirs(join(TMP, 'bam_sorted'), exist_ok=True)

# Mapping shotgun reads.
rule sg_align_split:
    input:
        R1 = join(TMP, 'split_reads', '{sg_library}{reseq}_R1', '{sg_library}{reseq}_R1.{split}.fq.gz'),
        R2 = join(TMP, 'split_reads', '{sg_library}{reseq}_R2', '{sg_library}{reseq}_R2.{split}.fq.gz'),
        index_flag = lambda w: join(TMP, 'ref', f'{samples.loc[w.sg_library, "ref"]}' + '.bt2_index.done'),
    params:
        index = lambda w: join(TMP, 'ref', samples.loc[w.sg_library, 'ref'] + '.genome'),
    output: 
        bam = join(TMP, 'sg_align', '{sg_library}{reseq}', '{sg_library}{reseq}_{split}.bam'),
    threads: config['threads']
    conda: "../envs/metator_binning.yaml"
    shell:
        """
        bowtie2 \
            -p {threads} \
            -x {params.index} \
            --very-sensitive-local \
            -1 {input.R1} \
            -2 {input.R2} |
            samtools sort -n -@ {threads} -o {output}
        """


# Merge splits from individual mapping jobs into one bam file per library
rule merge_split_sg_alignments:
    input:
        expand(
          join(TMP, 'sg_align', '{{sg_library}}{{reseq}}', '{{sg_library}}{{reseq}}_{split}.bam'),
          split=split_names,
        )
    output: join(TMP, 'sg_align', '{sg_library}{reseq}', '{sg_library}{reseq}.bam'), 
    threads: config['threads']
    conda: "../envs/metator_binning.yaml"
    shell: 'samtools merge -n -l 9 -@ {threads} {output} {input}'


# Merge resequencing.
rule merge_reseq_sg_alignments:
    input:
        lambda w: [join(TMP, 'sg_align', i.split('_R')[0], i.split('_R')[0] + '.bam') for i in list(filter(lambda x:w.sg_library in x, [k for k in os.listdir(FASTQ_DIR) if 'R1' in k]))],
    params:
        tmp = join(TMP, 'sg_align', '{sg_library}.bam'),
    output: 
        bam = join(OUT_DIR, '{ref}', 'sg_bam', '{sg_library}.bam'),
        idx = join(OUT_DIR, '{ref}', 'sg_bam', '{sg_library}.bam.bai'),
    threads: config['threads']
    conda: "../envs/metator_binning.yaml"
    shell: 
      """
      nb=$(ls {input} | wc -l)
      if [[ $nb -eq 1 ]] 
      then
          samtools sort -l 9 -@ {threads} -o {output.bam} {input}
      else
          samtools merge -n -@ {threads} {params.tmp} {input}
          samtools sort -l 9 -@ {threads} -o {output.bam} {params.tmp}
      fi
      samtools index -@ {threads} {output.bam}
      """

# Selecting shotgun reads mapping on our contigs.
rule sg_extract_reads:
    input:
        bins = join(
            OUT_DIR, '{ref}', 'metavir', 'checkV_bins', "quality_summary.tsv"
        ),
        bam = expand(
          join(OUT_DIR, '{{ref}}', 'sg_bam', '{sg_library}.bam'), 
          sg_library=samples.library[(samples.type == 'sg')]
        ),
    params:
        script = join(
            BASE_DIR, 'meta_analysis', 'scripts', 'extract_reads_sg.sh'
        ),
        ref_dir = join(OUT_DIR, '{ref}', 'metavir', 'contact_map'),
        tmp_dir = join(TMP, '{ref}', 'metavir', 'contact_map'),
        adapters = join(BASE_DIR, 'meta_analysis', 'external', 'adapters.fa'),
        threshold = config['min_phage_length'],
        minlen = lambda w: config[w.ref]['readlen'],
        trimqual = lambda w: config[w.ref]['trimqual'],
        threads = 4,
        jobid = join(TMP, '{ref}', 'jobid_sg.txt'),
    output:
        join(TMP, '{ref}', 'metavir', 'contact_map', 'sg_extract_reads.done')
    conda: "../envs/metator_binning.yaml"
    threads: 1
    shell:
        """
        bam=$(echo {input.bam} | sed 's/ /___/g')
        rm -f {params.jobid}
        for data in $(cut -f1-2 {input.bins} | grep MetaVIR | sed 's/\t/_/') 
            do i=$(echo $data | cut -f2 -d '_')
            k=$(printf "%05d\n" $i)
            length=$(echo $data | cut -f3 -d '_')
            threshold={params.threshold}
            if [ "$length" -gt "$threshold" ] ; then
                jobid=$(sbatch -J SG_$k -c {params.threads} --mem 8G \
                    {params.script} \
                        $k \
                        $i \
                        $bam \
                        {params.adapters} \
                        {params.minlen} \
                        {params.trimqual} \
                        {params.tmp_dir} \
                        {params.ref_dir} \
                        {params.threads} | \
                    cut -f4 -d ' ')
                echo $jobid >> {params.jobid}
            fi
        done
        let 'wait=2'
        cat {params.jobid}
        while [ $wait -gt 1 ]
            let 'wait=3'
            do for job_id in $(cat {params.jobid})
                do
                countR=$(sacct -j $job_id | egrep -c 'JobID|RUNNING|PENDING')
                countF=$(sacct -j $job_id | egrep -c 'JobID|FAILED')
                if [ $countR -gt 1 ] ; then
                    let 'wait=2'
                fi
                if [ $countF -gt 1 ] ; then
                    for job_id in $(cat {params.jobid})
                        do scancel $job_id
                    done
                    let 'wait=1'
                    echo $job_id" have failed. Cancel the others."
                fi
            done
            if [ $wait -eq 3 ] ; then
                touch {output}
                break
            fi
            if [ $wait -eq 1 ] ; then
                break
            fi
        done
        rm -f {params.jobid}
        """

# Sort bam:
rule sorted_hic_bam:
    input:
        bam = join(TMP, 'bam', '{hic_library}{reseq}_R{end}.bam'),
    output:
        sorted_bam = join(TMP, 'bam_sorted', '{hic_library}{reseq}_R{end}.bam'),
        index_bam = join(
            TMP, 'bam_sorted', '{hic_library}{reseq}_R{end}.bam.bai'
        ),
        check = touch(join(TMP, 'bam', '{hic_library}{reseq}_R{end}.sorted.done')),
    conda: '../envs/metator_binning.yaml'
    threads: config['threads']
    shell:
        """
        samtools sort -l 9 -@ {threads} -o {output.sorted_bam} {input.bam}
        samtools index -@ {threads} {output.sorted_bam}
        """


# Extract reads from HiC samples:
rule hic_extract_reads:
    input:
        bins = join(
            OUT_DIR, '{ref}', 'metavir', 'checkV_bins', "quality_summary.tsv"
        ),
        check = expand(
            join(TMP, 'bam', 'MM33_NVQ_R{end}.sorted.done'),
            end = ['1', '2'],
        ),
        fq_for = join(FASTQ_DIR, 'MM33_NVQ_R1.fq.gz'),
        fq_rev = join(FASTQ_DIR, 'MM33_NVQ_R2.fq.gz'),
    params:
        script = join(
            BASE_DIR, 'meta_analysis', 'scripts', 'extract_reads_hic.sh'
        ),
        bam = lambda w: list(itertools.chain(
            *[[
                join(TMP, 'bam_sorted', bam) for bam in os.listdir(
                    join(TMP, 'bam_sorted')
                ) if ((lib in bam) and ('.bai' not in bam))
            ] for lib in samples.library[np.logical_and(
                list(samples.type == 'hic'),
                list(samples['ref'] == w.ref)
            )]]
        )),
        ref_dir = join(OUT_DIR, '{ref}', 'metavir', 'contact_map'),
        tmp_dir = join(TMP, '{ref}', 'metavir', 'contact_map'),
        threshold = config['min_phage_length'],
        adapters = join(BASE_DIR, 'meta_analysis', 'external', 'adapters.fa'),
        enzyme = lambda w: config[w.ref]['enzyme'],
        minlen = lambda w: config[w.ref]['readlen'],
        trimqual = lambda w: config[w.ref]['trimqual'],
        threads = 4,
        jobid = join(TMP, '{ref}', 'jobid_hic.txt'),
    output:
        join(TMP, '{ref}', 'metavir', 'contact_map', 'hic_extract_reads.done'),
    conda: "../envs/metator_binning.yaml"
    threads: config['threads']
    shell:
        """
        rm -f {params.jobid}
        bam=$(echo {params.bam} | sed 's/ /___/g')
        for data in $(cut -f1-2 {input.bins} | grep MetaVIR | sed 's/\t/_/') 
            do i=$(echo $data | cut -f2 -d '_')
            k=$(printf "%05d\n" $i)
            length=$(echo $data | cut -f3 -d '_')
            threshold={params.threshold}
            if [ "$length" -gt "$threshold" ] ; then
                jobid=$(sbatch -J HC_$k -c {params.threads} --mem 8G \
                    {params.script} \
                        $k \
                        $i \
                        $bam \
                        {input.fq_for} \
                        {input.fq_rev} \
                        {params.adapters} \
                        {params.enzyme} \
                        {params.minlen} \
                        {params.trimqual} \
                        {params.tmp_dir} \
                        {params.ref_dir} \
                        {params.threads} | \
                    cut -f4 -d ' ')
                echo $jobid >> {params.jobid}
            fi
        done
        let 'wait=2'
        cat {params.jobid}
        while [ $wait -gt 1 ]
            let 'wait=3'
            do for job_id in $(cat {params.jobid})
                do
                countR=$(sacct -j $job_id | egrep -c 'JobID|RUNNING|PENDING')
                countF=$(sacct -j $job_id | egrep -c 'JobID|FAILED')
                if [ $countR -gt 1 ] ; then
                    let 'wait=2'
                fi
                if [ $countF -gt 1 ] ; then
                    for job_id in $(cat {params.jobid})
                        do scancel $job_id
                    done
                    let 'wait=1'
                    echo $job_id" have failed. Cancel the others."
                fi
            done
            if [ $wait -eq 3 ] ; then
                touch {output}
                break
            fi
            if [ $wait -eq 1 ] ; then
                break
            fi
        done
        rm -f {params.jobid}
        """

# Do a new assembly.
rule assemble_phage:
    input:
        sg_reads_extracted = join(
            TMP, '{ref}', 'metavir', 'contact_map', 'sg_extract_reads.done'
        ),
        hic_reads_extracted = join(
            TMP, '{ref}', 'metavir', 'contact_map', 'hic_extract_reads.done'
        ),
        bins = join(
            OUT_DIR, '{ref}', 'metavir', 'checkV_bins', "quality_summary.tsv"
        ),
    output:
        touch(join(
            TMP, '{ref}', 'metavir', 'contact_map', 'phage_reassembly.done')
        ),
    params:
        script = join(
            BASE_DIR, 'meta_analysis', 'scripts', 'spades_reassembly.sh'
        ),
        mem = config['spades_mem'],
        ref_dir = join(OUT_DIR, '{ref}', 'metavir', 'contact_map'),
        tmp_dir = join(TMP, '{ref}', 'metavir', 'contact_map'),
        threshold = config['min_phage_length'],
        threads = 8,
        jobid = join(TMP, '{ref}', 'jobid_spa.txt'),
    threads: config['threads']
    conda: '../envs/spades.yaml'
    shell:
        """
        rm -f {params.jobid}
        for data in $(cut -f1-2 {input.bins} | grep MetaVIR | sed 's/\t/_/') 
            do i=$(echo $data | cut -f2 -d '_')
            k=$(printf "%05d\n" $i)
            length=$(echo $data | cut -f3 -d '_')
            threshold={params.threshold}
            if [ "$length" -gt "$threshold" ] ; then
                jobid=$(sbatch -J HC_$k -c {params.threads} --mem 8G \
                    {params.script} \
                    $k \
                    $i \
                    {params.tmp_dir} \
                    {params.ref_dir} \
                    {params.threads} \
                    {params.mem} | \
                cut -f4 -d ' ')
                echo $jobid >> {params.jobid}
            fi
        done
        let 'wait=2'
        cat {params.jobid}
        while [ $wait -gt 1 ]
            let 'wait=3'
            do for job_id in $(cat {params.jobid})
                do
                countR=$(sacct -j $job_id | egrep -c 'JobID|RUNNING|PENDING')
                countF=$(sacct -j $job_id | egrep -c 'JobID|FAILED')
                if [ $countR -gt 1 ] ; then
                    let 'wait=2'
                fi
                if [ $countF -gt 1 ] ; then
                    for job_id in $(cat {params.jobid})
                        do scancel $job_id
                    done
                    let 'wait=1'
                    echo $job_id" have failed. Cancel the others."
                fi
            done
            if [ $wait -eq 3 ] ; then
                touch {output}
                break
            fi
            if [ $wait -eq 1 ] ; then
                break
            fi
        done
        rm -f {params.jobid}
        """
