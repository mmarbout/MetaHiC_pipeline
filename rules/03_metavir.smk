#!/bin/env snakemake -s

# Rules to run metavir which works in two steps:
#   - The first one is to associate the host MAGs with the phage contig.
#   - The second is to bin the phage contigs together if they belongs to the
#     same species.
#
# In order to run metavir some phage contigs annotation is necessary. For
# that we use virsorter2 and checkV. The methods is adapted from this tutorial:
# https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v3

# Setup annotation pipelines. Not possible on the cluster as they are no 
# internet connection. In local you can add the locations of the databases as 
# arguments for visrsorte2 and checkV 
rule setup_virsorter:
    output: 
        touch(join(DB_DIR, 'virsorter_db.done'))
    params: 
        db_dir = DB_DIR
    conda: '../envs/phage_annotation.yaml'
    threads: config['threads_small']
    shell: 'virsorter setup -d {params.db_dir}/db-vs2 -j {threads}'

rule setup_checkv:
    output: 
        touch(join(DB_DIR, 'checkv_db.done'))
    params: 
        db_dir = DB_DIR
    conda: '../envs/phage_annotation.yaml'
    threads: 1
    shell: 'checkv download_database {params.db_dir}/checkv'

rule setup_genomad:
    output: 
        touch(join(DB_DIR, 'genomad_db.done'))
    params: 
        db_dir = DB_DIR
    conda: '../envs/phage_annotation.yaml'
    threads: 1
    shell: 'genomad download-database {params.db_dir}'

# Virsorter2 first pass to extract viral contigs.
rule run_virsorter:
    input:
        assembly = join(TMP, 'ref', '{ref}.fa'),
        # db_flag = join(DB_DIR, 'virsorter_db.done'),
    output:
        fasta = join(OUT_DIR, '{ref}', "annotation", "virsorter", "final-viral-combined.fa"),
        score = join(OUT_DIR, '{ref}', "annotation", "virsorter", "final-viral-score.tsv"),
    params:
        out_dir = join(OUT_DIR, '{ref}', "annotation", "virsorter"),
        # db_dir = join(DB_DIR, 'db-vs2'),
    conda: '../envs/phage_annotation.yaml'
    threads: config['threads_large']
    shell:
        """
        module add R hmmer prodigal blast+ muscle mcl diamond VirSorter/2.2.3
        rm -r {params.out_dir}
        virsorter run \
            --keep-original-seq \
            -i {input.assembly} \
            -w {params.out_dir} \
            --include-groups dsDNAphage,ssDNA \
            --min-length 1000 \
            --min-score 0.3 \
            -j {threads} \
            --rm-tmpdir \
            all
        """

# Launch genomad.
rule run_genomad:
    input:
        assembly = join(TMP, 'ref', '{ref}.fa'),
        # db_flag = join(DB_DIR, 'genomad_db.done'),
    output:
        fasta = join(
            OUT_DIR, '{ref}', "annotation", "genomad", "{ref}_summary", "{ref}_virus.fna"
        ),
        score = join(
            OUT_DIR, '{ref}', "annotation", "genomad", "{ref}_summary", "{ref}_virus_summary.tsv"
        ),
        plasmid = join(
            OUT_DIR, '{ref}', "annotation", "genomad", "{ref}_summary", "{ref}_plasmid_summary.tsv"
        ),
    params:
        outdir = join(OUT_DIR, '{ref}', "annotation", "genomad"),
        dbdir = join(DB_DIR, 'genomad_db'),
    conda: '../envs/phage_annotation.yaml'
    threads: config['threads']
    shell:
        """
        rm -r {params.outdir}
        genomad end-to-end \
            --cleanup \
            --restart \
            --composition metagenome \
            --threads {threads} \
            {input.assembly} \
            {params.outdir} \
            {params.dbdir} \
        """   

# CheckV to clean non-viral sequences.
rule run_checkv_virsorter:
    input:
        phage_contigs = join(
            OUT_DIR,
            '{ref}',
            "annotation",
            "virsorter",
            "final-viral-combined.fa"
        ),
        # db_flag = join(DB_DIR, 'checkv_db.done'),
    output:
        fasta = join(OUT_DIR, '{ref}', "annotation", "checkv_virsorter", "combined.fa"),
        summary = join(OUT_DIR, '{ref}', "annotation", "checkv_virsorter", "quality_summary.tsv"),
    params:
        out_dir = join(OUT_DIR, '{ref}', "annotation", "checkv_virsorter"),
        db_dir = join(DB_DIR, 'checkv-db-v1.5'),
    conda: '../envs/phage_annotation.yaml'
    threads: config['threads']
    shell:
        """
        checkv end_to_end \
            {input.phage_contigs} \
            {params.out_dir} \
            -d {params.db_dir} \
            -t {threads} 
        cat \
            {params.out_dir}/proviruses.fna \
            {params.out_dir}/viruses.fna > \
            {output.fasta}
        """

# CheckV to clean non-viral sequences.
rule run_checkv_genomad:
    input:
        phage_contigs = join(
            OUT_DIR, '{ref}', "annotation", "genomad", "{ref}_summary", "{ref}_virus.fna"
        ),
        # db_flag = join(DB_DIR, 'checkv_db.done'),
    output:
        fasta = join(OUT_DIR, '{ref}', "annotation", "checkv_genomad", "combined.fa"),
        summary = join(OUT_DIR, '{ref}', "annotation", "checkv_genomad", "quality_summary.tsv"),
    params:
        out_dir = join(OUT_DIR, '{ref}', "annotation", "checkv_genomad"),
        db_dir = join(DB_DIR, 'checkv-db-v1.5'),
    conda: '../envs/phage_annotation.yaml'
    threads: config['threads']
    shell:
        """
        checkv end_to_end \
            {input.phage_contigs} \
            {params.out_dir} \
            -d {params.db_dir} \
            -t {threads} 
        cat \
            {params.out_dir}/proviruses.fna \
            {params.out_dir}/viruses.fna > \
            {output.fasta}
        """

# The two last step of the tutorial are not run as we don't care about
rule phage_curation_virsorter:
    input:
        virsorter = join(OUT_DIR, '{ref}', "annotation", "virsorter", "final-viral-score.tsv"),
        checkv = join(OUT_DIR, '{ref}', "annotation", "checkv_virsorter", "quality_summary.tsv"),
    params:
        sort_virsorter = join(OUT_DIR, '{ref}', "annotation", "virsorter_sorted.tsv"),
        sort_checkv = join(OUT_DIR, '{ref}', "annotation", "checkv_sorted.tsv"),
    output:
        table = join(OUT_DIR, '{ref}', "annotation", "phage_table_virsorter.tsv"),
        liste = join(OUT_DIR, '{ref}', "annotation", "phage_list_virsorter.txt"),
    conda: '../envs/phage_annotation.yaml'
    threads: 1
    shell:
        """
        # Join virsorter and checkV output. Remove partial contig, provirus, 
        # contigs below 1kb, and the ones with no viral genes.
        grep -v 'seqname' {input.virsorter} | sort -k 1,1 > {params.sort_virsorter}
        grep -v contig_id {input.checkv} | sort -k 1,1 > {params.sort_checkv}
        join {params.sort_virsorter} {params.sort_checkv} | \
            sed 's/ /\t/g' | \
            grep -v '_partial' | \
            grep -v 'Yes' | \
            awk '$10 > 1000 {{print $0}}' | \
            awk '($14 > 0 || ($15 == 0 || $4 >= 0.95 || $7 > 2)) {{print $0}}' > \
            {output.table}
            cut -f1 -d '|' {output.table} > {output.liste}
        rm {params.sort_virsorter} {params.sort_checkv}
        """

# The two last step of the tutorial are not run as we don't care about
res_finder_dir = config["resfinderDir"]
rule phage_curation_genomad:
    input:
        genomad = join(
            OUT_DIR, '{ref}', "annotation", "genomad", "{ref}_summary", "{ref}_virus_summary.tsv"
        ),
        plasmid = join(
            OUT_DIR, '{ref}', "annotation", "genomad", "{ref}_summary", "{ref}_plasmid_summary.tsv"
        ),
        resfinder = join(res_finder_dir, 'ResFinder_{ref}/results_tab.txt'),
        checkv = join(OUT_DIR, '{ref}', "annotation", "checkv_genomad", "quality_summary.tsv"),
    params:
        sort_genomad = join(OUT_DIR, '{ref}', "annotation", "genomad_sorted.tsv"),
        sort_checkv = join(OUT_DIR, '{ref}', "annotation", "checkv_sorted.tsv"),
        liste = join(OUT_DIR, '{ref}', "annotation", "tmp.MGE_list_genomad.txt"),
    output:
        table = join(OUT_DIR, '{ref}', "annotation", "phage_table_genomad.tsv"),
        liste = join(OUT_DIR, '{ref}', "annotation", "MGE_list_genomad.txt"),
    conda: '../envs/phage_annotation.yaml'
    threads: 1
    shell:
        """
        # Join genomad and checkV output. Remove partial contig, provirus, 
        # contigs below 1kb, and the ones with no viral genes.
        set +eu
        grep -v 'seq_name' {input.genomad} | sort -k 1,1 > {params.sort_genomad}
        grep -v contig_id {input.checkv} | sort -k 1,1 > {params.sort_checkv}
        join {params.sort_genomad} {params.sort_checkv} | \
            sed 's/No terminal repeats/No_terminal_repeats/' | \
            sed 's/ /\t/g' | \
            grep -v 'Provirus' | \
            grep -v 'Yes' | \
            awk '$12 > 1000 {{print $0}}' | \
            awk '($16 > 0 || ($17 == 0 || $6 >= 0.7 || $9 > 2)) {{print $0}}' > \
            {output.table}
            cut -f1 -d '|' {output.table} > {params.liste}
        grep -v 'seq_name' {input.plasmid} | \
            sort -k 1,1 | \
            sed 's/No terminal repeats/No_terminal_repeats/' | \
            awk '$2 > 1000 {{print $0}}' | \
            awk '$4 > 1 {{print $0}}' | \
            awk '$6 > 0.5 {{print $0}}' | \
            cut -f1 >> {params.liste}
        grep -v 'Resistance' {input.resfinder} | \
            cut -f6 >> {params.liste}
        sort -u {params.liste} > {output.liste}
        rm {params.sort_genomad} {params.sort_checkv} 
        """

# Binning of the phage contigs.
rule metavir_run:
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
        phages = join(OUT_DIR, '{ref}', "annotation", "MGE_list_genomad.txt"),
        binning = join(OUT_DIR, '{ref}', "metator", "binning.txt"),
        contig = join(OUT_DIR, '{ref}', "metator", "contig_data_final.txt"),
        network = join(OUT_DIR, '{ref}', "metator", "network.txt"),
        # db_flag = join(DB_DIR, 'checkv_db.done'),
    params:
        pairs = lambda w: ','.join([
            join(TMP, 'pairs', f'{hic_library}_{split}.valid_idx_pcrfree_sorted.pairs.gz') 
            for hic_library in samples.library[np.logical_and(
                list(samples.type == 'hic'),
                list(samples.ref == w.ref),
            )] for split in split_names]),
        outdir = join(OUT_DIR, '{ref}', 'metavir_MGE'),
        tmpdir =  join(TMP, '{ref}', 'metavir'),
        checkv_db = join(DB_DIR, 'checkv-db-v1.5'),
        threshold = config['metavir']['threshold'],
    output:
        phages_data = join(OUT_DIR, '{ref}', 'metavir_MGE', 'phages_bin_summary.tsv'),
        bins = join(
            OUT_DIR, '{ref}', 'metavir_MGE', 'checkV_bins', "quality_summary.tsv"
        ),
    threads: config['threads']
    conda: '../envs/phage_annotation.yaml'
    shell:
        """
        module add diamond/2.0.4 hmmer/3.3 prodigal-gv/2.9.0 CheckV/1.0.1
        metavir binning \
            --binning {input.binning} \
            --contigs-data {input.contig} \
            --network {input.network} \
            --pairs {params.pairs} \
            --fasta {input.assembly} \
            --outdir {params.outdir} \
            --phages {input.phages} \
            --threshold {params.threshold} \
            --plot \
            --threads {threads} \
            --checkv-db {params.checkv_db} \
            --tmpdir {params.tmpdir}
        """

# Plot the contact map of large bins.
rule phage_contact_map:
    input:
        assembly = join(TMP, 'ref', '{ref}.fa'),
        contigs = join(OUT_DIR, '{ref}', 'metavir', "phages_data_final.tsv"),
        bins = join(
            OUT_DIR, '{ref}', 'metavir', 'checkV_bins', "quality_summary.tsv"
        ),
        pairs_idx = lambda w: [join(
            OUT_DIR,
            w.ref,
            'pairs',
            f'{hic_lib}_sorted.pairs.gz.px2'
        ) for hic_lib in samples.library[np.logical_and(list(samples.type == 'hic'), list(samples['ref'] == w.ref))]],
    output:
        fasta = join(OUT_DIR, '{ref}', 'metavir', 'contact_map', 'metavir_00001', '1.fa'),
    params:
        sorted_pairs = lambda w: ','.join([join(
            OUT_DIR,
            w.ref,
            'pairs',
            f'{hic_lib}_sorted.pairs.gz'
        ) for hic_lib in samples.library[np.logical_and(list(samples.type == 'hic'), list(samples['ref'] == w.ref))]]),
        enzyme = lambda w: samples.enzyme[samples['ref'] == w.ref][0],
        min_contig_size = config['min_contig_size'],
        res = config['res_phage'],
        threshold = config['min_phage_length'],
        balance_args = config['balance_args'],
        out_dir = join(OUT_DIR, '{ref}', 'metavir', 'contact_map'),
        tmp_dir = join(TMP, '{ref}', 'metavir_contact_map'),
        base_dir = BASE_DIR,
    threads: 1
    conda: "../envs/metator_binning.yaml"
    shell:
        """
        mkdir -p {params.out_dir}/0_plot/
        for data in $(cut -f1-2 {input.bins} | grep MetaVIR | sed 's/\t/_/') 
            do i=$(echo $data | cut -f2 -d '_')
            k=$(printf "%05d\n" $i)
            length=$(echo $data | cut -f3 -d '_')
            threshold={params.threshold}
            if [ "$length" -gt "$threshold" ] ; then
                metator contactmap \
                    -a {input.assembly} \
                    -c {input.contigs} \
                    -e {params.res} \
                    -n $i \
                    -p {params.sorted_pairs} \
                    --force \
                    --pcr-dup \
                    -m cool \
                    -s {params.min_contig_size} \
                    -O MetaVir_bin \
                    --tmpdir={params.tmp_dir}/metavir_"$k" \
                    -o {params.out_dir}/metavir_"$k" ;
                if test -f {params.out_dir}/metavir_"$k"/"$i".cool ; then
                    cooler balance \
                        {params.out_dir}/metavir_"$k"/"$i".cool \
                        --min-nnz 2 \
                        --ignore-diags 1
                    python3 {params.base_dir}/meta_analysis/scripts/plot_phage.py \
                        {params.out_dir}/metavir_"$k"/"$i".cool \
                        {params.res} \
                        {params.out_dir}/0_plot/"$k".png
                fi
            fi
        done
        """
