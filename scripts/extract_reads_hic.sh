#! /bin/bash

# Set variables
bin=$1
bin_simple=$2
bams=$3
fq_for=$4
fq_rev=$5
adapters=$6
enzyme=$7
min_len=$8
trimqual=$9
tmp_dir=${10}
out_dir=${11}
threads=${12}

echo "MetaVIR bin $bin in progress..."

# Check tmp directory
rm -rf "$tmp_dir"/metavir_hic_"$bin"
mkdir -p "$tmp_dir"/metavir_hic_"$bin"

# For each bam do the same.
echo 'Extract reads from bam...'
for bam in $(echo $bams | sed 's/___/ /g')
do
    contigs=$(grep '>' "$out_dir"/metavir_"$bin"/"$bin_simple".fa | \
        cut -f2 -d '>' | \
        cut -f1 -d ' ')
    for contig in $contigs
    do
        lib=$(basename $bam)
        # Extract reads.
        samtools view -@ "$threads" -u -P $bam $contig | \
            samtools sort \
            -n \
            -l 9 \
            -@ "$threads" \
            -o "$tmp_dir"/metavir_hic_"$bin"/"$contig"_"$lib" - 
    done
done

# Merge from all bam.
echo 'Merge bam...'
samtools merge \
    -@ "$threads" \
    -n \
    -o "$tmp_dir"/metavir_hic_"$bin"/merge.bam \
    "$tmp_dir"/metavir_hic_"$bin"/*_*.bam

# Extract fastq IDs.
echo 'Extract reads IDs...'
samtools view "$tmp_dir"/metavir_hic_"$bin"/merge.bam | \
    cut -f1 | \
    sed 's/:[0-9]*$//' | \
    sort -u > \
    "$tmp_dir"/metavir_hic_"$bin"/list_id.txt

# Extract undigested fastq.
echo 'Extract forward reads from fastq...'
seqtk subseq \
    "$fq_for" \
    "$tmp_dir"/metavir_hic_"$bin"/list_id.txt > \
    "$tmp_dir"/metavir_hic_"$bin"/hic_"$bin"_R1.fq

echo 'Extract reverse reads from fastq...'
seqtk subseq \
    "$fq_rev" \
    "$tmp_dir"/metavir_hic_"$bin"/list_id.txt > \
    "$tmp_dir"/metavir_hic_"$bin"/hic_"$bin"_R2.fq 

# Compressed fastq.
echo 'Compressed fastq...'
gzip "$tmp_dir"/metavir_hic_"$bin"/hic_"$bin"_R1.fq
gzip "$tmp_dir"/metavir_hic_"$bin"/hic_"$bin"_R2.fq

# Digest fastq.
echo 'Digest reads...'
hicstuff cutsite \
    -1 "$tmp_dir"/metavir_hic_"$bin"/hic_"$bin"_R1.fq.gz \
    -2 "$tmp_dir"/metavir_hic_"$bin"/hic_"$bin"_R2.fq.gz \
    -e $enzyme \
    -m all \
    -s 0 \
    -t "$threads" \
    -p "$tmp_dir"/metavir_hic_"$bin"/hc_"$bin"

# Clean the reads.
echo 'Clean reads...'
bbduk.sh \
    in1="$tmp_dir"/metavir_hic_"$bin"/hc_"$bin"_R1.fq.gz \
    in2="$tmp_dir"/metavir_hic_"$bin"/hc_"$bin"_R2.fq.gz \
    out1="$out_dir"/metavir_"$bin"/hc_cleaned_"$bin"_R1.fq.gz \
    out2="$out_dir"/metavir_"$bin"/hc_cleaned_"$bin"_R2.fq.gz \
    ref="$adapters" \
    ktrim=r \
    hdist=1 \
    tpe \
    tbo \
    minlen=$min_len \
    qtrim=rl \
    trimq=$trimqual
    
# Clean tmp directory
# rm -rf "$tmp_dir"/metavir_hic_"$bin"
