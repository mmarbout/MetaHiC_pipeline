#! /bin/bash

# Set variables
bin=$1
bin_simple=$2
bams=$3
adapters=$4
min_len=$5
trimqual=$6
tmp_dir=$7
out_dir=$8
threads=$9

echo "MetaVIR bin $bin in progress..."

# Check tmp directory
rm -rf "$tmp_dir"/metavir_"$bin"
mkdir -p "$tmp_dir"/metavir_"$bin"

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
            -o "$tmp_dir"/metavir_"$bin"/"$contig"_"$lib".bam - 
    done
done

# Merge from all bam.
echo 'Merge bam...'
samtools merge \
    -@ "$threads" \
    -n \
    -o "$tmp_dir"/metavir_"$bin"/merge.bam \
    "$tmp_dir"/metavir_"$bin"/*_*.bam

# Transform them to fastq.
echo 'Extract reads from fastq...'
bedtools bamtofastq \
    -i "$tmp_dir"/metavir_"$bin"/merge.bam \
    -fq "$tmp_dir"/metavir_"$bin"/sg_"$bin"_R1.fq \
    -fq2 "$tmp_dir"/metavir_"$bin"/sg_"$bin"_R2.fq
gzip "$tmp_dir"/metavir_"$bin"/sg_"$bin"_R1.fq
gzip "$tmp_dir"/metavir_"$bin"/sg_"$bin"_R2.fq

# Clean the reads.
echo $min_len
echo 'Clean reads...'
bbduk.sh \
    in1="$tmp_dir"/metavir_"$bin"/sg_"$bin"_R1.fq.gz \
    in2="$tmp_dir"/metavir_"$bin"/sg_"$bin"_R2.fq.gz \
    out1="$out_dir"/metavir_"$bin"/sg_cleaned_"$bin"_R1.fq.gz \
    out2="$out_dir"/metavir_"$bin"/sg_cleaned_"$bin"_R2.fq.gz \
    ref="$adapters" \
    ktrim=r \
    hdist=1 \
    tpe \
    tbo \
    minlen=$min_len \
    qtrim=rl \
    trimq=$trimqual

# Clean tmp directory
rm -rf "$tmp_dir"/metavir_"$bin"
