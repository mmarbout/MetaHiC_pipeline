#! /bin/bash

# Set variables
bin=$1
bin_simple=$2
tmp_dir=$3
out_dir=$4
threads=$5
mem=$6

echo "MetaVIR bin $bin in progress..."

# Check tmp directory
rm -rf "$tmp_dir"/metavir_spades_"$bin"
mkdir -p "$tmp_dir"/metavir_spades_"$bin"

# Run spades
echo 'Run spades...'
spades.py \
    -o "$out_dir"/metavir_"$bin"/spades/ \
    -1 "$out_dir"/metavir_"$bin"/sg_cleaned_"$bin"_R1.fq.gz \
    -2 "$out_dir"/metavir_"$bin"/sg_cleaned_"$bin"_R2.fq.gz \
    -1 "$out_dir"/metavir_"$bin"/hc_cleaned_"$bin"_R1.fq.gz \
    -2 "$out_dir"/metavir_"$bin"/hc_cleaned_"$bin"_R2.fq.gz \
    --trusted-contigs "$out_dir"/metavir_"$bin"/"$bin_simple".fa \
    -t "$threads" \
    -m "$mem" \
    --only-assembler \
    --tmp-dir "$tmp_dir"/metavir_spades_"$bin"/

# Clean tmp directory
rm -rf "$tmp_dir"/metavir_spades_"$bin"
