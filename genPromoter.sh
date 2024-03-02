#!/bin/bash

set -o pipefail

promoter_length=1000

echo "Enter species name:"
read species_name

fasta_file="Caenorhabditis_Genome_Sequences/${species_name}.fa"
gff_file="Caenorhabditis_Gene_Annotations/${species_name}.gff3"

genePattern="g1039|g106|g1096|g11140|g12655|g12657|g12736|g13088|g13525|g1652|g166|g212|g269|g274|g5238|g543|g568|g623|g634|g6886|g6889|g6892"

echo "Gene Pattern: $genePattern"

# Extract information for genes from input gff file
awk -F'\t' -v pattern="$genePattern" '$3 == "gene" && !($9 ~ "Name=" pattern || $9 ~ "ID=" pattern || $9 ~ "sequence_name=" pattern)' "$gff_file" > "${species_name}AllGenes.gff"

sortBed -i "${species_name}AllGenes.gff" | gff2bed > "${species_name}AllGenes.bed"

# Generate index for fasta file
samtools faidx "$fasta_file"

# Create a table (from index file) that contains contig/chromosome sizes
cut -f1-2 "${fasta_file}.fai" > "${species_name}sizes.chr"

# Create bed file that contains locations of promoters
bedtools flank -i "${species_name}AllGenes.bed" -g "${species_name}sizes.chr" -l $promoter_length -r 0 -s > "${species_name}AllUpstreamRegions.bed"

# Extract promoter regions from fasta file with contigs/chromosomes using bed file
bedtools getfasta -s -fi "$fasta_file" -bed "${species_name}AllUpstreamRegions.bed" -fo "${species_name}AllUpstreamRegions.fa" -name

echo "Results are in ${species_name}AllUpstreamRegions.fa file"

# Count promoters
result=$(grep -c '^>' "${species_name}AllUpstreamRegions.fa")

echo "Number of obtained promoters: $result"
