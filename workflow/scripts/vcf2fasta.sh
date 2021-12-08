#!/bin/bash
# VCF to fasta

# Source environment with all necessary software.
#module load anaconda; source activate snakemake
module load bcftools

# read from command line.
ref=$1
vcf=$2
sample_name=$3
bed=$4
fasta=$5

# Get sample name for correct genotype
samp=$(bcftools query -l ${vcf} )

# Make consensus masked sequence & rename fasta header. 
bcftools consensus --include 'TYPE!="indel"' --mask ${bed} --fasta-ref ${ref} --sample ${samp} --absent 'N' --missing 'N' ${vcf}  | \
    seqtk rename - ${sample_name}  > ${fasta}