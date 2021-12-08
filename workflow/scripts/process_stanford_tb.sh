#########################################
##### Process Stanford Mtb seq data #####
#########################################
ssh -X smsh11dsu-srcf-d15-35.scg.stanford.edu

#### Set up work-space
qlogin -m 50 -c 12

# Set up environmemnt
ref=/labs/jandr/walter/tb/data/refs/H37Rv.fa
SCRIPTS_DIR=/labs/jandr/walter/tb/scripts/
DATA_DIR=/labs/jandr/walter/tb/data/Stanford/
cd /labs/jandr/walter/tb
module load bcftools 

module load anaconda
source activate snakemake
cd /labs/jandr/walter/tb/capture/

######################################
#### Download data from BaseSpace ####
######################################

cd ${DATA_DIR}

# List completed runs
$HOME/bin/bs list project

# Download runs
$HOME/bin/bs download project --name=MT01_MtB_Baits-2021-10-01rerun2 -o MT01_MtB_Baits-2021-10-01rerun2
$HOME/bin/bs download project --name=MT01_MtB_Baits-2021-10-06rerun2 -o MT01_MtB_Baits-2021-10-06rerun2
$HOME/bin/bs download project --name=MT01_MtB_Baits-2021-09-17 -o MT01_MtB_Baits-2021-09-17   
$HOME/bin/bs download project --name=MT02_MTB_2021-10-29 -o MT02_MTB_2021-10-29

$HOME/bin/bs list biosample --project-name=MT01_MtB_Baits-2021-10-06rerun2 


# Move FASTQs out of subdirectories for easier processing
for dir in *; do 
 echo $dir
 mv ${dir}/*/*fastq.gz ${dir}/
done

##############################
#### Process capture runs ####
##############################

# Run snakemake on these samples; combine fastqs for each. 

# Run pipeline. 
snakemake --cores 8 process/bams/MT01_MtB_Baits-2021-09-17_7056_S19_bwa_MTB_ancestor_reference_rg_sorted.bam
snakemake -np

# Run Snakemake on the cluster. 
today=$(date +"%Y-%m-%d")
nohup snakemake -j 500 -k --cluster-config cluster.yaml --cluster "sbatch -A {cluster.account} --mem={cluster.mem} -t {cluster.time} --cpus-per-task {threads} " > runs/snakemake_${today}.out &

# Errors: 
runs/snakemake_${today}.out

# Try dry-run locally: 
snakemake -np process/vars/10561-Room_S8_bwa_MTB_ancestor_reference_gatk_vqsrfilt.vcf.gz
# Try for single sample, vcf
snakemake -np process/vars/10561-Room_S8_bwa_MTB_ancestor_reference_gatk.vcf.gz

# Try on cluster for single sample
nohup snakemake process/fasta/4514-food_S4_bwa_H37Rv_gatk_qfilt.fa -j 10 -k --cluster-config cluster.yaml --cluster "sbatch -A {cluster.account} --mem={cluster.mem} -t {cluster.time} --cpus-per-task {threads} " > runs/snakemake_${today}.out &

nohup snakemake process/fasta/10561-Room_S8_bwa_H37Rv_gatk_qfilt.fa -j 8 -k --cluster-config cluster.yaml --cluster "sbatch -A {cluster.account} --mem={cluster.mem} -t {cluster.time} --cpus-per-task {threads} " > runs/snakemake_${today}.out &

# snakemake question
# https://stackoverflow.com/questions/69578275/snakemake-combine-inputs-over-one-wildcard

snakemake -R combine_stats -np

# Combine coverage thus far. 
sed -n '1'p  $(ls process/stats/*combined_stats.csv | head -n1) > combined_cov_stats.csv
for file in *process/stats/*combined_stats.csv ; 
  do sed -n '2'p $file 
done >> combined_cov_stats.csv

# Copy coverage stats locally. 
rclone copy combined_cov_stats.csv box:Box/TB/seq_data/cov/ 

# Add mosdepth coverage plot. 
for bam in process/bams/*.merged.rmdup.bam; do
  prefix=${bam/.merged.rmdup.bam}
  mosdepth --by 2000 -Q 30 ${prefix} ${bam}
done


###############################
#### Environmental samples ####
###############################

## Create MSA of 7 environmental samples
cat $ref process/fasta/{4514,10561,8328}*qfilt.fa > results/enviro_msa.fa

## Get snps fasta. 
snp-sites -m results/enviro_msa.fa -o results/enviro_msa_snps.fa


# Do the same for vqsr snps
cat $ref process/fasta/{4514,10561,8328}*vqsrfilt.fa > results/enviro_vqsr_msa.fa
snp-sites -m results/enviro_vqsr_msa.fa -oresults/enviro_vqsr_msa_snps.fa


# Check fastas are same length
for fasta in process/fasta/{4514,10561,8328}*qfilt.fa; do 
  seqtk comp $fasta
  grep -v '>' $fasta | wc -c
done
# true length: 4411532

# Substantial differences between sequences from same room. 
cat $ref process/fasta/{4514,10561,8328}*vqsrfilt.fa > results/enviro_msa_vqsr.fa

## Get snps fasta. 
snp-sites -m results/enviro_msa_vqsr.fa -o results/enviro_msa_vqsr_snps.fa

# Get positions of difference
## Get positions of differences for 10561-Food_S71 and 10561-BC_S61 
cat process/fasta/4514*qfilt.fa > results/4514_msa.fa
snp-sites -v results/4514_msa.fa -o a

#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  4514-food_S41   4514-lighter-384_S51
1       146985  .       T       C       .       .       .       GT      0       1
1       200746  .       T       G       .       .       .       GT      0       1
1       384253  .       T       C       .       .       .       GT      0       1
1       537019  .       A       G       .       .       .       GT      0       1
1       1570727 .       C       T       .       .       .       GT      0       1
1       1698774 .       T       C       .       .       .       GT      0       1
1       1944482 .       C       A       .       .       .       GT      0       1
1       2266487 .       C       G       .       .       .       GT      0       1
1       2266508 .       T       A       .       .       .       GT      0       1
1       2266517 .       C       T       .       .       .       GT      0       1
1       2266550 .       T       G       .       .       .       GT      0       1
1       2266553 .       G       C       .       .       .       GT      0       1
1       2283994 .       T       C       .       .       .       GT      0       1
1       2389429 .       C       T       .       .       .       GT      0       1
1       3796462 .       T       C       .       .       .       GT      0       1
1       4331685 .       G       A       .       .       .       GT      0       1
1       4408898 .       C       G       .       .       .       GT      0       1
(END)


# Test if pairwise diffs are in PPE genes
bedtools intsersect -a results/10561_snps.vcf -b $bed # No.

# Then test where the pairwise SNPs occur.

# Look at het sites in H37Rv
vcf_ref=process/vars/H37RV_S1_bwa_H37Rv_gatk_qfilt.vcf.gz
bcftools filter --include 'AD[0:0] >= 10 & AD[0:1] >= 10 & TYPE != "indel"' $vcf_ref # No het sites.




#### Look into duplicates and other sequencing errors. ####

# Try for a single bam file from the culture samples. 

# Look at total reads output 
data_dir=/labs/jandr/walter/tb/data/Stanford/MT02_MTB_2021-10-29/
for fastq in ${data_dir}/*.fastq.gz ; do 
  reads=$(zcat $fastq | wc -l); printf "${fastq}\t${reads}\n";
done > MT02_rawreads.tsv


