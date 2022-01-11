#########################################
##### Process Stanford Mtb seq data #####
#########################################

#### Set up work-space
qlogin -m 100 -c 12

# Set up environmemnt
ref=/labs/jandr/walter/tb/data/refs/H37Rv.fa
DATA_DIR=/labs/jandr/walter/tb/data/Stanford/
PROCESS_DIR=/labs/jandr/walter/tb/mtb/
SNPEFF_DIR=/labs/jandr/walter/tb/data/refs/snpEff/
module load bcftools 
module load anaconda
source activate snakemake

## Push to github
git add -A 
git commit -m 'adding scripts'
git push/git pull

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
$HOME/bin/bs download project --name=MT03_MTB_2021-11-24 -o MT03_MTB_2021-11-24 
$HOME/bin/bs download project --name=MT04_A10-12_MTB_2021-12-10 -o MT04_A10-12_MTB_2021-12-10

# List samples within a single project
$HOME/bin/bs list biosample --project-name=MT01_MtB_Baits-2021-10-06rerun2 

# Move FASTQs out of subdirectories for easier processing
for dir in *; do 
 echo $dir
 mv ${dir}/*/*fastq.gz ${dir}/
done

##############################
#### Process capture runs ####
##############################
cd ${PROCESS_DIR}

# Run Snakemake on the cluster. 
today=$(date +"%Y-%m-%d")
nohup snakemake -j 500 -k --cluster-config config/config.yaml --cluster "sbatch -A {cluster.account} --mem={cluster.mem} -t {cluster.time} --cpus-per-task {threads} -o {cluster.output}" > runs/snakemake_${today}.out & 

# Look at run errors. 
runs/snakemake_${today}.out

nohup snakemake $file -j 500 -k --cluster-config config/config.yaml --cluster \
 'sbatch -A {cluster.account} --mem={cluster.mem} -t {cluster.time} --cpus-per-task {threads} --output "runs/slurm/slurm-%j.out"' > runs/snakemake_${today}.out &
 
# Try dry-run locally: 
snakemake -np process/vars/10561-Room_S8_bwa_MTB_ancestor_reference_gatk_vqsrfilt.vcf.gz

# Try on cluster for single sample
nohup snakemake results/bams/4514-food_S4_bwa_H37Rv.merged.rmdup.bam -j 10 -k --cluster-config cluster.yaml --cluster "sbatch -A {cluster.account} --mem={cluster.mem} -t {cluster.time} --cpus-per-task {threads} " > runs/snakemake_${today}.out &

# snakemake question
# https://stackoverflow.com/questions/69578275/snakemake-combine-inputs-over-one-wildcard

# Dry run of a single rule. 
snakemake -R per_samp_stats -np

# Combine coverage thus far. 
sed -n '1'p  $(ls process/stats/*combined_stats.csv | head -n1) > combined_cov_stats.csv
for file in *process/stats/*combined_stats.csv ; 
  do sed -n '2'p $file 
done >> combined_cov_stats.csv

# Copy coverage stats locally. 
rclone copy combined_cov_stats.csv box:Box/TB/seq_data/cov/ 

# Buid DAG of jobs *need to comment out print command; this causes issue with dot command.
snakemake --dag  results/fasta/4514-*fa | dot -Tsvg > process/plots/snakemake_dag.svg

###############################
#### Environmental samples ####
###############################
cd ${PROCESS_DIR}

## Create MSA of 7 environmental samples + 1 culture sample
cat $ref results/{4514,10561,8328}*/fasta/*qfilt.fa > enviro/enviro_msa.fa

## Get snps fasta and VCF
snp-sites -m enviro/enviro_msa.fa -o enviro/enviro_msa_snps.fa
snp-sites -v enviro/enviro_msa.fa -o enviro/enviro.vcf

# Replace chromosome name and zip
cat enviro/enviro.vcf | awk '{gsub(/^1/,"NC_000962.3"); print}' | bgzip > enviro/enviro.vcf.gz
tabix -f enviro/enviro.vcf.gz

# Get sites of fixed differences
bcftools query -f 'NC_000962.3\t%POS\n' enviro/enviro.vcf.gz > enviro/fixed_diffs.txt

# For all environmental VCFs & H37Rv, get het sites and sites included in the SNP positions, output as iSNV VCF. 
for vcf in results/{4514,10561,8328, H37RV_S1}*/vars/*qfilt.vcf.gz ; do 
  vcf_filt=enviro/per_samp/$(basename ${vcf/.vcf.gz})_isnv.vcf.gz
  vcf_fixed_diffs=enviro/per_samp/$(basename ${vcf/.vcf.gz})_fixed.vcf.gz
  bcftools filter --include 'AD[0:0] >= 5 & AD[0:1] >= 5 & TYPE != "indel"' $vcf | bgzip > ${vcf_filt}
  # Also get a vcf with fixed differences for each VCF. 
  bcftools view -T 'enviro/fixed_diffs.txt' ${vcf} | bgzip > ${vcf_fixed_diffs}
done

# Read this into R. /labs/jandr/walter/tb/mtb/scripts/plot_coverage.Rmd

###########################
#### Process TGen runs ####
###########################
DATA_DIR=/labs/jandr/walter/tb/data/

sample_list=config/tgen_samples.tsv
echo 'sample,fastq_1,fastq_2','batch','run'> ${sample_list}

for file in ${DATA_DIR}IS-1000/Walter-TB/*R1_001.fastq.gz; do
  fq1=${file}
  fq2=${file/R1_001.fastq.gz/R2_001.fastq.gz}
  samp=$(basename $file)
  samp=${samp/_S*_R1_001.fastq.gz}
  dir=$(dirname $file)
  run=${dir##*/}
  bat=$(dirname $(echo ${dir/"${DATA_DIR}"}))
  echo ${samp},${fq1},${fq2},${bat},${run} >> ${sample_list}
done

# Run Snakemake on the cluster. 
nohup snakemake -j 500 -k --cluster-config config/cluster_config.yaml --cluster "sbatch -A {cluster.account} --mem={cluster.mem} -t {cluster.time} --cpus-per-task {threads} " > runs/snakemake_${today}.out & 

# Update SLURM profile here: https://eriqande.github.io/eca-bioinf-handbook/managing-workflows-with-snakemake.html#using-snakemake-on-a-computing-cluster
######################
#### snpEff set-up ####
#######################

module load snpeff

snpEff="/scg/apps/software/snpeff/4.3t/snpEff/snpEff.jar"
bed='/labs/jandr/walter/varcal/data/refs/ppe_gagneux_0based_fmt.bed'
snpeff_bed='/labs/jandr/walter/varcal/data/refs/ppe_gagneux_0based_snpeff.bed'
gff=/labs/jandr/walter/tb/data/refs/H37Rv.gff.gz

# Download TB database locally : java -jar -Xmx50g ${snpEff} download Mycobacterium_tuberculosis_h37rv -dataDir /labs/jandr/walter/tb/data/refs/snpEff/
# config file is updated with snpEff datadir.
java -jar -Xmx50g ${snpEff} download Mycobacterium_tuberculosis_h37rv -dataDir /labs/jandr/walter/tb/data/refs/snpEff/

# Create a header line for adding PPE annotation
echo -e '##FORMAT=<ID=PPE,Number=1,Type=String,Description="Located within PE/PPE genes">' > /labs/jandr/walter/tb/data/refs/ppe_hdr.txt

# Is the bed file annoting correctly? Confirmed bed file is 0/1-based.  
# Chromosome      33581   33794   Rv0031 # From Bed file

### From gff (1-based)
# zcat $gff | grep Rv0031
# Chromosome	ena	gene	33582	33794	.	+	.	ID=gene:Rv0031;biotype=protein_coding;description=Possible remnant of a transposase;gene_id=Rv0031;logic_name=ena
# Chromosome	ena	mRNA	33582	33794	.	+	.	ID=transcript:CCP42753;Parent=gene:Rv0031;biotype=protein_coding;transcript_id=CCP42753
#     
