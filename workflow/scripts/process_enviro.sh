##########################
##### Process Enviro #####
##########################

#### Set up work-space
qlogin -m 100 -c 12

# Set up environmemnt
ref=/labs/jandr/walter/tb/data/refs/H37Rv.fa
DATA_DIR=/labs/jandr/walter/tb/data/
PROCESS_DIR=/labs/jandr/walter/tb/enviro/
SNPEFF_DIR=/labs/jandr/walter/tb/data/refs/snpEff/
SCRIPTS_DIR=${PROCESS_DIR}/workflow/scripts/

module add anaconda/3_2022.05
module add IQ-TREE/2.2.0 
mamba activate snakemake

## Push to github for full pipeline.
git add -A 
git commit -m 'adding scripts'
git push -u origin master

##############################
#### Process capture runs #### old (run in mtb directory). updated Snakemake file outputs: mtb_tgen
##############################
cd ${PROCESS_DIR}

# Run Snakemake on the cluster. 
today=$(date +"%Y-%m-%d")
nohup snakemake -j 500 -k --cluster-config config/cluster_config.yaml  \
--cluster "sbatch -A {cluster.account} --mem={cluster.memory} -t {cluster.time} --cpus-per-task {threads} " > runs/snakemake_${today}.out &

# Look at run errors. 
runs/snakemake_${today}.out

nohup snakemake $file -j 500 -k --cluster-config config/config.yaml --cluster \
 'sbatch -A {cluster.account} --mem={cluster.mem} -t {cluster.time} --cpus-per-task {threads} --output "runs/slurm/slurm-%j.out"' > runs/snakemake_${today}.out &
 
# Try dry-run locally: 
snakemake -np process/vars/10561-Room_S8_bwa_MTB_ancestor_reference_gatk_vqsrfilt.vcf.gz

# Try on cluster for single sample
nohup snakemake results/8328-BC_S9/vars/8328-BC_S9_bwa_H37Rv_gatk_ann.vcf.gz -j 1 -k --cluster-config config/cluster_config.yaml  --cluster "sbatch -A {cluster.account} --mem={cluster.memory} -t {cluster.time} --cpus-per-task {threads} " > runs/snakemake_${today}.out &

# Test run. 
snakemake results/MT01_MtB_Baits-2021-10-06rerun2/4514-lighter-384/stats/4514-lighter-384_read_counts.txt  -j 1 -k --cluster-config config/cluster_config.yaml --cluster "sbatch -A {cluster.account} --mem={cluster.memory} -t {cluster.time} --cpus-per-task {threads} " > runs/snakemake_${today}.out &


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
for vcf in results/{4514,10561,8328,H37RV_S1}*/vars/*qfilt.vcf.gz ; do 
  vcf_filt=enviro/per_samp/$(basename ${vcf/.vcf.gz})_isnv.vcf.gz
  vcf_fixed_diffs=enviro/per_samp/$(basename ${vcf/.vcf.gz})_fixed.vcf.gz
  #bcftools filter --include 'AD[0:0] >= 5 & AD[0:1] >= 5 & TYPE != "indel"' $vcf | bgzip > ${vcf_filt}
  # Also get a vcf with fixed differences for each VCF. 
  bcftools view -T 'enviro/fixed_diffs.txt' ${vcf} | bgzip > ${vcf_fixed_diffs}
done

# Read this into R. /labs/jandr/walter/tb/mtb/scripts/plot_coverage.Rmd

# Collate all tb-profiler runs (directory needs to be named results)
cp results/*/stats/results/*json results/stats/results/
cd results/stats/
tb-profiler collate 

# Also collect json files from IS-1000 run
cp  ../mtb_tgen/results/IS-1000/*/stats/results/*json results/stats/results/
cd results/stats/
tb-profiler collate 

# Copy environmental plots locally. 
for file in process/plots/* ; do 
  rclone copy $file  box:Box/TB/Environmental/plots/
done

## Filter enviro VCFs to exclude MAF >5% and < 95%
for vcf in results/{4514,10561,8328,H37RV_S1}*/vars/*qfilt.vcf.gz ; do 
 echo $vcf
 vcf_filt=enviro/per_samp/$(basename ${vcf/.vcf.gz})_nohets.vcf.gz
# bcftools filter --exclude 'AD[0:0]/FORMAT/DP > 0.1 & AD[0:0]/FORMAT/DP < 0.9' --set-GTs . $vcf | bgzip > ${vcf_filt}
 bcftools filter -e 'AD[0:0]/(FORMAT/DP) > .05 && AD[0:1]/(FORMAT/DP) > 0.05' --set-GTs . $vcf | bgzip > ${vcf_filt}

done

# VCF to fasta and get pairwise differences. 
bed='/labs/jandr/walter/varcal/data/refs/ppe_gagneux_0based_snpeff.bed.gz'
for vcf in enviro/per_samp/*_nohets.vcf.gz; do
  tabix $vcf
  sample_name=$(basename ${vcf/_bwa_H37Rv_gatk_qfilt_fixed.vcf.gz})
  echo $sample_name
  fasta=enviro/per_samp/${sample_name}_fixed.fa
  ${SCRIPTS_DIR}vcf2fasta.sh $ref $vcf $sample_name $bed $fasta
done

## Get snps fasta and VCF
for samp in enviro/per_samp/*_bwa_H37Rv_gatk_qfilt_nohets.vcf.gz_fixed.fa; do
  mv $samp ${samp/_bwa_H37Rv_gatk_qfilt_nohets.vcf.gz_fixed.fa/_nohets.fa}
done

cat $ref enviro/per_samp/*nohets.fa > enviro/enviro_fixed_msa.fa
snp-sites -m enviro/enviro_fixed_msa.fa -o enviro/enviro_fixed_msa_snps.fa

## Compare to unmasked fastas.
cat results/{4514,10561,8328,H37RV_S1}*/fasta/*qfilt.fa $ref > enviro/enviro_msa.fa
snp-sites -m enviro/enviro_msa.fa -o enviro/enviro_msa_snps.fa

# Look at variants in H37Rv. 
cat $ref ../results/H37RV_S1/fasta/H37RV_S1_bwa_H37Rv_gatk_qfilt.fa > pos_control.fa
snp-sites pos_control.fa 

## Create MSA of 7 environmental samples + 1 culture sample. Add culture sample of additional Rif+ case.
cat $ref results/{4514,10561,8328}*/fasta/*qfilt.fa ../mtb_tgen/results/Stanford/11080-7651-envculture/fasta/11080-7651-envculture_bwa_H37Rv_gatk_qfilt.fa \
> enviro/enviro_culture_msa.fa

## Get snps fasta and VCF
snp-sites -m enviro/enviro_culture_msa.fa -o enviro/enviro_culture_msa_snps.fa

######################
#### Enviro Trees ####
######################
#### Enviro Trees ####
cd /labs/jandr/walter/tb/mtb/

# A. Create an MSA with all environmental sequences. ## Use updated RAXML-ng!
cat enviro/enviro_fixed_msa.fa ../mtb_tgen/results/IS-1000/*/fasta/*qfilt.fa > enviro/stanford_msa.fa 
snp-sites -m enviro/stanford_msa.fa -o enviro/stanford_msa_snps.fa

#### Updating RaxML jobs here ####
cd /labs/jandr/walter/tb/mtb/enviro/
ref=/labs/jandr/walter/varcal/data/refs/H37Rv.fa
msa=stanford_msa_snps.fa
prefix=trees/stanford_is-1000_workers
threads=8
pin=thread-nopin
brlen=mod
model=GTR
/labs/jandr/walter/tb/scripts/raxml_ng_parallel.sh \
$msa $ref $model $prefix $threads $pin $brlen
# Submit job. 
sbatch --account=jandr -t 48:00:00 /labs/jandr/walter/tb/scripts/raxml_ng_parallel.sh \
$msa $ref $model $prefix $threads $pin $brlen -o trees/stanford_is-1000.out \
 -e trees/stanford_is-1000.er
 

# Testing now with increased # of workers.
# the "local" run is the fastest. 

# Extract reads that map to rpoB point mutation. 
bam=results/10561-Food_S7/bams/10561-Food_S7_bwa_H37Rv.merged.rmdup.bam

samtools view -h ${bam} NC_000962.3:761155 | samtools fasta - # Errors with IO

# files are blocked by large tarball:
/labs/jandr/walter/tb/samps.tar.tgz

########################
#### Enviro upload ####
#######################
cd $PROCESS_DIR
module load aspera
# Copy to SRA
key=/labs/jandr/walter/lambda/metadata/aspera.openssh
file_path=/tmp/kwalter/sra_upload
mkdir -p $file_path

# Copy BAM files of interest to filepath
cp results/{10561,4514,8328}*/bams/*rmdup.bam ${file_path}

# Rename incorrectly named sample.
cd ${file_path}
mv 8328-Food_S10_bwa_H37Rv.merged.rmdup.bam 4514-Food_S10_bwa_H37Rv.merged.rmdup.bam
ascp -i ${key} -QT -l100m -k1 -d ${file_path} subasp@upload.ncbi.nlm.nih.gov:uploads/kwalter_stanford.edu_hG59bCBT

