# Mtb variant calling 

Pipeline for *M. tuberculosis* variant identification from short-read data.

## Usage
```
# Navigate to root directory. 
# Identify files that have not been produced (no run)
snakemake -np

# Run snakemake
snakemake
```

## Directory structure

```
├── .gitignore
├── README.md
├── LICENSE.md
├── workflow
│   ├── rules
│   ├── envs
|   │   ├── tool1.yaml
│   ├── scripts
|   │   ├──process_stanford_tb.sh
|   │   ├──trim_reads.sh
|   │   ├──run_kraken.sh
|   │   ├──map_reads.sh
|   │   ├──cov_stats.sh
|   │   ├──mykrobe_predict.sh
|   │   ├──call_varsk.sh
|   │   ├──vqsr.sh
|   │   ├──vcf2fasta.sh
│   ├── notebooks
│   ├── report
|   └── Snakefile
├── config
│   ├── config.yaml
├── results
│   ├── trim
│   ├── bams
│   ├── vars
│   ├── stats
│   ├── fasta
└── resources
```