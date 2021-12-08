#!/bin/sh
# properties = {"type": "single", "rule": "gatk_call", "local": false, "input": ["/labs/jandr/walter/tb/data/refs/H37Rv.fa", "process/bams/8328-BC_S9_bwa_H37Rv.merged.rmdup.bam"], "output": ["process/vars/8328-BC_S9_bwa_H37Rv_gatk.vcf.gz"], "wildcards": {"samp": "8328-BC_S9", "mapper": "bwa", "ref": "H37Rv"}, "params": {"ploidy": "1"}, "log": ["process/vars/8328-BC_S9_bwa_H37Rv_gatk.log"], "threads": 1, "resources": {}, "jobid": 439, "cluster": {"time": 7200, "mem": "100G", "account": "jandr"}}
 cd /oak/stanford/scg/lab_jandr/walter/tb/capture && \
PATH='/scg/apps/software/snakemake/5.27.4/bin':$PATH /scg/apps/software/snakemake/5.27.4/bin/python3.9 \
-m snakemake process/vars/8328-BC_S9_bwa_H37Rv_gatk.vcf.gz --snakefile /oak/stanford/scg/lab_jandr/walter/tb/capture/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /oak/stanford/scg/lab_jandr/walter/tb/capture/.snakemake/tmp.au998iez /labs/jandr/walter/tb/data/refs/H37Rv.fa process/bams/8328-BC_S9_bwa_H37Rv.merged.rmdup.bam --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler ilp \
\
\
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules gatk_call --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /oak/stanford/scg/lab_jandr/walter/tb/capture/.snakemake/tmp.au998iez/439.jobfinished || (touch /oak/stanford/scg/lab_jandr/walter/tb/capture/.snakemake/tmp.au998iez/439.jobfailed; exit 1)

