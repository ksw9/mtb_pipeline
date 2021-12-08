#!/bin/sh
# properties = {"type": "single", "rule": "vqsr_filter", "local": false, "input": ["process/vars/L5005_S22_bwa_H37Rv_gatk_vqsr.vcf.gz"], "output": ["process/vars/L5005_S22_bwa_H37Rv_gatk_vqsr_snps.vcf.gz", "process/vars/L5005_S22_bwa_H37Rv_gatk_vqsrfilt.vcf.gz"], "wildcards": {"samp": "L5005_S22", "mapper": "bwa", "ref": "H37Rv"}, "params": {}, "log": ["process/vars/L5005_S22_bwa_H37Rv_gatk_vqsrfilt.log"], "threads": 1, "resources": {}, "jobid": 665, "cluster": {"time": 7200, "mem": "2G", "account": "jandr"}}
 cd /oak/stanford/scg/lab_jandr/walter/tb/capture && \
PATH='/scg/apps/software/snakemake/5.27.4/bin':$PATH /scg/apps/software/snakemake/5.27.4/bin/python3.9 \
-m snakemake process/vars/L5005_S22_bwa_H37Rv_gatk_vqsrfilt.vcf.gz --snakefile /oak/stanford/scg/lab_jandr/walter/tb/capture/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /oak/stanford/scg/lab_jandr/walter/tb/capture/.snakemake/tmp.au998iez process/vars/L5005_S22_bwa_H37Rv_gatk_vqsr.vcf.gz --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler ilp \
\
\
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules vqsr_filter --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /oak/stanford/scg/lab_jandr/walter/tb/capture/.snakemake/tmp.au998iez/665.jobfinished || (touch /oak/stanford/scg/lab_jandr/walter/tb/capture/.snakemake/tmp.au998iez/665.jobfailed; exit 1)

