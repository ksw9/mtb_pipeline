#!/bin/sh
# properties = {"type": "single", "rule": "vqsr", "local": false, "input": ["/labs/jandr/walter/tb/data/refs/H37Rv.fa", "process/vars/L5639_S9_bwa_H37Rv_gatk.vcf.gz"], "output": ["process/vars/L5639_S9_bwa_H37Rv_gatk_vqsr.vcf.gz"], "wildcards": {"samp": "L5639_S9", "mapper": "bwa", "ref": "H37Rv"}, "params": {}, "log": ["process/vars/L5639_S9_bwa_H37Rv_gatk_vqsr.log"], "threads": 1, "resources": {}, "jobid": 530, "cluster": {"time": 7200, "mem": "2G", "account": "jandr"}}
 cd /oak/stanford/scg/lab_jandr/walter/tb/capture && \
PATH='/scg/apps/software/snakemake/5.27.4/bin':$PATH /scg/apps/software/snakemake/5.27.4/bin/python3.9 \
-m snakemake process/vars/L5639_S9_bwa_H37Rv_gatk_vqsr.vcf.gz --snakefile /oak/stanford/scg/lab_jandr/walter/tb/capture/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /oak/stanford/scg/lab_jandr/walter/tb/capture/.snakemake/tmp.au998iez /labs/jandr/walter/tb/data/refs/H37Rv.fa process/vars/L5639_S9_bwa_H37Rv_gatk.vcf.gz --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler ilp \
\
\
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules vqsr --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /oak/stanford/scg/lab_jandr/walter/tb/capture/.snakemake/tmp.au998iez/530.jobfinished || (touch /oak/stanford/scg/lab_jandr/walter/tb/capture/.snakemake/tmp.au998iez/530.jobfailed; exit 1)

