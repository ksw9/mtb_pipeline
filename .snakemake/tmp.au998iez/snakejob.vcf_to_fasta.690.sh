#!/bin/sh
# properties = {"type": "single", "rule": "vcf_to_fasta", "local": false, "input": ["/labs/jandr/walter/tb/data/refs/H37Rv.fa", "process/vars/L5639_S9_bwa_H37Rv_gatk_qfilt.vcf.gz", "/labs/jandr/walter/varcal/data/refs/ppe_gagneux_0based_fmt.bed"], "output": ["process/fasta/L5639_S9_bwa_H37Rv_gatk_qfilt.fa", "process/fasta/L5639_S9_bwa_H37Rv_gatk_qfilt_ppemask.fa"], "wildcards": {"samp": "L5639_S9", "mapper": "bwa", "ref": "H37Rv", "filt": "qfilt"}, "params": {}, "log": ["process/fasta/L5639_S9_bwa_H37Rv_gatk_qfilt.log"], "threads": 1, "resources": {}, "jobid": 690, "cluster": {"time": 7200, "mem": "2G", "account": "jandr"}}
 cd /oak/stanford/scg/lab_jandr/walter/tb/capture && \
PATH='/scg/apps/software/snakemake/5.27.4/bin':$PATH /scg/apps/software/snakemake/5.27.4/bin/python3.9 \
-m snakemake process/fasta/L5639_S9_bwa_H37Rv_gatk_qfilt.fa --snakefile /oak/stanford/scg/lab_jandr/walter/tb/capture/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /oak/stanford/scg/lab_jandr/walter/tb/capture/.snakemake/tmp.au998iez /labs/jandr/walter/tb/data/refs/H37Rv.fa process/vars/L5639_S9_bwa_H37Rv_gatk_qfilt.vcf.gz /labs/jandr/walter/varcal/data/refs/ppe_gagneux_0based_fmt.bed --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler ilp \
\
\
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules vcf_to_fasta --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /oak/stanford/scg/lab_jandr/walter/tb/capture/.snakemake/tmp.au998iez/690.jobfinished || (touch /oak/stanford/scg/lab_jandr/walter/tb/capture/.snakemake/tmp.au998iez/690.jobfailed; exit 1)

