#!/bin/sh
# properties = {"type": "single", "rule": "realign_indels", "local": false, "input": ["/labs/jandr/walter/tb/data/refs/MTB_ancestor_reference.fasta", "process/bams/5720_S13_bwa_MTB_ancestor_reference.merged.rmdup.bam"], "output": ["process/bams/5720_S13_bwa_MTB_ancestor_reference.merged.realn.bam"], "wildcards": {"samp": "5720_S13", "mapper": "bwa", "ref": "MTB_ancestor_reference"}, "params": {}, "log": ["process/bams/5720_S13_bwa_MTB_ancestor_reference_realn_indels.log"], "threads": 1, "resources": {}, "jobid": 655, "cluster": {"time": 7200, "mem": "2G", "account": "jandr"}}
cd /oak/stanford/scg/lab_jandr/walter/tb/capture && \
/home/kwalter/.conda/envs/snakemake/bin/python3.6 \
-m snakemake process/bams/5720_S13_bwa_MTB_ancestor_reference.merged.realn.bam --snakefile /oak/stanford/scg/lab_jandr/walter/tb/capture/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /oak/stanford/scg/lab_jandr/walter/tb/capture/.snakemake/tmp.9oxjlqnb /labs/jandr/walter/tb/data/refs/MTB_ancestor_reference.fasta process/bams/5720_S13_bwa_MTB_ancestor_reference.merged.rmdup.bam --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules realign_indels --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/oak/stanford/scg/lab_jandr/walter/tb/capture/.snakemake/tmp.9oxjlqnb/655.jobfinished" || (touch "/oak/stanford/scg/lab_jandr/walter/tb/capture/.snakemake/tmp.9oxjlqnb/655.jobfailed"; exit 1)

