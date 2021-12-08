#!/bin/sh
# properties = {"type": "single", "rule": "gatk_call", "local": false, "input": ["/labs/jandr/walter/tb/data/refs/MTB_ancestor_reference.fasta", "process/bams/10561-Room_S8_bwa_MTB_ancestor_reference.merged.realn.bam"], "output": ["process/vars/10561-Room_S8_bwa_MTB_ancestor_reference_gatk.vcf.gz"], "wildcards": {"samp": "10561-Room_S8", "mapper": "bwa", "ref": "MTB_ancestor_reference"}, "params": {"ploidy": "1"}, "log": ["process/vars/10561-Room_S8_bwa_MTB_ancestor_reference_gatk.log"], "threads": 1, "resources": {}, "jobid": 589, "cluster": {"time": 7200, "mem": "2G", "account": "jandr"}}
cd /oak/stanford/scg/lab_jandr/walter/tb/capture && \
/home/kwalter/.conda/envs/snakemake/bin/python3.6 \
-m snakemake process/vars/10561-Room_S8_bwa_MTB_ancestor_reference_gatk.vcf.gz --snakefile /oak/stanford/scg/lab_jandr/walter/tb/capture/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /oak/stanford/scg/lab_jandr/walter/tb/capture/.snakemake/tmp.9oxjlqnb /labs/jandr/walter/tb/data/refs/MTB_ancestor_reference.fasta process/bams/10561-Room_S8_bwa_MTB_ancestor_reference.merged.realn.bam --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules gatk_call --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/oak/stanford/scg/lab_jandr/walter/tb/capture/.snakemake/tmp.9oxjlqnb/589.jobfinished" || (touch "/oak/stanford/scg/lab_jandr/walter/tb/capture/.snakemake/tmp.9oxjlqnb/589.jobfailed"; exit 1)

