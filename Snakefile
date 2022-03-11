import os
import sys
import glob
import pandas as pd
import plotly as plt
import plotly.graph_objects as go
from snakemake.utils import min_version

min_version("5.11.0")

include: "scripts/common.py"
configfile: "config/config.yaml"

st = pd.read_table('config/samplesheet.tsv').set_index('sample')
sample_file_dict = gather_files()

treated_samples = list(st.index)
control_samples = list(st["control"])
all_samples = list(set(treated_samples + control_samples))

rule all:
	input:
		# quality control -------------------------------------------------------------------------
		# expand("data/fastqc/{read}_fastqc.html", read = all_samples),
		# expand("data/fastq_screen/{read}_screen.txt", read = reads),
		expand("data/preseq/estimates_{sample}.txt", sample = all_samples),
		expand("data/preseq/lcextrap_{sample}", sample = all_samples),
		# read alignment --------------------------------------------------------------------------
		expand("data/bowtie2/{sample}.bam", sample = all_samples),
		expand("data/ban/{sample}.ban.sorted.markd.bam", sample = all_samples),
		expand("data/tracks/{sample}.bw", sample = all_samples),
		# peak calling ----------------------------------------------------------------------------
		expand("data/macs2/{sample}_peaks.xls", sample = treated_samples)

rule fastqc:
	input:
		"data/raw/{read}.fastq.gz"
	output:
		"data/fastqc/{read}_fastqc.html"
	conda:
		"../envs/fastqc.yaml"
	log:
		"data/logs/fastqc_{read}.log"
	threads: 4
	shell:
		"fastqc -t {threads} --outdir data/fastqc {input} > {log} 2>&1"

rule fastq_screen:
	input:
		fastq = "data/fastp/{read}.fastq.gz",
		config = config["FASTQ_SCREEN_CONFIG"]
	output:
		"data/fastq_screen/{read}_screen.txt"
	conda:
		"../envs/fastq_screen.yaml"
	log:
		"data/logs/fastq_screen_{read}.txt"
	threads: 8
	shell:
		"fastq_screen --aligner bowtie2 --threads {threads} --outdir data/fastq_screen "
		"--conf {input.config} --force {input.fastq} > {log} 2>&1"

# rule single_fastp:
# 	input:
# 		"data/raw/{sample}.fastq.gz"
# 	output:
# 		"data/fastp/{sample}.fastq.gz"
# 	conda:
# 		"envs/fastp.yaml"
# 	threads: 4
# 	shell:
# 		"fastp -i {input} -o {output} "
# 		"--detect_adapter_for_pe --thread {threads} -j {log} -h /dev/null"

# rule paired_fastp:
# 	input:
# 		r1 = "data/raw/{sample}_R1.fastq.gz",
# 		r2 = "data/raw/{sample}_R2.fastq.gz"
# 	output:
# 		r1 = "data/fastp/{sample}_R1.fastq.gz",
# 		r2 = "data/fastp/{sample}_R2.fastq.gz"
# 	conda:
# 		"envs/fastp.yaml"
# 	log:
# 		"data/logs/{sample}.fastp.json"
# 	threads: 4
# 	shell:
# 		"fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
# 		"--detect_adapter_for_pe --thread {threads} -j {log} -h /dev/null"

rule single_bowtie2:
	input:
		"data/raw/{sample}.fastq.gz"
	output:
		"data/bowtie2/{sample}.bam"
	log:
		err="data/logs/single_bowtie2_{sample}.err"
	conda:
		"envs/bowtie2.yaml"
	threads: 8
	shell:
		"bowtie2 --local -p {threads} -x {config[GENOME]} "
		"-U {input} 2>{log.err} | samtools view -@ {threads} -Sbh - > {output}"

rule paired_bowtie2:
	input:
		r1 = "data/raw/{sample}_R1.fastq.gz",
		r2 = "data/raw/{sample}_R2.fastq.gz"
	output:
		"data/bowtie2/{sample}.bam"
	log:
		err="data/logs/paired_bowtie2_{sample}.err"
	conda:
		"envs/bowtie2.yaml"
	threads: 8
	shell:
		"bowtie2 --local --very-sensitive-local "
		"--no-unal --no-mixed --threads {threads} "
		"-I 10 -X 700 -x {config[GENOME]} "
		"-1 {input.r1} -2 {input.r2} 2>{log.err} | samtools view -@ {threads} -Sbh - > {output}"

rule sortbam:
	input:
		"data/bowtie2/{sample}.bam"
	output:
		"data/bowtie2/{sample}.sorted.bam"
	conda:
		"envs/bowtie2.yaml"
	threads: 4
	shell:
		"samtools sort -@ {threads} -o {output} {input}"

rule markd:
	input:
		"data/bowtie2/{sample}.sorted.bam"
	output:
		"data/markd/{sample}.sorted.markd.bam"
	conda:
		"envs/sambamba.yaml"
	threads: 4
	shell:
		"sambamba markdup -r --tmpdir data/markd {input} {output}"

rule banlist:
	input:
		"data/markd/{sample}.sorted.markd.bam"
	output:
		"data/ban/{sample}.ban.sorted.markd.bam"
	conda:
		"envs/bedtools.yaml"
	shell:
		"bedtools intersect -v -a {input} -b {config[BANLIST]} > {output}"

rule index:
	input:
		"data/ban/{sample}.ban.sorted.markd.bam"
	output:
		"data/ban/{sample}.ban.sorted.markd.bam.bai"
	conda:
		"envs/bowtie2.yaml"
	threads: 8
	shell:
		"samtools index -@ {threads} {input}"

rule tracks:
	input:
		"data/ban/{sample}.ban.sorted.markd.bam",
		"data/ban/{sample}.ban.sorted.markd.bam.bai"
	output:
		"data/tracks/{sample}.bw"
	conda:
		"envs/deeptools.yaml"
	threads: 8
	shell:
		"bamCoverage -b {input[0]} -o {output} -p {threads} --binSize 10 --smoothLength 50 --normalizeUsing CPM"

rule preseq:
	input:
		"data/ban/{sample}.ban.sorted.markd.bam"
	output:
		"data/preseq/estimates_{sample}.txt"
	conda:
		"envs/preseq.yaml"
	resources:
		defect_mode = defect_mode
	log:
		"data/logs/preseq_{sample}.log"
	shell:
		"preseq c_curve -B {resources.defect_mode} -l 1000000000 -P -o {output} {input} > {log} 2>&1"

rule preseq_lcextrap:
	input:
		"data/ban/{sample}.ban.sorted.markd.bam"
	output:
		"data/preseq/lcextrap_{sample}"
	conda:
		"envs/preseq.yaml"
	resources:
		defect_mode = defect_mode
	log:
		"data/logs/preseq_lcextrap_{sample}.log"
	shell:
		"preseq lc_extrap -B {resources.defect_mode} -l 1000000000 -P -e 1000000000 -o {output} {input} > {log} 2>&1"

rule macs2:
	input:
		sample = "data/ban/{sample}.ban.sorted.markd.bam"
	output:
		"data/macs2/{sample}_peaks.xls"
	params:
		control = get_control,
		peak = get_peak,
		seed = 0
	conda:
		"envs/macs2.yaml"
	shell:
		"macs2 callpeak -t {input} {params.control} -n {wildcards.sample} "
		"--outdir data/macs2 --seed {params.seed}"

