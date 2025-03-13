rule all:
    input:
        expand("010.fastqc/{sample}_fastqc.zip", sample=glob_wildcards("000.fastq/{sample}.fastq").sample)

rule fastqc:
    input:
        "000.fastq/{sample}.fastq"
    output:
        "010.fastqc/{sample}_fastqc.zip"
    shell:
        "mkdir -p 010.fastqc && fastqc {input} -o 010.fastqc/"

