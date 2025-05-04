samples = glob_wildcards("000.fastq/{sample}.GRCh38DH.exome.chr21.fq.gz").sample

rule all:
    input:
        "high_impact.tsv",
        "annotated_snps.vcf",
        "snps.sqlite",
        "snakefile",
        "upload_done.txt"

rule fastqc:
    input:
        expand("000.fastq/{sample}.GRCh38DH.exome.chr21.fq.gz", sample=samples)
    output:
        expand("fastqc_output/{sample}_fastqc.html", sample=samples)
    shell:
        "fastqc 000.fastq/*.fq.gz -o fastqc_output"

rule align:
    input:
        ref="/staging/leuven/stg_00079/teaching/hg38_9/chr9.fa",
        fq="000.fastq/{sample}.GRCh38DH.exome.chr21.fq.gz"
    output:
        bam="{sample}.bam"
    shell:
        "bwa mem {input.ref} {input.fq} | samtools sort -o {output.bam}"

rule index_bam:
    input:
        bam="{sample}.bam"
    output:
        bai="{sample}.bam.bai"
    shell:
        "samtools index {input.bam}"

rule call_variants:
    input:
        ref="/staging/leuven/stg_00079/teaching/hg38_9/chr9.fa",
        bamN="HG00101.bam", 
        bamT="HG00111.bam"
    output:
        "raw_snps.vcf"
    shell:
        "bcftools mpileup -Ou -f {input.ref} {input.bamN} {input.bamT} | \
         bcftools call -mv -Ov -o {output}"

rule clean_vcf:
    input:
        vcf="raw_snps.vcf",
        ref="/staging/leuven/stg_00079/teaching/hg38_9/chr9.fa"
    output:
        "clean_snps.vcf"
    shell:
        "vt decompose {input.vcf} | \
         vt normalize -n -r {input.ref} - | \
         vt uniq - | \
         vt view -f \"QUAL>20\" -h - > {output} "

rule annotate:
    input:
        vcf="clean_snps.vcf"
    output:
        "annotated_snps.vcf"
    params:
        snpeff_jar="/lustre1/project/stg_00079/teaching/I0U19a_conda_2025/share/snpeff-5.2-1/snpEff.jar",
        db="hg38",
        data_dir="/staging/leuven/stg_00079/teaching/snpeff_db"
    shell:
        "java -Xmx3400m -jar {params.snpeff_jar} eff {params.db} -dataDir {params.data_dir} \
         {input.vcf} > {output}"

rule download:
    output:
        directory("000.fastq")
    shell:
       """
        mkdir -p 000.fastq
        iget -r /gbiomed/home/large_omics_course/fastq/group_1 ./
        mv group_1/*.fq.gz 000.fastq/
        rmdir group_1
        """

rule upload:
    input:
        vcf="annotated_snps.vcf",
        db="snps.sqlite",
        snakefile="snakefile"
    output:
        touch("upload_done.txt")
    shell:
        """
        imkdir -p /gbiomed/home/large_omics_course/output/r1032561/ || true
        iput -f {input.vcf} /gbiomed/home/large_omics_course/output/r1032561/
        iput -f {input.db} /gbiomed/home/large_omics_course/output/r1032561/
        iput -f {input.snakefile} /gbiomed/home/large_omics_course/output/r1032561/
        touch {output}
        """

rule extract_high_impact:
    input:
        "annotated_snps.vcf"
    output:
        "high_impact.tsv"
    shell:
        "bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT\\t]\\t%INFO/ANN\\n' {input} | \\\n"
        "awk '$5 != $6 && $7 ~ /HIGH/' > {output}"

rule upload_to_db:
    input:
        vcf="annotated_snps.vcf"
    output:
        "snps.sqlite"
    shell:
        """
        bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/ANN\\n' {input.vcf} > snp_rows.tsv &&
        python3 -c "import sqlite3; import pandas as pd; \\
        df = pd.read_csv('snp_rows.tsv', sep='\\t', names=['chrom','pos','ref','alt','ann']); \\
        conn = sqlite3.connect('snps.sqlite'); \\
        conn.execute(\\\"CREATE TABLE IF NOT EXISTS snps (chrom TEXT, pos INT, ref TEXT, alt TEXT, ann TEXT)\\\"); \\
        df.to_sql('snps', conn, if_exists='append', index=False); \\
        conn.commit(); conn.close()"
        """





