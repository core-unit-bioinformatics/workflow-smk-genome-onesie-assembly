
rule samtools_index_fasta:
    input:
        "{filepath}.fasta"
    output:
        "{filepath}.fasta.fai"
    conda:
        DIR_ENVS.joinpath("biotools.yaml")
    resources:
        mem_mb=lambda wc, attempt: 2048 * attempt,
        time_hrs=lambda wc, attempt: 1 * attempt
    shell:
        "samtools faidx {input}"


rule index_pacbio_bam_file:
    input:
        bam = "{filepath}.bam"
    output:
        pbi = "{filepath}.bam.pbi"
    conda:
        DIR_ENVS.joinpath("pbtools.yaml")
    threads: CPU_LOW
    resources:
        mem_mb=lambda wc, attempt: 2048 * attempt,
        time_hrs=lambda wc, attempt: 1 * attempt
    shell:
        "pbindex --num-threads {threads} {input.bam}"
