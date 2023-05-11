
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


rule dump_bam_to_fastq:
    input:
        bam = lambda wc:
            MAP_SAMPLE_TO_INPUT_FILES[(wc.sample, wc.read_type)][("bam", wc.path_id)],
        pbi = lambda wc:
            pathlib.Path(
                MAP_SAMPLE_TO_INPUT_FILES[(wc.sample, wc.read_type)][("bam", wc.path_id)]
            ).with_suffix(".bam.pbi")
    output:
        fastq = DIR_PROC.joinpath(
            "00-prepare", "dump_fastq", "{sample}_{read_type}.{path_id}.fastq.gz"
        )
    log:
        DIR_LOG.joinpath(
            "00-prepare", "dump_fastq", "{sample}_{read_type}.{path_id}.dump.log"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "00-prepare", "dump_fastq", "{sample}_{read_type}.{path_id}.dump.rsrc"
        )
    conda:
        DIR_ENVS.joinpath("pbtools.yaml")
    threads: CPU_LOW
    resources:
        mem_mb=lambda wc, attempt: 4096 * attempt,
        time_hrs=lambda wc, attempt: 23 * attempt
    params:
        prefix=lambda wc, output: str(output.fastq).rsplit(".", 2)[0]
    shell:
        "bam2fastq --output {params.prefix} --num-threads {threads} "
            "-c 9 {input.bam} &> {log}"
