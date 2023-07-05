

localrules: pacbio_create_bam_fofn_file
rule pacbio_create_bam_fofn_file:
    input:
        reads = lambda wildcards:
            MAP_SAMPLE_TO_INPUT_FILES[(wildcards.sample, "clr")][("bam", "all")]
    output:
        fofn = DIR_PROC.joinpath(
            "20-postprocess", "bam_fofn", "{sample}_bam_clr.fofn"
        )
    run:
        import pathlib as pl

        assert len(input.reads) > 0

        file_paths = []
        for read_file in input.reads:
            assert pl.Path(read_file).is_file()
            file_paths.append(str(read_file))

        with open(output.fofn, "w") as fofn_dump:
            _ = fofn_dump.write("\n".join(sorted(file_paths)) + "\n")
    # END OF RUN BLOCK


rule pbmm2_produce_read_assembly_alignments:
    """
    This alignment step can either be used to
    run one step of assembly polishing, or to
    first split a flye CLR into two haplotypes
    with the hapdup tool.
    """
    input:
        asm=DIR_PROC.joinpath(
            "10-assemble", "flye", "{sample}_clr.wd", "assembly.fasta"
        ),
        fai=DIR_PROC.joinpath(
            "10-assemble", "flye", "{sample}_clr.wd", "assembly.fasta.fai"
        ),
        reads_fofn = DIR_PROC.joinpath(
            "20-postprocess", "bam_fofn", "{sample}_bam_clr.fofn"
        )
    output:
        bam=DIR_PROC.joinpath(
            "20-postprocess", "pacbio_clr", "read_asm_align",
            "{sample}_clr.sort.bam"
        ),
        bai=DIR_PROC.joinpath(
            "20-postprocess", "pacbio_clr", "read_asm_align",
            "{sample}_clr.sort.bam.bai"
        )
    log:
        DIR_LOG.joinpath(
            "20-postprocess", "pacbio_clr", "read_asm_align",
            "{sample}_clr.pbmm2.log"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "20-postprocess", "pacbio_clr", "read_asm_align",
            "{sample}_clr.pbmm2.rsrc"
        )
    conda:
        DIR_ENVS.joinpath("pbtools.yaml")
    threads: CPU_HIGH
    resources:
        mem_mb=lambda wc, attempt: int((32 + 32 * attempt) * 1024),
        time_hrs=lambda wc, attempt: 11 * attempt,
    params:
        sort_mem=4096,
        sort_threads=CPU_LOW,
        tempdir=lambda wildcards: DIR_PROC.joinpath(
            "tmp", "20-postprocess", "pacbio_clr", "read_asm_align",
            f"{wildcards.sample}_clr.sort.wd"
        )
    shell:
        "mkdir -p {params.tempdir} && TMPDIR={params.tempdir} "
        "pbmm2 align --sort --sort-memory {params.sort_mem}M --sort-threads {params.sort_threads} "
            "--num-threads {threads} --log-level DEBUG --log-file {log} "
            "--preset SUBREAD --sample {wildcards.sample} --bam-index BAI "
            "{input.asm} {input.reads_fofn} {output.bam}"
        " ; rm -rfd {params.tempdir}"
