

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


rule pbmm2_produce_polishing_alignments:
    """
    This type of polishing can only be done for CLR
    assemblies. Hence, everything can be fixed to
    that scenario.
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


rule gcpp_assembly_polishing_pass1:
    input:
        asm=DIR_PROC.joinpath(
            "10-assemble", "flye", "{sample}_clr.wd", "assembly.fasta"
        ),
        fai=DIR_PROC.joinpath(
            "10-assemble", "flye", "{sample}_clr.wd", "assembly.fasta.fai"
        ),
        bam=DIR_PROC.joinpath(
            "20-postprocess", "pacbio_clr", "read_asm_align",
            "{sample}_clr.sort.bam"
        ),
        pbi=DIR_PROC.joinpath(
            "20-postprocess", "pacbio_clr", "read_asm_align",
            "{sample}_clr.sort.bam.pbi"
        )
    output:
        asm=DIR_PROC.joinpath(
            "20-postprocess", "pacbio_clr", "asm_polishing", "{sample}_clr.pass1.wd",
            "{sample}_clr.pass1.fasta"
        ),
        gff=DIR_PROC.joinpath(
            "20-postprocess", "pacbio_clr", "asm_polishing", "{sample}_clr.pass1.wd",
            "{sample}_clr.pass1.gff"
        ),
        vcf=DIR_PROC.joinpath(
            "20-postprocess", "pacbio_clr", "asm_polishing", "{sample}_clr.pass1.wd",
            "{sample}_clr.pass1.vcf"
        )
    log:
        DIR_LOG.joinpath(
            "20-postprocess", "pacbio_clr", "asm_polishing", "{sample}_clr.pass1.gcpp.log"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "20-postprocess", "pacbio_clr", "asm_polishing", "{sample}_clr.pass1.gcpp.rsrc"
        )
    conda:
        DIR_ENVS.joinpath("pbtools.yaml")
    threads: CPU_HIGH
    resources:
        mem_mb=lambda wc, attempt: int((64 + 32 * attempt) * 1024),
        time_hrs=lambda wc, attempt: 71 * attempt,
    shell:
        "gcpp --num-threads {threads} --reference {input.asm} --algorithm arrow "
            "--log-level DEBUG --log-file {log} "
            "--output {output.asm},{output.gff},{output.vcf} {input.bam}"
