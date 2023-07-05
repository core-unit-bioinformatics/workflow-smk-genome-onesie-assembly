
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
        bai=DIR_PROC.joinpath(
            "20-postprocess", "pacbio_clr", "read_asm_align",
            "{sample}_clr.sort.bam.bai"
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
