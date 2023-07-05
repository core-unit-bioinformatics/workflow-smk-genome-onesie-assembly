
rule pacbio_clr_flye_hapdup_assembly:
    """
    NB: hapdup github readme claims compatibility
    with both HiFi and CLR for PacBio, but only a
    Hifi preset exists, and the underlying variant
    calling models for PEPPER/MARGIN are also only
    hard-coded for ONT or HiFi.
    """
    input:
        asm = DIR_PROC.joinpath(
            "10-assemble", "flye", "{sample}_clr.wd", "assembly.fasta",
        ),
        idx = DIR_PROC.joinpath(
            "10-assemble", "flye", "{sample}_clr.wd", "assembly.fasta.fai",
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
        chk = DIR_PROC.joinpath(
            "20-postprocess", "pacbio_clr", "flye_hapdup",
            "{sample}_clr.hapdup.ok"
        )
    log:
        DIR_LOG.joinpath(
            "20-postprocess", "pacbio_clr", "flye_hapdup",
            "{sample}_clr.hapdup.log"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "20-postprocess", "pacbio_clr", "flye_hapdup",
            "{sample}_clr.hapdup.rsrc"
        )
    singularity:
        "/gpfs/project/projects/medbioinf/container/hapdup_0.12.sif"
        #TODO move to env config
    threads: 24
    resources:
        time_hrs = lambda wildcards, input, attempt: 47 * attempt,
        mem_mb = lambda wildcards, input, attempt: 32768 + 32768 * attempt,
    params:
        preset = 'ont',  # see rule docstring / closest to CLR?
        outdir = lambda wildcards, output: pathlib.Path(output.chk).with_suffix(".wd"),
    shell:
        'hapdup --assembly {input.asm} --bam {input.bam} '
            '--out-dir {params.outdir} -t {threads} --rtype {params.preset}'
            ' &> {log} && touch {output.chk}'


rule run_pacbio_clr_flye_hapdup_assembly:
    input:
        hapdup_chk = expand(
            DIR_PROC.joinpath(
                "20-postprocess", "pacbio_clr", "flye_hapdup",
                "{sample}_clr.hapdup.ok"
            ),
            sample=CLR_SAMPLES,
        )
