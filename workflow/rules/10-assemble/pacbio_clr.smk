
rule flye_assemble_pacbio_clr:
    input:
        reads = lambda wildcards:
            MAP_SAMPLE_TO_INPUT_FILES[(wildcards.sample, "clr")][("fastq", "all")]
    output:
        check = DIR_PROC.joinpath(
            "10-assemble", "flye", "{sample}_clr.ok"
        )
    log:
        DIR_LOG.joinpath(
            "10-assemble", "flye", "{sample}_clr.assm.log"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "10-assemble", "flye", "{sample}_clr.assm.rsrc"
        )
    conda:
        DIR_ENVS.joinpath("assembler", "flye.yaml")
    threads: CPU_MAX
    resources:
        mem_mb=lambda wildcards, attempt: int((768 + 288 * attempt) * 1024),
        time_hrs=lambda wildcards, attempt: 71 * attempt
    params:
        wd=lambda wildcards: pathlib.Path(output.check).with_suffix(".wd"),
        gsize=FLYE_GENOME_SIZE_PARAM,
        keep_hap=FLYE_KEEP_HAPLOTYPES_PARAM,
        no_alt=FLYE_NO_ALT_CONTIGS_PARAM
    shell:
        "flye --pacbio-raw {input.reads} --genome-size {params.gsize} "
            "--threads {threads} --iterations 2 --debug "
            "{params.no_alt} "
            "{params.keep_hap} "
            "--out-dir {params.wd} &> {log}"
                " && "
            "touch {output.check}"


rule run_flye_pacbio_clr_assemblies:
    input:
        checks = expand(
            DIR_PROC.joinpath(
                "assemblies", "flye", "{sample}_clr.ok"
            ),
            sample=CLR_SAMPLES
        )
