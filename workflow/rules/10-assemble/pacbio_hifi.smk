

rule hifiasm_assemble_pacbio_hifi:
    input:
        reads = lambda wildcards:
            MAP_SAMPLE_TO_INPUT_FILES[(wildcards.sample, "hifi")][("fastq", "all")]
    output:
        check = DIR_PROC.joinpath(
            "10-assemble", "hifiasm", "{sample}_hifi.ok"
        ),
    log:
        DIR_LOG.joinpath(
            "10-assemble", "hifiasm", "{sample}_hifi.assm.log"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "10-assemble", "hifiasm", "{sample}_hifi.assm.rsrc"
        )
    conda:
        DIR_ENVS.joinpath("assembler", "hifiasm.yaml")
    threads: CPU_HIGH
    resources:
        mem_mb=lambda wildcards, attempt: int((160 * attempt) * 1024),
        time_hrs=lambda wildcards, attempt: 71 * attempt
    params:
        prefix=lambda wildcards, output: pathlib.Path(output.check).with_suffix(".wd").joinpath(wildcards.sample),
    shell:
        "hifiasm -t {threads} -o {params.prefix} {input.reads} &> {log}"
            " && "
        "touch {output.check}"


rule run_hifiasm_pacbio_hifi_assemblies:
    input:
        assemblies = expand(
            rules.hifiasm_assemble_pacbio_hifi.output.check,
            sample=HIFI_SAMPLES,
        )
