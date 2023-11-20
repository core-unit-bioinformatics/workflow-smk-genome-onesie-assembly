

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
    threads: lambda wildcards, input: hifiasm_cpu_cores(input.size_mb)
    resources:
        mem_mb=lambda wildcards, attempt, input: hifiasm_memory_mb(input.size_mb) * attempt,
        time_hrs=lambda wildcards, attempt: 71 * attempt
    params:
        prefix=lambda wildcards, output: pathlib.Path(output.check).with_suffix(".wd").joinpath(wildcards.sample),
        wd=lambda wildcards, output: pathlib.Path(output.check).with_suffix(".wd"),
    shell:
        "mkdir -p {params.wd}"
            " && "
        "hifiasm -t {threads} -o {params.prefix} {input.reads} &> {log}"
            " && "
        "touch {output.check}"


rule run_hifiasm_pacbio_hifi_assemblies:
    input:
        assemblies = expand(
            rules.hifiasm_assemble_pacbio_hifi.output.check,
            sample=HIFI_SAMPLES,
        )
