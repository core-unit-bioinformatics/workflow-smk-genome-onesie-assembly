# TODO
# this module needs to be updated/extende
# to include the different phasing options for hifi


HIFIASM_UNPHASED_OUTPUT_EXTENSIONS = [
    ".bp.r_utg.gfa", ".bp.p_utg.gfa",
    ".bp.r_utg.noseq.gfa", ".bp.p_utg.noseq.gfa",
    ".ovlp.reverse.bin", ".ovlp.source.bin", ".ec.bin"
]

if HIFIASM_DUMP_EC_READS:
    HIFIASM_UNPHASED_OUTPUT_EXTENSIONS.append(".ec.fa")

if HIFIASM_DUMP_READ_OVERLAPS:
    HIFIASM_UNPHASED_OUTPUT_EXTENSIONS.append(".ovlp.paf")

if HIFIASM_PRIMARY_ALT_ASSEMBLY:
    HIFIASM_UNPHASED_OUTPUT_EXTENSIONS.append(".bp.p_ctg.gfa")
    HIFIASM_UNPHASED_OUTPUT_EXTENSIONS.append(".bp.a_ctg.gfa")
    HIFIASM_UNPHASED_MAIN_WILDCARDS = ["prm", "alt"]
    HIFIASM_UNPHASED_MAIN_EXT = [".bp.p_ctg.gfa", ".bp.a_ctg.gfa"]
else:
    HIFIASM_UNPHASED_OUTPUT_EXTENSIONS.append(".bp.hap1.p_ctg.gfa")
    HIFIASM_UNPHASED_OUTPUT_EXTENSIONS.append(".bp.hap2.p_ctg.gfa")
    HIFIASM_UNPHASED_MAIN_WILDCARDS = ["hap1", "hap2"]
    HIFIASM_UNPHASED_MAIN_EXT = [".bp.hap1.p_ctg.gfa", ".bp.hap2.p_ctg.gfa"]

HIFIASM_OUTPUT_FILES = [
    DIR_PROC.joinpath(
        "10-assemble", "hifiasm", "{sample}_hifi.wd",
        f"{{sample}}{ext}"
    ) for ext in HIFIASM_UNPHASED_OUTPUT_EXTENSIONS
]


rule hifiasm_assemble_pacbio_hifi_unphased:
    input:
        reads = lambda wildcards:
            MAP_SAMPLE_TO_INPUT_FILES[(wildcards.sample, "hifi")][("fastq", "all")]
    output:
        check = DIR_PROC.joinpath(
            "10-assemble", "hifiasm", "{sample}_hifi.ok"
        ),
        outfiles = HIFIASM_OUTPUT_FILES
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
        time_hrs=lambda wildcards, attempt: 23 * attempt
    params:
        prefix=lambda wildcards, output: pathlib.Path(output.check).with_suffix(".wd").joinpath(wildcards.sample),
        wd=lambda wildcards, output: pathlib.Path(output.check).with_suffix(".wd"),
    shell:
        "mkdir -p {params.wd}"
            " && "
        "hifiasm -t {threads} -o {params.prefix} "
        f"{HIFIASM_DUAL_SCAFFOLDING_PARAM} {HIFIASM_DUMP_EC_READS_PARAM} "
        f"{HIFIASM_DUMP_READ_OVERLAPS_PARAM} {HIFIASM_PRIMARY_ALT_ASSEMBLY_PARAM} "
        f"-l {HIFIASM_PURGE_DUPS_LEVEL} {HIFIASM_TELOMERE_MOTIF_PARAM} "
        "{input.reads} &> {log}"
            " && "
        "touch {output.check}"


rule hifiasm_dump_main_assembly_to_fasta:
    input:
        gfa1 = DIR_PROC.joinpath(
            "10-assemble", "hifiasm", "{sample}_hifi.wd",
            "{sample}" + f"{HIFIASM_UNPHASED_MAIN_EXT[0]}"
        ),
        gfa2 = DIR_PROC.joinpath(
            "10-assemble", "hifiasm", "{sample}_hifi.wd",
            "{sample}" + f"{HIFIASM_UNPHASED_MAIN_EXT[1]}"
        ),
    output:
        fasta1 = DIR_RES.joinpath(
            "assemblies", "hifiasm", "{sample}_hifi",
            "{sample}_hifi.hifiasm.asm-" + f"{[HIFIASM_UNPHASED_MAIN_WILDCARDS[0]]}.fasta.gz"
        ),
        fasta2 = DIR_RES.joinpath(
            "assemblies", "hifiasm", "{sample}_hifi",
            "{sample}_hifi.hifiasm.asm-" + f"{[HIFIASM_UNPHASED_MAIN_WILDCARDS[1]]}.fasta.gz"
        )
    conda:
        DIR_ENVS.joinpath("biotools.yaml")
    threads: CPU_LOW
    resources:
        mem_mb=lambda wildcards, attempt: 2048 * attempt
    shell:
        "gfatools -l 0 {input.gfa1} | bgzip -c --threads {threads} > {output.fasta1}"
            " && "
        "samtools faidx {output.fasta1}"
            " && "
        "gfatools -l 0 {input.gfa2} | bgzip -c --threads {threads} > {output.fasta2}"
            " && "
        "samtools faidx {output.fasta2}"



HIFIASM_ASSEMBLY_RESULT_FILES = []


if HIFIASM_DUMP_EC_READS:

    rule compress_hifiasm_ec_reads:
        input:
            fasta = DIR_PROC.joinpath(
                "10-assemble", "hifiasm", "{sample}_hifi.wd",
                "{sample}.ec.fa"
            )
        output:
            fagz = DIR_RES.joinpath(
                "assemblies", "hifiasm",
                "{sample}_hifi", "aux", "{sample}_hifi.hifiasm-ec-reads.fa.gz"
            )
        conda:
            DIR_ENVS.joinpath("biotools.yaml")
        threads: CPU_LOW
        shell:
            "pigz -p {threads} -c {input.fasta} > {output.fagz}"

    HIFIASM_ASSEMBLY_RESULT_FILES.append(
        rules.compress_hifiasm_ec_reads.output.fagz
    )


if HIFIASM_DUMP_READ_OVERLAPS:

    rule compress_hifiasm_read_overlaps:
        input:
            paf = DIR_PROC.joinpath(
                "10-assemble", "hifiasm", "{sample}_hifi.wd",
                "{sample}.ovlp.paf"
            )
        output:
            paf = DIR_RES.joinpath(
                "assemblies", "hifiasm",
                "{sample}_hifi", "aux", "{sample}_hifi.hifiasm-read-overlaps.paf.gz"
            )
        conda:
            DIR_ENVS.joinpath("biotools.yaml")
        threads: CPU_LOW
        shell:
            "pigz -p {threads} -c {input.paf} > {output.paf}"

    HIFIASM_ASSEMBLY_RESULT_FILES.append(
        rules.compress_hifiasm_read_overlaps.paf
    )




rule run_hifiasm_pacbio_hifi_assemblies:
    input:
        assemblies = expand(
            rules.hifiasm_assemble_pacbio_hifi.output.check,
            sample=HIFI_SAMPLES,
        ),
        main_fasta = rules.hifiasm_dump_main_assembly_to_fasta.output,
        hifiasm_files = HIFIASM_ASSEMBLY_RESULT_FILES
