"""
Module that implements unphased and phased
assembly for PacBio HiFi reads; currently,
only hifiasm is supported as assembler.
"""

rule hifiasm_assemble_pacbio_hifi_unphased:
    input:
        reads = lambda wildcards:
            MAP_SAMPLE_TO_INPUT_FILES[(wildcards.sample, "hifi")][("fastq", "all")]
    output:
        check = DIR_PROC.joinpath(
            "10-assemble", "hifiasm", f"{HIFIASM_UNPHASED_STEM}.ok"
        ),
        outfiles = HIFIASM_OUTPUT_FILES
    log:
        DIR_LOG.joinpath(
            "10-assemble", "hifiasm", f"{HIFIASM_UNPHASED_STEM}.assm.log"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "10-assemble", "hifiasm", f"{HIFIASM_UNPHASED_STEM}.assm.rsrc"
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


rule hifiasm_assemble_pacbio_hifi_triophased:
    input:
        reads = lambda wildcards:
            MAP_SAMPLE_TO_INPUT_FILES[(wildcards.sample, "hifi")][("fastq", "all")],
        pat_db = lambda wildcards:
            MAP_SAMPLE_TO_INPUT_FILES[(wildcards.sample, "hifi")][("yak", "pat")],
        mat_db = lambda wildcards:
            MAP_SAMPLE_TO_INPUT_FILES[(wildcards.sample, "hifi")][("yak", "mat")],
    output:
        check = DIR_PROC.joinpath(
            "10-assemble", "hifiasm", f"{HIFIASM_TRIO_STEM}.ok"
        ),
        outfiles = HIFIASM_TRIO_OUTPUT_FILES
    log:
        DIR_LOG.joinpath(
            "10-assemble", "hifiasm", f"{HIFIASM_TRIO_STEM}.assm.log"
        )
    benchmark:
        DIR_RSRC.joinpath(
            "10-assemble", "hifiasm", f"{HIFIASM_TRIO_STEM}.assm.rsrc"
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
        "-1 {input.pat_db} -2 {input.mat_db} "
        f"{HIFIASM_DUAL_SCAFFOLDING_PARAM} {HIFIASM_DUMP_EC_READS_PARAM} "
        f"{HIFIASM_DUMP_READ_OVERLAPS_PARAM} "
        f"-l {HIFIASM_TRIO_PURGE_DUPS_LEVEL} {HIFIASM_TELOMERE_MOTIF_PARAM} "
        "{input.reads} &> {log}"
            " && "
        "touch {output.check}"


rule hifiasm_dump_main_assembly_to_fasta:
    """Convert one of the main assembly GFA files
    of any hifiasm run mode to (bgzipped) FASTA.
    """
    input:
        gfa = lambda wildcards: hifiasm_wd_file(
            wildcards.sample, wildcards.hifiasm_mode,
            HIFIASM_MAIN_ASSEMBLIES[wildcards.hifiasm_mode][wildcards.asm_label]
        )
    output:
        fasta = DIR_RES.joinpath(
            "assemblies", "hifiasm", "{sample}.hifi",
            "{sample}.hifi.hifiasm-{hifiasm_mode}.{asm_label}.fasta.gz"
        ),
        fai = DIR_RES.joinpath(
            "assemblies", "hifiasm", "{sample}.hifi",
            "{sample}.hifi.hifiasm-{hifiasm_mode}.{asm_label}.fasta.gz.fai"
        ),
        gzi = DIR_RES.joinpath(
            "assemblies", "hifiasm", "{sample}.hifi",
            "{sample}.hifi.hifiasm-{hifiasm_mode}.{asm_label}.fasta.gz.gzi"
        ),
    wildcard_constraints:
        sample=CONSTRAINT_HIFI_SAMPLES,
        hifiasm_mode=CONSTRAINT_HIFIASM_RUN_MODES,
        asm_label=CONSTRAINT_HIFIASM_ASM_LABELS,
    conda:
        DIR_ENVS.joinpath("biotools.yaml")
    threads: CPU_LOW
    resources:
        mem_mb=lambda wildcards, attempt: 8192 * attempt
    shell:
        "gfatools gfa2fa -l 0 {input.gfa} | bgzip -c --threads {threads} > {output.fasta}"
            " && "
        "samtools faidx {output.fasta}"


# Samples are either trio-phased (paternal and maternal yak
# k-mer databases given in the sample sheet) or assembled with
# the configured unphased run mode (unps or rdps), not both.
HIFI_TRIO_SAMPLES = sorted(
    sample for (sample, read_type) in TRIO_PHASED_SAMPLES
    if read_type == "hifi"
)
HIFI_UNPHASED_SAMPLES = sorted(
    sample for (sample, read_type) in UNPHASED_SAMPLES
    if read_type == "hifi"
)


def hifiasm_sample_run_mode(sample):
    if sample in HIFI_TRIO_SAMPLES:
        return HIFIASM_TRIO_RUN_MODE
    return HIFIASM_UNPHASED_RUN_MODE


# The error-corrected reads and the read overlaps are identical
# for all run modes (hifiasm can reuse the *.bin files across
# run modes), hence they are only saved once per sample,
# taken from whichever run mode is used for the sample.

HIFIASM_ASSEMBLY_RESULT_FILES = []


if HIFIASM_DUMP_EC_READS:

    rule compress_hifiasm_ec_reads:
        input:
            fasta = lambda wildcards: hifiasm_wd_file(
                wildcards.sample, hifiasm_sample_run_mode(wildcards.sample), ".ec.fa"
            )
        output:
            fagz = DIR_RES.joinpath(
                "assemblies", "hifiasm",
                "{sample}.hifi", "aux", "{sample}.hifi.hifiasm-ec-reads.fa.gz"
            )
        wildcard_constraints:
            sample=CONSTRAINT_HIFI_SAMPLES,
        conda:
            DIR_ENVS.joinpath("biotools.yaml")
        threads: CPU_LOW
        resources:
            time_hrs=lambda wildcards, attempt: attempt * attempt
        shell:
            "pigz -p {threads} -c {input.fasta} > {output.fagz}"

    HIFIASM_ASSEMBLY_RESULT_FILES.append(
        rules.compress_hifiasm_ec_reads.output.fagz
    )


if HIFIASM_DUMP_READ_OVERLAPS:

    rule compress_hifiasm_read_overlaps:
        input:
            paf = lambda wildcards: hifiasm_wd_file(
                wildcards.sample, hifiasm_sample_run_mode(wildcards.sample), ".ovlp.paf"
            )
        output:
            paf = DIR_RES.joinpath(
                "assemblies", "hifiasm",
                "{sample}.hifi", "aux", "{sample}.hifi.hifiasm-read-overlaps.paf.gz"
            )
        wildcard_constraints:
            sample=CONSTRAINT_HIFI_SAMPLES,
        conda:
            DIR_ENVS.joinpath("biotools.yaml")
        threads: CPU_LOW
        resources:
            time_hrs=lambda wildcards, attempt: attempt * attempt
        shell:
            "pigz -p {threads} -c {input.paf} > {output.paf}"

    HIFIASM_ASSEMBLY_RESULT_FILES.append(
        rules.compress_hifiasm_read_overlaps.output.paf
    )


rule run_hifiasm_pacbio_hifi_assemblies:
    input:
        assemblies = expand(
            rules.hifiasm_assemble_pacbio_hifi_unphased.output.check,
            sample=HIFI_UNPHASED_SAMPLES,
        ),
        trio_assemblies = expand(
            rules.hifiasm_assemble_pacbio_hifi_triophased.output.check,
            sample=HIFI_TRIO_SAMPLES,
        ),
        main_fasta = expand(
            rules.hifiasm_dump_main_assembly_to_fasta.output.fasta,
            sample=HIFI_UNPHASED_SAMPLES,
            hifiasm_mode=[HIFIASM_UNPHASED_RUN_MODE],
            asm_label=HIFIASM_MAIN_WILDCARDS,
        ),
        trio_main_fasta = expand(
            rules.hifiasm_dump_main_assembly_to_fasta.output.fasta,
            sample=HIFI_TRIO_SAMPLES,
            hifiasm_mode=[HIFIASM_TRIO_RUN_MODE],
            asm_label=HIFIASM_TRIO_MAIN_WILDCARDS,
        ),
        hifiasm_files = expand(
            HIFIASM_ASSEMBLY_RESULT_FILES,
            sample=HIFI_SAMPLES
        )
