import pathlib

DATA_ROOT = pathlib.Path(config.get("data_root", "/"))

RUN_FLYE = config.get("run_flye", False)
assert isinstance(RUN_FLYE, bool)
RUN_FLYE_HAPDUP = config.get("run_flye_hapdup", False)
assert isinstance(RUN_FLYE_HAPDUP, bool)

if RUN_FLYE:
    FLYE_GENOME_SIZE_PARAM = config["flye_genome_size"]
    assert isinstance(FLYE_GENOME_SIZE_PARAM, str)
else:
    FLYE_GENOME_SIZE_PARAM = ""

FLYE_KEEP_HAPLOTYPES = config.get("flye_keep_haplotypes", False)
assert isinstance(FLYE_KEEP_HAPLOTYPES, bool)
if FLYE_KEEP_HAPLOTYPES:
    FLYE_KEEP_HAPLOTYPES_PARAM = " --keep-haplotypes "
else:
    FLYE_KEEP_HAPLOTYPES_PARAM = ""

FLYE_NO_ALT_CONTIGS = config.get("flye_no_alt_contigs", False)
assert isinstance(FLYE_NO_ALT_CONTIGS, bool)
if FLYE_NO_ALT_CONTIGS:
    FLYE_NO_ALT_CONTIGS_PARAM = " --no-alt-contigs "
else:
    FLYE_NO_ALT_CONTIGS_PARAM = ""
