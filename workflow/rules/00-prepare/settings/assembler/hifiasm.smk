import re


HIFIASM_PRIMARY_ALT_ASSEMBLY=config.get("hifiasm_primary_alt", False)
assert isinstance(HIFIASM_PRIMARY_ALT_ASSEMBLY, bool)
if HIFIASM_PRIMARY_ALT_ASSEMBLY:
    HIFIASM_PRIMARY_ALT_ASSEMBLY_PARAM="--primary"
else:
    HIFIASM_PRIMARY_ALT_ASSEMBLY_PARAM=""

HIFIASM_PURGE_DUPS_LEVEL=config.get("hifiasm_purge_dups", 1)
assert isinstance(HIFIASM_PURGE_DUPS_LEVEL, int)
assert 0 <= HIFIASM_PURGE_DUPS_LEVEL <= 3

HIFIASM_TELOMERE_MOTIF=config.get("hifiasm_telomere_motif", None)
if HIFIASM_TELOMERE_MOTIF is not None:
    assert isinstance(HIFIASM_TELOMERE_MOTIF, str)
    if HIFIASM_TELOMERE_MOTIF in ["human", "hsa", "homo_sapiens"]:
        HIFIASM_TELOMERE_MOTIF = "CCCTAA"
    # unclear if hifiasm would trip over this ...
    HIFIASM_TELOMERE_MOTIF = HIFIASM_TELOMERE_MOTIF.upper()
    assert re.match("^[ACGT]+$", HIFIASM_TELOMERE_MOTIF) is not None
    HIFIASM_TELOMERE_MOTIF_PARAM="--telo-m " + HIFIASM_TELOMERE_MOTIF
else:
    HIFIASM_TELOMERE_MOTIF_PARAM=""

HIFIASM_DUAL_SCAFFOLDING=config.get("hifiasm_dual_scaffolding", False)
assert isinstance(HIFIASM_DUAL_SCAFFOLDING, bool)
if HIFIASM_DUAL_SCAFFOLDING:
    HIFIASM_DUAL_SCAFFOLDING_PARAM="--dual-scaf"
else:
    HIFIASM_DUAL_SCAFFOLDING_PARAM=""

HIFIASM_DUMP_READ_OVERLAPS=config.get("hifiasm_dump_read_overlaps", False)
assert isinstance(HIFIASM_DUMP_READ_OVERLAPS, bool)
if HIFIASM_DUMP_READ_OVERLAPS:
    HIFIASM_DUMP_READ_OVERLAPS_PARAM="--write-paf"
else:
    HIFIASM_DUMP_READ_OVERLAPS_PARAM=""

HIFIASM_DUMP_EC_READS=config.get("hifiasm_dump_ec_reads", False)
assert isinstance(HIFIASM_DUMP_EC_READS, bool)
if HIFIASM_DUMP_EC_READS:
    HIFIASM_DUMP_EC_READS_PARAM="--write-ec"
else:
    HIFIASM_DUMP_EC_READS_PARAM=""
