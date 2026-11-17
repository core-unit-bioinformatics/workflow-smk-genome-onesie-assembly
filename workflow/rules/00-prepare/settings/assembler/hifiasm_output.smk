
"""
This module exists only to capture the "complexity"
of hifiasm output file name patterns.

The different variants are specified here (as of v0.25.0)
https://github.com/chhylp123/hifiasm/blob/ec9a8b222d149d25b7355e83765698640d59b189/Overlaps.cpp#L39140

Note that the infix is only added to the assembly graph / contig
output files (*.r_utg.*, *.p_utg.*, *.p_ctg.* etc.). The binary
files (*.bin) and the optional *.ec.fa / *.ovlp.paf dumps always
use the plain output prefix (no infix).

This workflow supports

1) the default run mode: creates two haplotype assemblies phased
using only the HiFi reads (no external phasing information),
i.e. a "bipartition" (guess for what 'bp' in the output stands for)
creates two somewhat balanced haplotypes. Expect (a lot of) phasing
errors / switches between the parental haplotypes. In addition to
hap1/hap2, a primary assembly (.bp.p_ctg.gfa) is created.
Influenced by the '-l' / purging parameter; note that '-l 0' (no purging)
is intended for assembling inbred/homozygous genomes (see hifiasm README),
for which a haplotype bipartition is meaningless. Hence, in the default
run mode, '-l 0' makes hifiasm fall back to primary/alt output (with unverified file
names, probably keeping the '.bp.' infix). The expected hap1/hap2 output
files would then be missing, so this combination is rejected at workflow
startup (see ValueError below).
- output files carry a '.bp.' infix
HIFIASM_BIPART_MODE

2) the primary/alt run mode: only activated via the '--primary' switch,
and creates a collapsed, 'reference-style' assembly plus the leftover
tigs in the alt set. Influenced by the '-l' parameter: the purge level
determines which haplotigs are moved from the primary to the alt set
(hifiasm default: 3).
- output files carry no infix
HIFIASM_PRMALT_MODE

3) trio-based phasing: create properly phased diploid assemblies with
hap1/paternal and hap2/maternal assemblies. Influenced by the '-l' parameter,
but only levels 0 and 1 are valid (hifiasm default: 0).
- output files carry a '.dip.' infix
HIFIASM_TRIO_MODE

4) (to be implemented): HiC-based phasing
- output files carry a '.hic.' infix
HIFIASM_HIC_MODE

"""

# modes 1) and 2) are mutually exclusive and are selected
# via the 'hifiasm_primary_alt' config option
HIFIASM_PRMALT_MODE = HIFIASM_PRIMARY_ALT_ASSEMBLY
HIFIASM_BIPART_MODE = not HIFIASM_PRIMARY_ALT_ASSEMBLY

if HIFIASM_BIPART_MODE and HIFIASM_PURGE_DUPS_LEVEL == 0:
    raise ValueError(
        "hifiasm default (bipartition) run mode with purge level 0 "
        "falls back to primary/alt output, which is not supported. "
        "Set 'hifiasm_purge_dups' > 0 or 'hifiasm_primary_alt: true'. "
        "For inbred/homozygous samples, use 'hifiasm_primary_alt: true'."
    )

# run mode labels used in output paths:
# {sample}.hifi.{hifiasm_mode}[.wd|.ok|...]
# unps = unphased (primary/alt)
# rdps = read-phased (bipartition),
# trio = trio-phased
if HIFIASM_PRMALT_MODE:
    HIFIASM_UNPHASED_RUN_MODE = "unps"
else:
    HIFIASM_UNPHASED_RUN_MODE = "rdps"
HIFIASM_TRIO_RUN_MODE = "trio"

# file name stems (no extension) for the assembly
# working directories, check files, logs etc.
HIFIASM_UNPHASED_STEM = f"{{sample}}.hifi.{HIFIASM_UNPHASED_RUN_MODE}"
HIFIASM_TRIO_STEM = f"{{sample}}.hifi.{HIFIASM_TRIO_RUN_MODE}"

# output files independent of run mode / no infix
HIFIASM_COMMON_OUTPUT_EXTENSIONS = [
    ".ovlp.reverse.bin", ".ovlp.source.bin", ".ec.bin"
]

if HIFIASM_DUMP_EC_READS:
    HIFIASM_COMMON_OUTPUT_EXTENSIONS.append(".ec.fa")

if HIFIASM_DUMP_READ_OVERLAPS:
    HIFIASM_COMMON_OUTPUT_EXTENSIONS.append(".ovlp.paf")


def _hifiasm_graph_extensions(infix):
    return [
        f"{infix}.r_utg.gfa", f"{infix}.p_utg.gfa",
        f"{infix}.r_utg.noseq.gfa", f"{infix}.p_utg.noseq.gfa",
    ]


# 1)
# entirely unphased / primary-alt mode
if HIFIASM_PRMALT_MODE:
    HIFIASM_MAIN_WILDCARDS = ["prm", "alt"]
    HIFIASM_MAIN_EXT = [".p_ctg.gfa", ".a_ctg.gfa"]
    HIFIASM_OUTPUT_EXTENSIONS = _hifiasm_graph_extensions("")
else:
# 2) read-based phasing / bipartition
    assert HIFIASM_BIPART_MODE
    HIFIASM_MAIN_WILDCARDS = ["bip1", "bip2"]
    HIFIASM_MAIN_EXT = [".bp.hap1.p_ctg.gfa", ".bp.hap2.p_ctg.gfa"]
    HIFIASM_OUTPUT_EXTENSIONS = _hifiasm_graph_extensions(".bp")
    HIFIASM_OUTPUT_EXTENSIONS.append(".bp.p_ctg.gfa")

HIFIASM_OUTPUT_EXTENSIONS.extend(HIFIASM_MAIN_EXT)
HIFIASM_OUTPUT_EXTENSIONS.extend(HIFIASM_COMMON_OUTPUT_EXTENSIONS)

HIFIASM_OUTPUT_FILES = [
    DIR_PROC.joinpath(
        "10-assemble", "hifiasm", f"{HIFIASM_UNPHASED_STEM}.wd",
        f"{{sample}}{ext}"
    ) for ext in HIFIASM_OUTPUT_EXTENSIONS
]


# 3) trio-phased assemblies
# in trio mode, hap1 is paternal and hap2 is maternal
HIFIASM_TRIO_MAIN_WILDCARDS = ["pat", "mat"]
HIFIASM_TRIO_MAIN_EXT = [".dip.hap1.p_ctg.gfa", ".dip.hap2.p_ctg.gfa"]

HIFIASM_TRIO_OUTPUT_EXTENSIONS = _hifiasm_graph_extensions(".dip")
HIFIASM_TRIO_OUTPUT_EXTENSIONS.extend(HIFIASM_TRIO_MAIN_EXT)
HIFIASM_TRIO_OUTPUT_EXTENSIONS.extend(HIFIASM_COMMON_OUTPUT_EXTENSIONS)

HIFIASM_TRIO_OUTPUT_FILES = [
    DIR_PROC.joinpath(
        "10-assemble", "hifiasm", f"{HIFIASM_TRIO_STEM}.wd",
        f"{{sample}}{ext}"
    ) for ext in HIFIASM_TRIO_OUTPUT_EXTENSIONS
]


# lookup tables for downstream rules (wildcard 'hifiasm_mode')
# only the active unphased run mode (unps or rdps) is included
HIFIASM_RUN_MODE_STEMS = {
    HIFIASM_UNPHASED_RUN_MODE: HIFIASM_UNPHASED_STEM,
    HIFIASM_TRIO_RUN_MODE: HIFIASM_TRIO_STEM,
}

# main assembly output per run mode: label -> file extension
HIFIASM_MAIN_ASSEMBLIES = {
    HIFIASM_UNPHASED_RUN_MODE: dict(zip(HIFIASM_MAIN_WILDCARDS, HIFIASM_MAIN_EXT)),
    HIFIASM_TRIO_RUN_MODE: dict(zip(HIFIASM_TRIO_MAIN_WILDCARDS, HIFIASM_TRIO_MAIN_EXT)),
}

CONSTRAINT_HIFIASM_RUN_MODES = "(" + "|".join(sorted(HIFIASM_RUN_MODE_STEMS.keys())) + ")"
CONSTRAINT_HIFIASM_ASM_LABELS = "(" + "|".join(
    sorted(set(HIFIASM_MAIN_WILDCARDS + HIFIASM_TRIO_MAIN_WILDCARDS))
) + ")"


