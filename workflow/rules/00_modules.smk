"""
Use this module to list all includes
required for your pipeline - do not
add your pipeline-specific modules
to "commons/00_commons.smk"
"""

include: "00-prepare/sample_table.smk"
include: "00-prepare/settings.smk"

include: "05-generics/indexing.smk"
include: "05-generics/convert.smk"

include: "10-assemble/pacbio_clr.smk"
include: "10-assemble/pacbio_hifi.smk"

include: "20-postprocess/pacbio_clr_readaln.smk"
include: "20-postprocess/pacbio_clr_hapdup.smk"
include: "20-postprocess/pacbio_clr_polish.smk"

include: "99-outputs/flye.smk"
