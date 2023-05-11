
FLYE_CLR_ASSEMBLIES = []
FLYE_HIFI_ASSEMBLIES = []
FLYE_ONT_ASSEMBLIES = []

if RUN_FLYE:
    if CLR_SAMPLES:
        FLYE_CLR_ASSEMBLIES.extend(
            rules.run_flye_pacbio_clr_assemblies.input.checks
        )
