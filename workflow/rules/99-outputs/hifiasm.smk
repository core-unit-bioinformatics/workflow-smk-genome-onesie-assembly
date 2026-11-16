HIFIASM_HIFI_ASSEMBLIES = []
HIFIASM_ONT_ASSEMBLIES = []

if RUN_HIFIASM:
    if HIFI_SAMPLES:
        HIFIASM_HIFI_ASSEMBLIES = rules.run_hifiasm_pacbio_hifi_assemblies.input
