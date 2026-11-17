
def hifiasm_memory_mb(input_size_mb):

    input_gb = input_size_mb / 1024

    if input_gb < 100:
        mem_gb = 160
    elif input_gb < 200:
        mem_gb = 320
    else:
        mem_gb = 960
    mem_mb = int(mem_gb * 1024)
    return mem_mb


def hifiasm_cpu_cores(input_size_mb):

    input_gb = input_size_mb / 1024

    if input_gb < 100:
        cores = min(24, CPU_MAX)
    elif input_gb < 200:
        cores = min(36, CPU_MAX)
    else:
        cores = min(48, CPU_MAX)
    return cores


def hifiasm_wd_file(sample, run_mode, ext):
    """Path to a file in the hifiasm working directory
    of the given sample and run mode, e.g. ext='.ec.fa'
    """
    stem = HIFIASM_RUN_MODE_STEMS[run_mode].format(sample=sample)
    return DIR_PROC.joinpath(
        "10-assemble", "hifiasm", f"{stem}.wd", f"{sample}{ext}"
    )
