
def hifiasm_memory_gb(input_size_mb):

    input_gb = input_size_mb / 1024

    if input_gb < 100:
        run_mem = 160
    elif input_gb < 200:
        run_mem = 320
    else:
        run_mem = 480
    return run_mem


def hifiasm_cpu_cores(input_size_mb):

    input_gb = input_size_mb / 1024

    if input_gb < 100:
        cores = min(24, CPU_MAX)
    elif input_gb < 200:
        cores = min(36, CPU_MAX)
    else:
        cores = min(48, CPU_MAX)
    return cores
