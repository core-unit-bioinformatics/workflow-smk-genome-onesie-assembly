import pathlib
import hashlib
import collections
import re

import pandas

SAMPLES = None
MAP_SAMPLE_TO_INPUT_FILES = None

CLR_SAMPLES = None
CONSTRAINT_CLR_SAMPLES = None
HIFI_SAMPLES = None
CONSTRAINT_HIFI_SAMPLES = None
ONT_SAMLES = None
CONSTRAINT_ONT_SAMPLES = None


def process_sample_sheet():

    SAMPLE_SHEET_FILE = pathlib.Path(config["samples"]).resolve(strict=True)
    SAMPLE_SHEET = pandas.read_csv(
        SAMPLE_SHEET_FILE,
        sep="\t",
        header=0
    )

    # step 1: each row is a sample,
    # just collect the input files
    sample_input, sample_lists = collect_input_files(SAMPLE_SHEET)
    num_samples = collections.Counter(
        sample for (sample, _) in sample_input.keys()
    ).total()
    assert num_samples == SAMPLE_SHEET.shape[0]
    all_samples = sorted(set(sample for (sample, _) in sample_input.keys()))

    global SAMPLES
    SAMPLES = all_samples

    global CLR_SAMPLES
    CLR_SAMPLES = sorted(sample_lists["clr"])
    global HIFI_SAMPLES
    HIFI_SAMPLES = sorted(sample_lists["hifi"])
    global ONT_SAMPLES
    ONT_SAMPLES = sorted(sample_lists["ont"])

    global MAP_SAMPLE_TO_INPUT_FILES
    MAP_SAMPLE_TO_INPUT_FILES = sample_input

    global CONSTRAINT_CLR_SAMPLES
    CONSTRAINT_CLR_SAMPLES = _build_constraint(CLR_SAMPLES)
    global CONSTRAINT_HIFI_SAMPLES
    CONSTRAINT_HIFI_SAMPLES = _build_constraint(HIFI_SAMPLES)
    global CONSTRAINT_ONT_SAMPLES
    CONSTRAINT_ONT_SAMPLES = _build_constraint(ONT_SAMPLES)

    return


def collect_input_files(sample_sheet):
    """
    The output of this function should
    be sufficient to run the workflow
    """
    sample_input = collections.defaultdict(dict)
    sample_lists = {
        "clr": [],
        "hifi": [],
        "ont": []
    }

    known_path_ids = dict()
    for row in sample_sheet.itertuples():
        sample = row.sample
        input_read_type = row.read_type
        sample_lists[input_read_type].append(sample)

        level_1_key = sample, input_read_type
        seq_input, input_hashes = collect_sequence_input(row.input)

        # need to distinguish Pacbio BAM
        # from reads already dumped to FASTQ
        file_map = {
            ("fastq", "all"): [],
            ("fastq", "path_hashes"): [],
            ("fastq", "path_ids"): [],
            ("bam", "all"): [],
            ("bam", "path_hashes"): [],
            ("bam", "path_ids"): [],

        }
        for seq_file, path_hash in zip(seq_input, input_hashes):
            # to improve readability, try using just a prefix
            # of full path hash
            path_id = path_hash[:8]
            if path_id in known_path_ids:
                assoc_path_hash = known_path_ids[path_id]
                if assoc_path_hash != path_hash:
                    # means: identical prefix, but different
                    # path hash, so cannot move on here
                    raise ValueError(
                        "Path prefix (ID) not unique:\n"
                        f"{assoc_path_hash}\n"
                        f"{path_hash}\n"
                    )
            known_path_ids[path_id] = path_hash
            if seq_file.name.lower().endswith("bam"):
                if input_read_type == "ont":
                    raise ValueError(
                        "Detected BAM sample input for read type ONT\n"
                        f"Ilegal input: {seq_file}\n"
                    )
                if "subread" in seq_file.name and input_read_type == "hifi":
                    raise ValueError(
                        f"BAM file appears to be a subreads file: {seq_file}\n"
                        "--> read type set to hifi in sample sheet."
                    )
                level_2_key = "bam", path_id
                file_map[level_2_key] = seq_file
                file_map[("bam", "path_ids")].append(path_id)
                file_map[("bam", "path_hashes")].append(path_hash)
                file_map[("bam", "all")].append(seq_file)

                fastq_file = f"{sample}_{input_read_type}.{path_id}.fastq.gz"
                dumped_fastq_path = DIR_PROC.joinpath(
                    "00-prepare", "dump_fastq", fastq_file
                )
                level_2_key = "fastq", path_id
                file_map[level_2_key] = dumped_fastq_path
                file_map[("fastq", "path_ids")].append(path_id)
                file_map[("fastq", "path_hashes")].append(path_hash)
                file_map[("fastq", "all")].append(dumped_fastq_path)
            else:  # assume FASTA/FASTQ etc.
                level_2_key = "fastq", path_id
                file_map[level_2_key] = seq_file
                file_map[("fastq", "path_ids")].append(path_id)
                file_map[("fastq", "path_hashes")].append(path_hash)
                file_map[("fastq", "all")].append(seq_file)

        sample_input[level_1_key] = file_map

    return sample_input, sample_lists


def _read_input_files_from_fofn(fofn_path):
    """Read input file listing from
    file of file names
    TODO: candidate for inclusion in template
    """

    input_files = []
    with open(fofn_path, "r") as listing:
        for line in listing:
            if not line.strip():
                continue
            try:
                file_path = pathlib.Path(line.strip()).resolve(strict=True)
            except FileNotFoundError:
                try:
                    file_path = DATA_ROOT.joinpath(line.strip()).resolve(strict=True)
                except FileNotFoundError:
                    err_msg = "\nERROR\n"
                    err_msg += f"Cannot find file: {line.strip}\n"
                    err_msg += f"Data root is set to: {DATA_ROOT}\n"
                    sys.stderr.write(err_msg)
                    raise
            input_files.append(file_path)

    return sorted(input_files)


def subset_path(full_path):
    """This helper exists to reduce
    the absolute path to a file
    to just the file name and its
    parent.
    TODO: should be codified as part
    of the template utilities to improve
    infrastructure portability of active
    workflows
    """
    folder_name = full_path.parent.name
    file_name = full_path.name
    subset_path = f"{folder_name}/{file_name}"
    # if it so happens that the file resides
    # in a root-level location, strip off
    # leading slash
    return subset_path.strip("/")


def collect_sequence_input(path_spec):
    """
    Generic function to collect HiFi or ONT/Nanopore
    input (read) files
    """
    input_files = []
    input_hashes = []
    for sub_input in path_spec.split(","):
        input_path = pathlib.Path(sub_input).resolve()
        if input_path.is_file() and input_path.name.endswith(".fofn"):
            fofn_files = _read_input_files_from_fofn(input_path)
            fofn_hashes = [
                hashlib.sha256(
                    subset_path(fp).encode("utf-8")
                ).hexdigest() for fp in fofn_files
            ]
            input_files.extend(fofn_files)
            input_hashes.extend(fofn_hashes)
        elif input_path.is_file():
            input_hash = hashlib.sha256(
                subset_path(input_path).encode("utf-8")
            ).hexdigest()
            input_files.append(input_path)
            input_hashes.append(input_hash)
        elif input_path.is_dir():
            collected_files = _collect_files(input_path)
            collected_hashes = [
                hashlib.sha256(
                    subset_path(fp).encode("utf-8")
                ).hexdigest() for fp in collected_files
            ]
            input_files.extend(collected_files)
            input_hashes.extend(collected_hashes)
        else:
            raise ValueError(f"Cannot handle input: {sub_input}")
    return input_files, input_hashes


def _build_constraint(values):
    escaped_values = sorted(map(re.escape, map(str, values)))
    constraint = "(" + "|".join(escaped_values) + ")"
    return constraint


def _collect_files(folder):

    all_files = set()
    for pattern in config["input_file_ext"]:
        pattern_files = set(folder.glob(f"**/*.{pattern}"))
        all_files = all_files.union(pattern_files)
    all_files = [f for f in sorted(all_files) if f.is_file()]
    if len(all_files) < 1:
        raise ValueError(f"No input files found underneath {folder}")
    return all_files

process_sample_sheet()
