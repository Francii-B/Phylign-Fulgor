def multiglob(patterns):
    files = []
    for pattern in patterns:
        files.extend(glob.glob(pattern))
    files = list(map(Path, files))
    return files


def get_all_query_filepaths():
    return multiglob(expand("input/*.{ext}", ext=extensions))


def get_all_query_filenames():
    return sorted([file.with_suffix("").name for file in get_all_query_filepaths()])


def get_batches():
    with open(BATCHES_FILE) as fin:
        return list(sorted(filter(len, map(str.strip, fin))))


def get_filename_for_all_queries():
    return "___".join(get_all_query_filenames())


def get_index_metadata(wildcards, input):
    batch = wildcards.batch
    decompressed_indexes_sizes_filepath = input.decompressed_indexes_sizes
    with open(decompressed_indexes_sizes_filepath) as decompressed_indexes_sizes_fh:
        for line in decompressed_indexes_sizes_fh:
            index_path, size_in_bytes, xz_decompress_RAM = line.strip().split()
            batch_for_index = Path(index_path).stem
            size_in_bytes = int(size_in_bytes)
            xz_decompress_RAM = int(xz_decompress_RAM)
            if batch == batch_for_index:
                return size_in_bytes, xz_decompress_RAM

    assert (
        False
    ), f"Error getting uncompressed batch size for batch {batch}: batch not found"


def get_uncompressed_batch_size(wildcards, input):
    return get_index_metadata(wildcards, input)[0]


def get_xz_decompress_RAM_in_MB(wildcards, input):
    xz_decompression_RAM_usage_in_bytes = get_index_metadata(wildcards, input)[1]
    xz_decompression_RAM_usage_in_MB = (
        int(xz_decompression_RAM_usage_in_bytes / 1024 / 1024) + 1
    )
    return xz_decompression_RAM_usage_in_MB


def get_uncompressed_batch_size_in_MB(wildcards, input, ignore_RAM, streaming):
    # if ignore_RAM:
    #    return 0
    # if streaming:
    #    # then we are decompressing and matching at the same time
    #    xz_decompression_RAM_usage_in_MB = get_xz_decompress_RAM_in_MB(wildcards, input)
    # else:
    # xz_decompression_RAM_usage_in_MB = 0
    size_in_bytes = get_uncompressed_batch_size(wildcards, input)
    size_in_MB = int(size_in_bytes / 1024 / 1024) + 1
    return size_in_MB  # + xz_decompression_RAM_usage_in_MB


def get_max_number_of_fulgor_threads_from_auto_string(auto_string):
    fulgor_threads = re.findall(r"auto\((\d+)\)", auto_string)
    parsing_was_successful = len(fulgor_threads) == 1
    assert parsing_was_successful, "Error parsing parameter fulgor_threads"
    fulgor_threads = int(fulgor_threads[0])
    return fulgor_threads


def get_number_of_fulgor_threads(wildcards, input, predefined_fulgor_threads, streaming):
    user_defined_nb_of_threads = not predefined_fulgor_threads.startswith("auto")
    if user_defined_nb_of_threads:
        return int(predefined_fulgor_threads)

    use_max_cores = predefined_fulgor_threads == "auto"
    if use_max_cores:
        max_number_of_fulgor_threads = workflow.cores
    else:
        max_number_of_fulgor_threads = get_max_number_of_fulgor_threads_from_auto_string(
            predefined_fulgor_threads
        )

    uncompressed_batch_size_in_MB = get_uncompressed_batch_size_in_MB(
        wildcards, input, ignore_RAM=False, streaming=streaming
    )
    max_RAM_MB = int(config["max_ram_gb"]) * 1024
    number_of_cores_to_use = round(
        uncompressed_batch_size_in_MB / max_RAM_MB * max_number_of_fulgor_threads
    )
    number_of_cores_to_use = max(number_of_cores_to_use, 1)
    number_of_cores_to_use = min(number_of_cores_to_use, max_number_of_fulgor_threads)
    is_using_more_than_half_of_the_cores = (
        number_of_cores_to_use > max_number_of_fulgor_threads / 2
    )
    if is_using_more_than_half_of_the_cores:
        # Usually in this situation we run just one Fulgor job simultaneously.
        number_of_cores_to_use = max_number_of_fulgor_threads
    return number_of_cores_to_use


def get_index_load_mode():
    allowed_index_load_modes = ["mem-stream", "mem-disk", "mmap-disk"]
    index_load_mode = config["index_load_mode"]
    assert (
        index_load_mode in allowed_index_load_modes
    ), f"index_load_mode must be one of {allowed_index_load_modes}"
    return index_load_mode


def get_atb_asms_urls():
    with open("data/atb_asms_download_link.tsv") as fin:
        return dict(line.strip().split("\t", 1) for line in fin)


def mfur_url_661k(batch):
    if batch <= "dustbin__15":
        return f"https://zenodo.org/record/14002973/files/{batch}.mfur"
    elif batch <= "mycobacterium_kansasii__01":
        return f"https://zenodo.org/record/14002975/files/{batch}.mfur"
    elif batch <= "salmonella_enterica__33":
        return f"https://zenodo.org/record/14006705/files/{batch}.mfur"
    else:
        return f"https://zenodo.org/record/14006707/files/{batch}.mfur"


def mfur_url_atb(batch):
    batch_id = int(batch.split(".")[-1])
    if batch_id <= 100:
        return f"https://zenodo.org/record/15994164/files/{batch}.mfur"
    elif batch_id <= 134:
        return f"https://zenodo.org/record/15994228/files/{batch}.mfur"
    elif batch_id <= 225:
        return f"https://zenodo.org/record/15994270/files/{batch}.mfur"
    elif batch_id <= 325:
        return f"https://zenodo.org/record/15994318/files/{batch}.mfur"
    elif batch_id <= 425:
        return f"https://zenodo.org/record/15994445/files/{batch}.mfur"
    elif batch_id <= 525:
        return f"https://zenodo.org/record/15994505/files/{batch}.mfur"
    elif batch_id <= 625:
        return f"https://zenodo.org/record/15994553/files/{batch}.mfur"
    else:
        return f"https://zenodo.org/record/15994624/files/{batch}.mfur"


def mfur_url_fct(wildcards):
    if DB == "661k":
        return mfur_url_661k(wildcards.batch)
    return mfur_url_atb(wildcards.batch)


def asms_url_661k(batch):
    asm_zenodo = 4602622
    asm_url = f"https://zenodo.org/record/{asm_zenodo}/files/{batch}.tar.xz"
    return asm_url


def asms_url_atb(batch):
    key = f"{batch}.tar.xz"
    assert key in ATB_ASMS_URLS, f"Missing ATB assembly URL for batch {batch}"
    return ATB_ASMS_URLS[key]


def asms_url_fct(wildcards):
    if DB == "661k":
        return asms_url_661k(wildcards.batch)
    return asms_url_atb(wildcards.batch)


def get_sleep_amount(attempt):
    return int(config["download_retry_wait"]) * (attempt - 1)


def get_query_file(wildcards):
    query_file = multiglob(expand(f"input/{wildcards.qfile}.{{ext}}", ext=extensions))
    assert len(query_file) == 1
    return query_file[0]
