import functools
import glob
from pathlib import Path

from snakemake.utils import min_version
import random
import re

##################################
## Helper functions
##################################


extensions = ["fa", "fasta", "fq", "fastq"]


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


##################################
## Initialization
##################################


configfile: "config.yaml"


min_version("6.2.0")
shell.prefix("set -euo pipefail")

DB_CONF = {
    "661k": {
        "default_batches": "data/batches_full.txt",
        "batch_regex": r".+__\d\d",
        "decompressed_indexes_sizes": "data/decompressed_indexes_sizes.txt",
        "accessions": "data/661k_batches.txt.xz",
    },
    "ATB": {
        "default_batches": "data/atb_batches_full.txt",
        "batch_regex": r".+\.batch\.\d{1,3}",
        "asm_url_tsv": "data/atb_asms_download_link.tsv",
        "decompressed_indexes_sizes": "data/decompressed_indexes_sizes_fulgor_ATB.txt",
        "accessions": "data/ATB_batches.txt.xz",
    },
}

DB = config.get("database", "661k")
assert DB in DB_CONF, f"database must be one of {sorted(DB_CONF)}"
BATCHES_FILE = config.get("batches", DB_CONF[DB]["default_batches"])
ATB_ASMS_URLS = get_atb_asms_urls() if DB == "ATB" else {}

batches = get_batches()
print(f"Batches: {batches}")

qfiles = get_all_query_filepaths()
print(f"Query files: {list(map(str, qfiles))}")

assemblies_dir = Path(f"{config['download_dir']}/asms")
mfur_dir = Path(f"{config['download_dir']}/mfur")

predefined_fulgor_threads = str(
    config["fulgor_threads"]
)

ignore_RAM = False
streaming = False
fulgor_is_an_io_heavy_job = False
index_load_mode = get_index_load_mode()

if index_load_mode == "mem-stream":
    streaming = True
elif index_load_mode == "mmap-disk":
    # we ignore RAM usage because the OS is responsible for controlling RAM usage in this case
    ignore_RAM = True
    # Fulgor becomes IO-heavy in mmap mode because the index is accessed on demand.
    fulgor_is_an_io_heavy_job = True


wildcard_constraints:
    batch=DB_CONF[DB]["batch_regex"],


##################################
## Top-level rules
##################################
rule all:
    """Run all
    """
    input:
        f"output/{get_filename_for_all_queries()}.sam_summary.gz",
        f"output/{get_filename_for_all_queries()}.sam_summary.stats",


rule download:
    """Download assemblies and meta-Fulgor indexes.
    """
    input:
        [f"{assemblies_dir}/{x}.tar.xz" for x in batches],
        [f"{mfur_dir}/{x}.mfur" for x in batches],


rule download_asms_batches:
    """Download assemblies.
    """
    input:
        [f"{assemblies_dir}/{x}.tar.xz" for x in batches],


rule download_mfur_batches:
    """Download meta-Fulgor indexes.
    """
    input:
        [f"{mfur_dir}/{x}.mfur" for x in batches],


rule match:
    """Match reads to the meta-Fulgor indexes.
    """
    input:
        all_matches=[
            f"intermediate/03_match/{batch}____{get_filename_for_all_queries()}.gz"
            for batch in batches
        ],


rule aggregate_matches:
    """Match reads to the meta-Fulgor indexes + aggregate the results
    """
    input:
        f"intermediate/04_filter/{get_filename_for_all_queries()}.fa",


rule map:
    """Map reads to the assemblies.
    """
    input:
        f"output/{get_filename_for_all_queries()}.sam_summary.gz",
        f"output/{get_filename_for_all_queries()}.sam_summary.stats",


rule fulgor_config:
    """Install Fulgor dependencies and compile
    """
    conda:
        "envs/fulgor_env.yml"
    shell:
        """
        ./scripts/submodule.sh
        """


##################################
## Download rules
##################################
rule download_asm_batch:
    """Download compressed assemblies
    """
    output:
        xz=f"{assemblies_dir}/{{batch}}.tar.xz",
    threads: 1
    resources:
        max_download_threads=1,
        mem_mb=200,
        # note: sleep_amount has to be defined as a resource
        # note: I tried a hack to route it to params, but it did not work, see https://github.com/snakemake/snakemake/issues/499
        sleep_amount=lambda wildcards, attempt: get_sleep_amount(attempt),
    params:
        url=asms_url_fct,
    shell:
        """
        scripts/download.sh {params.url} {output.xz} {resources.sleep_amount} .tar.xz
        """


rule download_mfur_batch:
    """Download uncompressed meta-Fulgor indexes
    """
    output:
        batch=f"{mfur_dir}/{{batch}}.mfur",
    threads: 1
    resources:
        max_download_threads=1,
        mem_mb=200,
        sleep_amount=lambda wildcards, attempt: get_sleep_amount(attempt),
    params:
        url=mfur_url_fct,
    shell:
        """
        scripts/download.sh {params.url} {output.batch} {resources.sleep_amount} .mfur
        """


##################################
## Processing rules
##################################
def get_query_file(wildcards):
    query_file = multiglob(expand(f"input/{wildcards.qfile}.{{ext}}", ext=extensions))
    assert len(query_file) == 1
    return query_file[0]


rule fix_query:
    """Normalize query to the matching input format: single-line FASTA with ACGT bases only.
    """
    output:
        fixed_query="intermediate/00_queries_preprocessed/{qfile}.fa",
    input:
        original_query=get_query_file,
    threads: 1
    resources:
        mem_mb=200,
    conda:
        "envs/seqtk.yaml"
    params:
        base_to_replace="A",
    shell:
        """
        seqtk seq -A -U -C {input.original_query} \\
                | awk '{{if(NR%2==1){{print $0;}}else{{gsub(/[^ACGT]/, \"{params.base_to_replace}\"); print;}}}}' \\
            > {output.fixed_query}
        """


rule concatenate_queries:
    """Concatenate queries so we run matching and alignment only once per batch.
    """
    output:
        concatenated_query=f"intermediate/01_queries_merged/{get_filename_for_all_queries()}.fa",
    input:
        all_queries=expand(
            "intermediate/00_queries_preprocessed/{qfile}.fa",
            qfile=get_all_query_filenames(),
        ),
    threads: 1
    resources:
        mem_mb=200,
    shell:
        """
        cat {input} > {output}
        """


# note: snakefmt makes incorrect breaks and spacing for threads; to keep the lines
#       short to prevent this behaviour, we use the following function
partial_fulgor_threads = functools.partial(
    get_number_of_fulgor_threads,
    predefined_fulgor_threads=predefined_fulgor_threads,
    streaming=streaming,
)


rule run_fulgor:
    """Run Fulgor for the matching stage against meta-Fulgor indexes.
    """
    output:
        raw_fulgor_output=temp("intermediate/03_match/{batch}____{qfile}-preprocessed.tsv"),
        match="intermediate/03_match/{batch}____{qfile}.gz",
    input:
        mfur_index=f"{mfur_dir}/{{batch}}.mfur",
        fa="intermediate/01_queries_merged/{qfile}.fa",
        decompressed_indexes_sizes=DB_CONF[DB]["decompressed_indexes_sizes"],
    resources:
        max_io_heavy_threads=int(fulgor_is_an_io_heavy_job),
        max_ram_mb=lambda wildcards, input: get_uncompressed_batch_size_in_MB(
            wildcards, input, ignore_RAM, streaming
        ),
        mem_mb=lambda wildcards, input: int(
            get_uncompressed_batch_size_in_MB(wildcards, input, ignore_RAM, streaming)
            + 1024
        ),
    threads: partial_fulgor_threads
    params:
        fulgor_threshold=config["fulgor_threshold"],
        nb_best_hits=config["nb_best_hits"],
    priority: 999
    # modified-Fulgor emits COBS-compatible text consumed by a downstream compatibility layer.
    shell:
        """
        ./scripts/benchmark.py --log logs/benchmarks/run_fulgor/{wildcards.batch}____{wildcards.qfile}.txt \\
            './external/modified-Fulgor/build/fulgor pseudoalign \\
                    --threshold {params.fulgor_threshold} \\
                    -t {threads} \\
                    -i {input.mfur_index} \\
                    -q {input.fa} --cobs \\
                    -o {output.raw_fulgor_output}; \\
                cat {output.raw_fulgor_output} | ./scripts/postprocess_kmer_matches.py -n {params.nb_best_hits} \\
                | gzip --fast \\
                > {output.match}'
        """


rule translate_matches:
    """Translate mfur matches.

    Output:
        ref - read - matches
    """
    output:
        fa="intermediate/04_filter/{qfile}.fa",
    input:
        fa="intermediate/01_queries_merged/{qfile}.fa",
        all_matches=[
            f"intermediate/03_match/{batch}____{{qfile}}.gz" for batch in batches
        ],
    conda:
        "envs/minimap2.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * 2 ** (attempt),  # 4GB, 8GB, 16GB, 32GB...
    log:
        "logs/04_filter/{qfile}.log",
    params:
        nb_best_hits=config["nb_best_hits"],
    shell:
        """
        ./scripts/benchmark.py --log logs/benchmarks/translate_matches/translate_matches___{wildcards.qfile}.txt \\
            './scripts/filter_queries.py \\
                    -n {params.nb_best_hits} \\
                    -q {input.fa} \\
                    {input.all_matches} \\
                > {output.fa} 2>{log}'
        """


rule batch_align_minimap2:
    output:
        sam="intermediate/05_map/{batch}____{qfile}.sam.gz",
    input:
        qfa="intermediate/04_filter/{qfile}.fa",
        asm=f"{assemblies_dir}/{{batch}}.tar.xz",
    log:
        log="logs/05_map/{batch}____{qfile}.log",
    params:
        minimap_preset=config["minimap_preset"],
        minimap_extra_params=config["minimap_extra_params"],
        pipe="--pipe" if config["prefer_pipe"] else "",
        refs_tmp="intermediate/05_map/{batch}____{qfile}.refs.tmp",
        accessions=DB_CONF[DB]["accessions"],
    conda:
        "envs/minimap2.yaml"
    threads: config["minimap_threads"]
    resources:
        mem_mb=lambda wildcards, attempt: 1000 * 2 ** (attempt),  # 1GB, 2GB, 4GB, 8GB...
    shell:
        """
        xzcat {params.accessions} \\
            | grep {wildcards.batch} \\
            | cut -f2 \\
            > {params.refs_tmp}

        ./scripts/benchmark.py --log logs/benchmarks/batch_align_minimap2/{wildcards.batch}____{wildcards.qfile}.txt \\
            './scripts/batch_align.py \\
                    --minimap-preset {params.minimap_preset} \\
                    --threads {threads} \\
                    --extra-params=\"{params.minimap_extra_params}\" \\
                    --accessions {params.refs_tmp} \\
                    {params.pipe} \\
                    {input.asm} \\
                    {input.qfa} \\
                2>{log} \\
                | {{ grep -Ev "^@" || true; }} \\
                | gzip --fast\\
                > {output.sam}'

        rm -f {params.refs_tmp}
        """


rule aggregate_sams:
    output:
        pseudosam="output/{qfile}.sam_summary.gz",
    input:
        sam=[f"intermediate/05_map/{batch}____{{qfile}}.sam.gz" for batch in batches],
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1000 * 2 ** (attempt),  # 1GB, 2GB, 4GB, 8GB...
    shell:
        """
        ./scripts/benchmark.py --log logs/benchmarks/aggregate_sams/aggregate_sams___{wildcards.qfile}.txt \\
            './scripts/aggregate_sams.sh {input.sam} \\
                > {output.pseudosam}'
        """


rule final_stats:
    output:
        stats="output/{qfile}.sam_summary.stats",
    input:
        pseudosam="output/{qfile}.sam_summary.gz",
        concatenated_query=f"intermediate/01_queries_merged/{get_filename_for_all_queries()}.fa",
    conda:
        "envs/minimap2.yaml"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1000 * 2 ** (attempt),  # 1GB, 2GB, 4GB, 8GB...
    shell:
        """
        ./scripts/benchmark.py --log logs/benchmarks/aggregate_sams/final_stats___{wildcards.qfile}.txt \\
            './scripts/final_stats.py {input.concatenated_query} {input.pseudosam} \\
                > {output.stats}'
        """
