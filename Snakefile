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
include: "rules/common.smk"


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
                    -o {output.raw_fulgor_output} \\
                 && \\
            ./scripts/postprocess_kmer_matches.py -n {params.nb_best_hits} < {output.raw_fulgor_output} \\
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
