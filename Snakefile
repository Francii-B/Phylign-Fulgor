import subprocess
import sys
from pathlib import Path
from snakemake.utils import min_version


##################################
## checks that gcc-11 is available for darwin
##################################
if sys.platform=="darwin":
    try:
        subprocess.check_call(["gcc-11", "--version"])
        subprocess.check_call(["g++-11", "--version"])
    except subprocess.CalledProcessError:
        print("Error: you are running on Mac OS X and gcc-11 or g++-11 were not detected. Try installing with brew install gcc@11")
        sys.exit(1)


##################################
## Initialization
##################################


configfile: "config.yaml"


min_version("6.2.0")
shell.prefix("set -euo pipefail")

batches = [x.strip() for x in open(config["batches"])]
print(f"Batches: {batches}")

qfiles = [x.with_suffix("").name for x in Path("queries").glob("*.fa")]
print(f"Query files: {qfiles}")


wildcard_constraints:
    batch=".+__\d\d",


##################################
## Download params
##################################
def cobs_url(wildcards):
    x = wildcards.batch
    if x >= "eubacterium":
        return f"https://zenodo.org/record/6849657/files/{x}.cobs_classic.xz"
    else:
        return f"https://zenodo.org/record/6845083/files/{x}.cobs_classic.xz"


asm_zenodo = 4602622
asms_url = f"https://zenodo.org/record/{asm_zenodo}/files"


##################################
## Top-level rules
##################################


rule all:
    """Run all
    """
    input:
        [f"output/{qfile}.sam_summary.xz" for qfile in qfiles],


rule download:
    """Download assemblies and COBS indexes.
    """
    input:
        [f"asms/{x}.tar.xz" for x in batches],
        [f"cobs/{x}.cobs_classic.xz" for x in batches],


rule match:
    """Match reads to the COBS indexes.
    """
    input:
        [f"intermediate/02_filter/{qfile}.fa" for qfile in qfiles],


rule map:
    """Map reads to the assemblies.
    """
    input:
        [f"output/{qfile}.sam_summary.xz" for qfile in qfiles],


##################################
## Download rules
##################################
rule download_asm_batch:
    """Download compressed assemblies
    """
    output:
        xz="asms/{batch}.tar.xz",
    params:
        url=asms_url,
    resources:
        download_thr=1,
    threads: 1
    shell:
        """
        curl "{params.url}/{wildcards.batch}.tar.xz"  > {output.xz}
        scripts/test_xz.py {output.xz}
        """


rule download_cobs_batch:
    """Download compressed cobs indexes
    """
    output:
        xz="cobs/{batch}.cobs_classic.xz",
    params:
        url=cobs_url,
    resources:
        download_thr=1,
    threads: 1
    shell:
        """
        curl "{params.url}"  > {output.xz}
        scripts/test_xz.py {output.xz}
        """


##################################
## Installation rules
##################################
rule install_cobs:
    output:
        "tools/cobs",
    threads: 1
    conda:
        f"envs/compile_cobs_{sys.platform}.yaml"
    shadow: "shallow"
    params:
        cobs_version = "0.2.0"
    log:
        "logs/install_cobs.log"
    shell: """
        ./scripts/install_cobs.sh {params.cobs_version} {output} >{log} 2>&1
    """


##################################
## Processing rules
##################################
rule decompress_cobs:
    """Decompress cobs indexes
    """
    output:
        cobs=temp("intermediate/00_cobs/{batch}.cobs_classic"),
    input:
        xz="cobs/{batch}.cobs_classic.xz",
    resources:
        decomp_thr=1,
    threads: config["cobs_thr"]  # The same number as of COBS threads to ensure that COBS is executed immediately after decompression
    shell:
        """
        xzcat "{input.xz}" > "{output.cobs}"
        """


rule fix_query:
    """Fix query to expected COBS format: single line fastas composed of ACGT bases only
    """
    output:
        fixed_query="intermediate/fixed_queries/{qfile}.fa",
    input:
        original_query="queries/{qfile}.fa",
    threads: 1
    conda:
        "envs/seqtk.yaml"
    params:
        base_to_replace="A"
    shell:
        """
        seqtk seq -A -U {input.original_query} | \
        awk '{{if(NR%2==1){{print $0;}}else{{gsub(/[^ACGT]/, \"{params.base_to_replace}\"); print;}}}}' \
        > {output.fixed_query}
        """


rule run_cobs:
    """Cobs matching
    """
    output:
        match=protected("intermediate/01_match/{batch}____{qfile}.xz"),
    input:
        cobs_executable = rules.install_cobs.output,
        cobs_index="intermediate/00_cobs/{batch}.cobs_classic",
        fa=rules.fix_query.output.fixed_query,
    threads: config["cobs_thr"]  # Small number in order to guarantee Snakemake parallelization
    params:
        kmer_thres=config["cobs_kmer_thres"],
    priority: 999
    shell:
        """
        {input.cobs_executable} query \\
            -t {params.kmer_thres} \\
            -T {threads} \\
            -i {input.cobs_index} \\
            -f {input.fa} \\
        | xz -v \\
        > {output.match}
        """


rule translate_matches:
    """Translate cobs matches.

    Output:
        ref - read - matches
    """
    output:
        fa="intermediate/02_filter/{qfile}.fa",
    input:
        fa=rules.fix_query.output.fixed_query,
        all_matches=[
            f"intermediate/01_match/{batch}____{{qfile}}.xz" for batch in batches
        ],
    conda:
        "envs/minimap2.yaml"
    threads: 1
    shell:
        """
        ./scripts/filter_queries.py -q {input.fa} {input.all_matches} \\
            > {output.fa}
        """


rule batch_align_minimap2:
    output:
        sam="intermediate/03_map/{batch}____{qfile}.sam",
    input:
        qfa="intermediate/02_filter/{qfile}.fa",
        asm="asms/{batch}.tar.xz",
    log:
        log="logs/03_map/{batch}____{qfile}.log",
    params:
        minimap_preset=config["minimap_preset"],
    conda:
        "envs/minimap2.yaml"
    threads: 1
    shell:
        """
        ((
        ./scripts/batch_align.py \\
            --minimap-preset {params.minimap_preset} \\
            {input.asm} \\
            {input.qfa} \\
            > {output.sam}
        ) 2>&1 ) > {log}
        """


rule aggregate_sams:
    output:
        pseudosam="output/{qfile}.sam_summary.xz",
    input:
        sam=[f"intermediate/03_map/{batch}____{{qfile}}.sam" for batch in batches],
    threads: 1
    shell:
        """
        head -n 9999999 {input.sam} \\
            | grep -v "@" \\
            | xz \\
            > {output.pseudosam}
        """
