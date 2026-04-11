.PHONY: \
	all test help clean cleanall fulgor_config\
	conda download download_asms download_mfur match aggregate_matches map \
	config report \
	cluster_slurm cluster_lsf cluster_lsf_test \
	format checkformat dryrun lint

SHELL := /usr/bin/env bash
.SHELLFLAGS := -eo pipefail -c
.DEFAULT_GOAL := all
DATETIME := $(shell date -u +"%Y_%m_%dT%H_%M_%S")

.SECONDARY:

.SUFFIXES:

THREADS := $(shell grep "^threads:" config.yaml | awk '{print $$2}')
MAX_DOWNLOAD_THREADS := $(shell grep "^max_download_threads" config.yaml | awk '{print $$2}')
DOWNLOAD_RETRIES := $(shell grep "^download_retries" config.yaml | awk '{print $$2}')
MAX_IO_HEAVY_THREADS := $(shell grep "^max_io_heavy_threads" config.yaml | awk '{print $$2}')
MAX_RAM_MB := $(shell grep "^max_ram_gb:" config.yaml | awk '{print $$2*1024}')
ifndef DATABASE
DATABASE := $(shell grep "^database:" config.yaml | awk '{print $$2}' | tr -d '"')
endif

ifeq ($(strip $(DATABASE)),)
DATABASE := 661k
endif

ifeq ($(DATABASE),ATB)
TEST_BATCHES := data/atb_batches_small.txt
TEST_EXPECTED := data/reads_1___reads_2___reads_3___reads_4-ATB.sam_summary.xz
else ifeq ($(DATABASE),661k)
TEST_BATCHES := data/batches_small.txt
TEST_EXPECTED := data/reads_1___reads_2___reads_3___reads_4.sam_summary.xz
else
$(error Unsupported DATABASE '$(DATABASE)'; expected 661k or ATB)
endif

SMK_DB_CFG := --config database=$(DATABASE)

ifeq ($(SMK_CLUSTER_ARGS),)
    # configure local run
    SMK_PARAMS := --cores ${THREADS} --rerun-incomplete --printshellcmds --keep-going --use-conda --resources max_download_threads=$(MAX_DOWNLOAD_THREADS) max_io_heavy_threads=$(MAX_IO_HEAVY_THREADS) max_ram_mb=$(MAX_RAM_MB)
else
    # configure cluster run
    SMK_PARAMS := --cores all --rerun-incomplete --printshellcmds --keep-going --use-conda --resources max_download_threads=10000000 max_io_heavy_threads=10000000 max_ram_mb=1000000000 $(SMK_CLUSTER_ARGS)
endif

DOWNLOAD_PARAMS := --cores $(MAX_DOWNLOAD_THREADS) -j $(MAX_DOWNLOAD_THREADS) --restart-times $(DOWNLOAD_RETRIES)
BENCHMARK_DIR := logs/benchmarks
REPORT_DIR := logs/reports
REPORT_FILE := $(REPORT_DIR)/report_$(DATABASE)_$(DATETIME).html


######################
## General commands ##
######################
all: ## Run everything (the default rule)
	$(MAKE) download DATABASE=$(DATABASE)
	$(MAKE) match DATABASE=$(DATABASE)
	$(MAKE) aggregate_matches DATABASE=$(DATABASE)
	$(MAKE) map DATABASE=$(DATABASE)

DIFF_CMD=diff -q <(gunzip --stdout output/reads_1___reads_2___reads_3___reads_4.sam_summary.gz | cut -f -3) <(xzcat $(TEST_EXPECTED) | cut -f -3)

test: ## Quick test using 3 batches
	snakemake download $(SMK_PARAMS) $(DOWNLOAD_PARAMS) $(SMK_DB_CFG) --config batches=$(TEST_BATCHES)  # download is not benchmarked
	mkdir -p $(BENCHMARK_DIR)
	scripts/benchmark.py --log $(BENCHMARK_DIR)/test_match_$(DATETIME).txt "snakemake match $(SMK_PARAMS) $(SMK_DB_CFG) --config batches=$(TEST_BATCHES) nb_best_hits=1"
	scripts/benchmark.py --log $(BENCHMARK_DIR)/test_map_$(DATETIME).txt   "snakemake map $(SMK_PARAMS) $(SMK_DB_CFG) --config batches=$(TEST_BATCHES) nb_best_hits=1"
	@if $(DIFF_CMD); then \
	    echo "Success! Test run produced the expected output."; \
	else \
		echo ""; \
		echo "ERROR. Test run DID NOT produce the expected output. Failed command:"; \
		echo ""; \
		echo "    $(DIFF_CMD)"; \
		echo ""; \
	    exit 1;\
	fi

help: ## Print help messages
	@echo -e "$$(grep -hE '^\S*(:.*)?##' $(MAKEFILE_LIST) \
		| sed \
			-e 's/:.*##\s*/:/' \
			-e 's/^\(.*\):\(.*\)/   \\x1b[36m\1\\x1b[m:\2/' \
			-e 's/^\([^#]\)/\1/g' \
			-e 's/: /:/g' \
			-e 's/^#\(.*\)#/\\x1b[90m\1\\x1b[m/' \
		| column -c2 -t -s : )"

dryrun: ## Show the planned Snakemake workflow without executing jobs
	snakemake --dry-run $(SMK_PARAMS) $(SMK_DB_CFG)

clean: ## Clean intermediate search files
	rm -fv intermediate/*/*
	rm -rfv logs
	rm -fv output/*
	mkdir -p .snakemake/old_log
	mv -v .snakemake/log/*.log .snakemake/old_log/ || true

cleanall: clean ## Clean all generated and downloaded files
	rm -f asms/*.xz{,.tmp}
	rm -f mfur/*.mfur

fulgor_config: ##Install Fulgor dependencies and compile
	git submodule update --init --recursive
	snakemake fulgor_config $(SMK_PARAMS)

####################
## Pipeline steps ##
####################

conda: ## Create the conda environments
	snakemake $(SMK_PARAMS) $(SMK_DB_CFG) --conda-create-envs-only

download: ## Download the assemblies and meta-Fulgor indexes
	snakemake download $(SMK_PARAMS) $(DOWNLOAD_PARAMS) $(SMK_DB_CFG)

download_asms: ## Download only the assemblies
	snakemake download_asms_batches $(SMK_PARAMS) $(DOWNLOAD_PARAMS) $(SMK_DB_CFG)

download_mfur: ## Download only the meta-Fulgor indexes
	snakemake download_mfur_batches $(SMK_PARAMS) $(DOWNLOAD_PARAMS) $(SMK_DB_CFG)

match: ## Match queries using Fulgor and select the best candidates per batch (queries -> candidates per batch)
	mkdir -p $(BENCHMARK_DIR)
	scripts/benchmark.py --log $(BENCHMARK_DIR)/match_$(DATETIME).txt "snakemake match $(SMK_PARAMS) $(SMK_DB_CFG)"

aggregate_matches: ## Select the best candidates across the entire reference collection (candidates per batch -> overall candidates)
	mkdir -p $(BENCHMARK_DIR)
	scripts/benchmark.py --log $(BENCHMARK_DIR)/aggregated_match_$(DATETIME).txt "snakemake aggregate_matches $(SMK_PARAMS) $(SMK_DB_CFG)"

map: ## Map candidates to assemblies (overall candidates -> alignments)
	mkdir -p $(BENCHMARK_DIR)
	scripts/benchmark.py --log $(BENCHMARK_DIR)/map_$(DATETIME).txt   "snakemake map $(SMK_PARAMS) $(SMK_DB_CFG)"

###############
## Reporting ##
###############

config: ## Print configuration without comments
	@cat config.yaml \
		| perl -pe 's/ *#.*//g' \
		| grep --color='auto' -E '.*\:'
	@#| grep -Ev ^$$

report: ## Generate Snakemake report
	mkdir -p $(REPORT_DIR)
	snakemake $(SMK_DB_CFG) --report $(REPORT_FILE)



#############
## Cluster ##
#############

cluster_slurm: ## Submit to a SLURM cluster
	sbatch \
        -c 10 \
        --mem=80GB \
        -t 0-08:00:00 \
        --wrap="make DATABASE=$(DATABASE)"

cluster_lsf: ## Submit to LSF cluster
	scripts/check_if_config_is_ok_for_cluster_run.py
	scripts/submit_lsf.sh DATABASE=$(DATABASE)

cluster_lsf_test: ## Submit the test pipeline to LSF cluster
	scripts/check_if_config_is_ok_for_cluster_run.py
	scripts/submit_lsf.sh test DATABASE=$(DATABASE)


####################
## For developers ##
####################

format: ## Reformat Python and Snakemake files
	yapf -i */*.py
	snakefmt Snakefile

checkformat: ## Check source code format
	snakefmt --check Snakefile
	yapf --diff */*.py

lint: ## Lint the Snakemake workflow
	snakemake --lint $(SMK_PARAMS) $(SMK_DB_CFG) || true
