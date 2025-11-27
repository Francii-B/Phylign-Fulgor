#! /usr/bin/env python3

import argparse
import gzip
from concurrent.futures import ProcessPoolExecutor


def open_maybe_gzip(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def process_match_file(match_fn):
    queries_in_file = set()
    matched_queries_in_file = set()
    genomes_in_file = set()
    pair_count = 0
    has_hit = False
    qname = None

    with open_maybe_gzip(match_fn) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("*"):
                parts = line[1:].split("\t")
                qname = parts[0].split(" ")[0]
                queries_in_file.add(qname)
            else:
                tmp_name, _ = line.split()
                _, ref = tmp_name.split("_", 1)

                has_hit = True
                pair_count += 1
                genomes_in_file.add(ref)
                matched_queries_in_file.add(qname)

    return {
        "queries_in_file": queries_in_file,
        "matched_queries_in_file": matched_queries_in_file,
        "genomes_count": len(genomes_in_file),
        "pair_count": pair_count,
        "has_hit": has_hit,
    }


def main():
    parser = argparse.ArgumentParser(description="Summarize COBS-like outputs.")
    parser.add_argument("match_fn", nargs="+", help="COBS match files (gzipped or plain).")
    parser.add_argument("-j", dest="n_jobs", type=int, default=8, help="Number of worker processes [8].")
    args = parser.parse_args()

    n_batches = len(args.match_fn)

    with ProcessPoolExecutor(max_workers=args.n_jobs) as ex:
        results = list(ex.map(process_match_file, args.match_fn))

    all_queries = set()
    matched_queries = set()
    total_pairs = 0
    total_target_genomes = 0
    total_target_batches = 0

    for r in results:
        all_queries.update(r["queries_in_file"])
        matched_queries.update(r["matched_queries_in_file"])
        total_pairs += r["pair_count"]
        total_target_genomes += r["genomes_count"]
        if r["has_hit"]:
            total_target_batches += 1

    print(f"queries\t{len(all_queries)}")
    print(f"matched_queries\t{len(matched_queries)}")
    print(f"distinct_genome_query_pairs\t{total_pairs}")
    print(f"target_genomes\t{total_target_genomes}")
    print(f"target_batches\t{total_target_batches}")
    print(f"processed_batches\t{n_batches}")


if __name__ == "__main__":
    main()
