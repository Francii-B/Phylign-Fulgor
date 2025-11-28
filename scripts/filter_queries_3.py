#! /usr/bin/env python3

import argparse
import atexit
import collections
import os
import re
import sys

from pathlib import Path
from pprint import pprint

from xopen import xopen

from concurrent.futures import ProcessPoolExecutor, as_completed

DEFAULT_KEEP = 0
"""
For every read we may want to know top N matches (N = keep_matches).

Approach:
    - for every read keep a buffer in memory
    - iterate over all translated cobs outputs
    - keep just top k scores when keep_matches > 0
    - when keep_matches == 0, keep all matches (no filtering)
"""


def cobs_iterator(cobs_matches_fn):
    """Iterator for cobs matches.

    Assumes that cobs ref names start with a random sorting prefix followed by
    an underscore (embedded by Leandro).

    Args:
        cobs_matches_fn (str): File name of cobs output.

    Returns:
        (qname, batch, matches): Qname, batch id, and list of assignments of
        the same query, in the form (ref, kmers)

    Todo:
        - if necessary in the future, add batch name from the file name
    """
    qname = None
    matches_buffer = []
    batch = os.path.basename(cobs_matches_fn).split("____")[0]
    print(f"Translating matches {cobs_matches_fn}", file=sys.stderr)
    with xopen(cobs_matches_fn) as f:
        for x in f:
            x = x.strip()
            if not x:
                continue
            if x[0] == "*":
                # HEADER
                # empty buffer
                if qname is not None:
                    yield qname, batch, matches_buffer
                    matches_buffer = []
                # parse header
                parts = x[1:].split("\t")
                qname = parts[0].split(" ")[0]  # remove fasta comments
                nmatches = int(parts[1])
            else:
                # MATCH
                tmp_name, kmers = x.split()
                rid, ref = tmp_name.split("_")
                matches_buffer.append((ref, kmers))
    # flush last buffer
    if qname is not None:
        yield qname, batch, matches_buffer


def readfq(fp):  # this is a generator function
    # From https://github.com/lh3/readfq/blob/master/readfq.py
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs)
                    # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break


class SingleQuery:
    """A simple optimized buffer for keeping top matches for a single read across batches.

    Args:
        keep_matches (int): The number of top matches to keep.
                            If 0, keep all matches (no filtering).

    Attributes:
        _matches (list): A list of (batch, ref, kmers)
    """

    def __init__(self, qname, seq, keep_matches):
        self._keep_matches = keep_matches
        # should be increased once the number of records > keep
        self._min_matching_kmers = 0
        self._matches = []
        self._qname = qname
        self._seq = seq

    def add_matches(self, batch, matches):
        """Add matches."""
        for mtch in matches:
            ref, kmers = mtch
            kmers = int(kmers)
            if kmers >= self._min_matching_kmers:
                self._matches.append((batch, ref, kmers))
        self._housekeeping()

    def _housekeeping(self):
        # When keep_matches <= 0, do not filter at all: keep everything.
        if self._keep_matches is None or self._keep_matches <= 0:
            return

        # 1. sort
        self._matches.sort(key=lambda x: (-x[2], x[0], x[1]))

        # 2. separate losers
        losers = self._matches[self._keep_matches:]
        self._matches = self._matches[:self._keep_matches]

        # 3. return back tie records from losers
        if losers:
            # get the tie value
            self._min_matching_kmers = self._matches[-1][2]
            for x in losers:
                if x[2] == self._min_matching_kmers:
                    self._matches.append(x)
                else:
                    break

    def fasta_record_matches(self):
        name = self._qname
        com = ",".join([x[1] for x in self._matches])
        seq = self._seq
        return f">{name} {com}\n{seq}"


class Sift:
    """Sifting class for all reported cobs assignments."""

    def __init__(self, query_fn, keep_matches):
        self._query_fn = query_fn
        self._keep_matches = keep_matches
        self._query_dict = self._create_query_dict(query_fn, keep_matches)

    @staticmethod
    def _create_query_dict(fx_fn, keep_matches):
        d = collections.OrderedDict()
        with xopen(fx_fn) as fo:
            for qname, seq, _ in readfq(fo):
                d[qname] = SingleQuery(qname=qname, keep_matches=keep_matches, seq=seq)
        return d

    def process_cobs_file(self, cobs_fn):
        """Sequential processing for a single cobs file (kept for compatibility / debugging)."""
        for i, (qname, batch, matches) in enumerate(cobs_iterator(cobs_fn)):
            print(f"Processing batch {batch} query #{i} ({qname})", file=sys.stderr)
            try:
                _ = self._query_dict[qname]
            except KeyError:
                # fallback: unknown query, no sequence
                self._query_dict[qname] = SingleQuery(qname=qname, seq="", keep_matches=self._keep_matches)
            self._query_dict[qname].add_matches(batch, matches)

    def print_tsv_summary(self):
        d = self._query_dict
        for q in d:
            for mtch in d[q]._matches:
                print(q, *mtch, sep="\t")

    def print_fa(self):
        d = self._query_dict
        for q in d:
            frm = d[q].fasta_record_matches()
            print(frm)


def _collect_cobs_records(cobs_fn):
    """Worker function for parallel processing.

    Returns:
        list of (qname, batch, matches) for this cobs file.
    """
    records = list(cobs_iterator(cobs_fn))
    return records


def process_files(query_fn, match_fns, keep_matches):
    # Build the query dictionary once
    sift = Sift(keep_matches=keep_matches, query_fn=query_fn)

    # Parallel parse all cobs files and collect their records
    with ProcessPoolExecutor() as executor:
        future_to_fn = {executor.submit(_collect_cobs_records, fn): fn for fn in match_fns}

        for future in as_completed(future_to_fn):
            fn = future_to_fn[future]
            try:
                records = future.result()
            except Exception as exc:
                print(f"Error while processing {fn}: {exc}", file=sys.stderr)
                continue

            print(f"Finished parsing {fn}", file=sys.stderr)
            # Merge all records into the single Sift instance
            for i, (qname, batch, matches) in enumerate(records):
                try:
                    q_obj = sift._query_dict[qname]
                except KeyError:
                    # fallback: unknown query, create with empty sequence
                    sift._query_dict[qname] = SingleQuery(qname=qname, seq="", keep_matches=keep_matches)
                    q_obj = sift._query_dict[qname]
                q_obj.add_matches(batch, matches)

    # Output final FASTA
    sift.print_fa()


def main():

    parser = argparse.ArgumentParser(description="Collect and threshold COBS matches (parallelized).")

    parser.add_argument(
        'match_fn',
        metavar='',
        nargs='+',
        help='COBS match file(s)',
    )

    parser.add_argument(
        '-q',
        metavar='str',
        dest='query_fn',
        required=True,
        help='query file (FASTA/FASTQ)',
    )

    parser.add_argument(
        '-n',
        metavar='int',
        dest='keep',
        type=int,
        default=DEFAULT_KEEP,
        help=f'no. of best hits to keep [{DEFAULT_KEEP}; 0 = keep all]',
    )

    args = parser.parse_args()
    process_files(args.query_fn, args.match_fn, args.keep)


if __name__ == "__main__":
    main()
