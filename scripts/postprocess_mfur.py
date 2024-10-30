#! /usr/bin/env python3

import argparse
import collections
import os
import re
import sys


# def get_nb_kmers(mfur_line):
#     p = mfur_line.split("\t")
#     nb_kmers = int(p[-1])
#     return nb_kmers


def keep_fname(mfur_fname):
    filename_with_extension = os.path.basename(mfur_fname)  # Remove the path
    filename_without_extension, _, _, = filename_with_extension.partition(".") # Remove extension (only 1 dot is expected)
    return "_" + filename_without_extension


def process_mfur_output(hits_to_keep): #format the output as 
    prev_query = -1
    for x in sys.stdin:
        query, filename, score = x.split("\t")
        
        # 1. Check if matches were found
        if filename == "NA": #no matches: skip to the next iteration (i.e. next query)
            print(x, end="")

        else:
           # 2. Check if the current match is related to the same query of the previous line 
            if prev_query != query: #re-assign if first line or new query, and restart counter + min score
                prev_query = query
                i = 0
                min_kmers = 0
            elif prev_query == query: #increase counter otherwise 
                i += 1

        fname_only = keep_fname(mfur_fname) #clean matched filename

        # 3. Extract the Top N results (+ties)
        if i < hits_to_keep:
            print(query + "\t" + fname_only + "\t" + score, end="")
        elif i == hits_to_keep:
            print(query + "\t" + fname_only + "\t" + score, end="")
            min_kmers = int(score)
        else:
            if int(score) == min_kmers:
                print(query + "\t" + fname_only + "\t" + score, end="")


def main():

    parser = argparse.ArgumentParser(
        description="Postprocess cobs output: keep top n hits (+ties) and remove random identifiers")

    parser.add_argument(
        '-n',
        metavar='int',
        dest='keep',
        required=True,
        type=int,
        help=f'no. of best hits to keep',
    )

    args = parser.parse_args()

    process_mfur_output(args.keep)


if __name__ == "__main__":
    main()
