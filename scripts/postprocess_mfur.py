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
    filename_without_extension, _, _, = filename_with_extension.partition(".") # Remove extension (first dot)
    return "_" + filename_without_extension


def process_mfur_output(hits_to_keep):
    """convert mfur output into COBS format

    Tempororary solution to avoid modifications in 'filter_queries.py'
    """
    prev_query = ""
    for x in sys.stdin:
        query, filename, score = x.split("\t")
        
        # 1. Check if the line is related to a different query
        if query != prev_query:
            print("*" + query, end="") #the n. of matches can be ignored (not used later)
            prev_query = query
            i = 0

        # 2. Check if matches were found for the query
        if filename != "NA": 
            i += 1 #increase match-counter
            fname_only = keep_fname(filename) #clean matched filename

            # 3. Extract the Top N results (+ties)
            if i < hits_to_keep:
                print(fname_only + "\t" + score, end="")
            elif i == hits_to_keep:
                print(fname_only + "\t" + score, end="")
                min_kmers = int(score)
            else:
                if int(score) == min_kmers:
                    print(fname_only + "\t" + score, end="")


def main():

    parser = argparse.ArgumentParser(
        description="Postprocess mfur output: keep top n hits (+ties) and remove random identifiers")

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
