#!/usr/bin/env bash

set -x

(
    for x in mfur/*.mfur; do
        printf "%s\t%d\t%d\n" "$x" $(cat $x | wc -c) $(cat $x | wc -c);
    done
) > data/decompressed_indexes_sizes.txt
