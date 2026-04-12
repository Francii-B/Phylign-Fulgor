#!/usr/bin/env python3
import yaml
import sys

with open("config.yaml", "r") as stream:
    try:
        config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc, file=sys.stderr)
        sys.exit(1)

# check if fulgor_threads is an int
try:
    fulgor_threads = config["fulgor_threads"]
    int(fulgor_threads)
except (TypeError, ValueError):
    print(
        "ERROR: to run Phylign in cluster mode, the parameter fulgor_threads in config.yaml MUST BE SET to a fixed "
        "int value. Aborting.",
        file=sys.stderr)
    sys.exit(1)
