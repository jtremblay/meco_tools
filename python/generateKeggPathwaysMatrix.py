#!/usr/bin/env python
 
"""Takes input in blastp against kegg database. (parsed output that is...)
and an abundance (gene count) table and output a KEGG-module abundance table.
Julien Tremblay - jtremblay514@gmail.com
"""
 
import os
import sys
import argparse
import re
import csv
import itertools
from collections import defaultdict
import operator
import numpy as np
import pprint

def main(arguments):

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-b', '--infile-blastp', help="Input file", type=argparse.FileType('r'), required=True)
    
    args = parser.parse_args(arguments)
    infile_blastp = os.path.abspath(args.infile_blastp.name)
    sys.stderr.write(infile_blastp + "\n")
    sys.stderr.write(str(args) + "\n")
   
    ######
    # Parse blastp 
    # First establish a hash of all modules in blastp output.
    fhand=open(infile_blastp, "r")
    reader=csv.DictReader(fhand, delimiter="\t")
    hash_modules = defaultdict(list)

    for row in reader:
        module_id_field    = row["PATHWAY_ID"]
        module_ids = module_id_field.split("==")
        for module_id in module_ids:
            hash_modules[module_id] = module_id
  
    hash_modules.pop('NULL', None)
    sys.stdout.write("gene_id")
    for key in sorted(hash_modules.iterkeys()):
        sys.stdout.write("\t" + key)
    sys.stdout.write("\n")
   
    # Repoen fhand
    fhand=open(infile_blastp, "r")
    reader=csv.DictReader(fhand, delimiter="\t")
    final_hash = defaultdict(list)
    for row in reader:
        gene_id            = row["#query"]
        module_id_field    = row["PATHWAY_ID"]

        # split MODULE_ID AND MODULE_DESC
        hash_ko = defaultdict(list)
        module_ids = module_id_field.split("==")
        for module_id in module_ids:
            hash_ko[gene_id].append(module_id)

        for gene_id in hash_ko.iterkeys():
            # get gene id
            module_ids = map(str, hash_ko[gene_id])

            sys.stdout.write(str(gene_id))
            for module in sorted(hash_modules.iterkeys()):
                if module in module_ids:
                    sys.stdout.write("\t1")
                else:
                    sys.stdout.write("\t0")
            sys.stdout.write("\n")

    sys.stderr.write("Done parsing blastp results\n")
    #pprint.pprint(final_hash)
    #pprint.pprint(hash_ko)
    #pprint.pprint(hash_modules)
    exit(0)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

