#!/usr/bin/env python
 
"""A simple python script template.

Julien Tremblay - jtremblay514@gmail.com
"""
 
import os
import sys
import argparse
import re
import csv
import itertools
from collections import defaultdict
from pprint import pprint 
 
def main(arguments):
 
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--infile', help="Input file", type=argparse.FileType('r'))
    parser.add_argument('-f', '--fields', help="fields to put in mapping file.", type=str, required=True)
     
    args = parser.parse_args(arguments)
    infile = os.path.abspath(args.infile.name)
    fields = args.fields

    #sys.stderr.write("infile: " + infile + "\n")
    #sys.stderr.write("fields: " + fields + "\n")
    #sys.stderr.write(args)

    # Split field into dict
    field_list = fields.split(",")
    #print field_list


    fh = open(infile)
    reader=csv.DictReader(fh, delimiter="\t")
    
    # Declare dict to store values.
    #hash = defaultdict(list)
    
    count = 0
    sys.stdout.write("#SampleID");
    k = 1;
    for field in field_list:
        field = field.replace("(", "")
        field = field.replace(")", "")
        field = field.replace(" ", "")
        field = field.replace("_", "")
        field = field.replace("-", "")
        field = field.replace("/", "")
        field = field.replace("\\", "")
        #if(k == 1):
        #    sys.stdout.write(field)
        #else:
        sys.stdout.write("\t" + field)

        #k = k+1

    sys.stdout.write("\n")
    for row in reader:
        sample_id = row['Sample_ID']
        #hash['Sample_ID'] = sample_id
        #sys.stdout.write(hash['Sample_ID'] + "\t")
        sample_id = sample_id.replace(" ", ".")
        sample_id = sample_id.replace("_", ".")
        sample_id = sample_id.replace("-", ".")
        sample_id = sample_id.replace("/", ".")
        sample_id = sample_id.replace("\\" ,".")
        sys.stdout.write(sample_id + "\t")

        for field in field_list:
            #hash[field] = row[field]
            element = row[field]
            element = element.replace(" ",".")
            element = element.replace("_",".")
            element = element.replace("-",".")
            element = element.replace("/",".")
            element = element.replace("\\",".")
            sys.stdout.write(element + "\t")
 
        sys.stdout.write("\n")
     
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


