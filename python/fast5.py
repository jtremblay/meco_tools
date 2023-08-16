#!/usr/bin/env python
 
"""A simple script to parse fast5 reads.

Julien Tremblay - 
"""
 
import os
import sys
import argparse
import re
import h5py
  
def main(arguments):
 
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--infile', help="Input file", type=argparse.FileType('r'))
    #iparser.add_argument('-o', '--outfile', help="Output file", default=sys.stdout, type=argparse.FileType('w'))
     
    args = parser.parse_args(arguments)
    infile = os.path.abspath(args.infile.name)
    print infile
     
    print args

    #fhand = open(infile)
    #count = 0
    #for line in fhand:
    #    count = count + 1
    #    row = line.split("\t")
    #    print 'Line Count:', count

    filename = 'my_file.fast5'
    datasetname = '/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq'
     
     with h5py.File(filename, 'r') as hdf:
             fq = hdf[datasetname][()]
                 print fq


     
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


