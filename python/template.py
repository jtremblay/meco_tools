#!/usr/bin/env python
 
"""A simple python script template.

Julien Tremblay - jtremblay514@gmail.com
"""
 
import os
import sys
import argparse
import re
  
def main(arguments):
 
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--infile', help="Input file", type=argparse.FileType('r'))
    #iparser.add_argument('-o', '--outfile', help="Output file", default=sys.stdout, type=argparse.FileType('w'))
     
    args = parser.parse_args(arguments)
    infile = os.path.abspath(args.infile.name)
    print infile
     
    print args

    fhand = open(infile)
    count = 0
    for line in fhand:
        count = count + 1
        row = line.split("\t")
        print 'Line Count:', count


     
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


