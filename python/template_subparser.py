#!/usr/bin/env python
 
"""A simple python script template.

Julien Tremblay - 
"""
 
import os
import sys
import argparse
import re
  
def main(arguments):
    
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(title='subcommands', description='valid subcommands', help='additional help', dest="command")
    parser_bd = subparsers.add_parser('betadiv')
    parser_bd.add_argument('-i', '--infile', help="Input file", type=argparse.FileType('r'))
    parser_bd.add_argument("-m", "--metric", help="Diversity metric (default: bray-curtis)", choices=["bray-curtis", "weighted-unifrac"], default="bray-curtis")
    #parser_bd.set_defaults(func=betadiv)
    parser_ts = subparsers.add_parser('taxsum')
    parser_ts.add_argument('-i', '--infile', help="Input file", type=argparse.FileType('r'))
    parser_ts.add_argument("-t", "--type", help="Summary type (default: absolute)", choices=["asolute", "relative"], default="absolute")
    #parser_bd.set_defaults(func=taxsum)


    parser.add_argument('-i', '--infile', help="Input file", type=argparse.FileType('r'))
    #parser.add_argument('-m', '--infile-matrix', help="Input file", type=argparse.FileType('r'))
    #iparser.add_argument('-o', '--outfile', help="Output file", default=sys.stdout, type=argparse.FileType('w'))
    

    args = parser.parse_args(arguments)
    
    if args.command == 'betadiv':
        print("betadiv")
    elif args.command == 'taxsum':
        print("taxsum")
    
    #infile = os.path.abspath(args.infile.name)
    #print(infile) #STDOUT
    #sys.stderr.write(infile + "\n") #STDERR
     
    print(args)

    #fhand = open(infile)
    #count = 0
    #for line in fhand:
    #    count = count + 1
    #    print('Line Count:', count)
    #    print(line)
    #    row = line.split("\t")
    #    print(row)
    #
    #print('Loop completed')

     
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


