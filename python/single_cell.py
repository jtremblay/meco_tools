#!/usr/bin/env python
 
"""Python utilities for single cell rna/dna sequencing data type.

Julien Tremblay - jtremblay514@gmail.com
"""

import os
import sys
import argparse
import re
import json
from pathlib import Path

cwd = os.getcwd()
sys.path.append(cwd + '/nrc')
from nrc.fastq_iterator import *
  
def main(arguments):
    
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(title='subcommands', description='valid subcommands', help='additional help', dest="command")
    parser_dmc1 = subparsers.add_parser('demultiplex_config_1')
    parser_dmc1.add_argument('-i', '--infile-rna3p', help="Input fastq file containing 3' RNA reads", type=argparse.FileType('r'))
    parser_dmc1.add_argument('-b', '--infile-barcodes', help="Input barcode file containing 3' RNA reads", type=argparse.FileType('r'))
    parser_dmc1.add_argument('-c', '--cell_index_coords',  nargs='+', action="store", default=[0,16], help='Index sequence corresponding to the single cell index. Usually the first 16 bp.')
    parser_dmc1.add_argument('-u', '--UMI_index_coords', nargs="+",  action="store", default=[16,26], help='Index sequence corresponding to the single cell index. Usually the first 10 bp following the cell index.')
    parser_dmc1.add_argument('-o', '--outdir', help="output directory", type=Path)
    parser_dmc1.add_argument('--dry_run', default=False, action=argparse.BooleanOptionalAction, help='Just print the stats without generating demultiplexed files.')

    #parser_dmc1.set_defaults(func=betadiv)
    parser_ts = subparsers.add_parser('blabla')
    parser_ts.add_argument('-i', '--infile', help="Input file", type=argparse.FileType('r'))
    parser_ts.add_argument("-t", "--type", help="Summary type (default: absolute)", choices=["asolute", "relative"], default="absolute")

    args = parser.parse_args(arguments)
    
    if args.command == 'demultiplex_config_1':
        print("demultiplex_config_1")
    elif args.command == 'taxsum':
        print("taxsum")
    
     
    print(args)
    infile_rna = args.infile_rna3p.name
    infile_barcodes = args.infile_barcodes.name
    outdir = args.outdir
    x = args.cell_index_coords[0]
    y = args.cell_index_coords[1]

    fastq_db_rna = Fastq(infile_rna)
    fastq_db_barcodes = Fastq(infile_barcodes)

    i = 0
    fhs = {}
    count = {}
    while(fastq_db_barcodes.next() == True):
        fastq_db_rna.next()
        cell_barcode = fastq_db_barcodes.seq_by_coords(args.cell_index_coords[0], args.cell_index_coords[1])
        if("N" in cell_barcode):
            print(cell_barcode, file=sys.stderr)
            continue

        umi_barcode = fastq_db_barcodes.seq_by_coords(args.UMI_index_coords[0], args.UMI_index_coords[1])
        
        """Open one file handler for each cell"""
        if(args.dry_run == True):
            if(cell_barcode not in count):
                count[cell_barcode] = 1
            else:
                count[cell_barcode] += 1
        else:
            if(cell_barcode not in fhs):
                fhs[cell_barcode] = open(os.path.join(outdir, cell_barcode + ".fastq"), "w")
                print(fastq_db_rna.header + " #" + cell_barcode + "-" + umi_barcode + "\n" + fastq_db_rna.seq + "\n+\n" + fastq_db_rna.qual, file=fhs[cell_barcode])
                count[cell_barcode] = 1
                fhs[cell_barcode].close()
            else:
                fhs[cell_barcode] = open(os.path.join(outdir, cell_barcode + ".fastq"), "a")
                print(fhs[cell_barcode], file=sys.stderr)
                print(fastq_db_rna.header + " #" + cell_barcode + "-" + umi_barcode + "\n" + fastq_db_rna.seq + "\n+\n" + fastq_db_rna.qual, file=fhs[cell_barcode])
                count[cell_barcode] += 1
                fhs[cell_barcode].close()

        i = i + 1
        #if(i > 1000000):
        #print(count, file=sys.stderr)
    print(json.dumps(count, indent=4, sort_keys=False), file=sys.stderr)


     
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


