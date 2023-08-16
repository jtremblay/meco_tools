#!/usr/bin/env python
 
"""Software to build and access metagenomic data in a hdf5 database.

Julien Tremblay - 
"""
import os
import sys
import argparse
import re
import h5py
import numpy as np 

cwd = os.path.dirname(__file__)
#print(cwd + '/nrc')
sys.path.insert(0, cwd + '/nrc')
#print(sys.path)
from nrc.shotgunmg import *
  
def main(arguments):
    
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    
    subparsers = parser.add_subparsers(title='subcommands', description='valid subcommands', help='additional help', dest="command")
    
    parser_build = subparsers.add_parser('build')
    parser_build.add_argument('-i', '--infile', help="Input file", type=argparse.FileType('r'))
    parser_build.add_argument('-t', '--type', help="Type of database to build (default: shotgunmg)", choices=["shotgunmg", "amplicons"], default="shotgunmg")
    parser_build.add_argument('-d', '--database-file', help="Database file", type=argparse.FileType('w'))
    parser_build.add_argument('-a', '--annotations-file', help="Annotations file from ShotgunMG pipeline output", type=argparse.FileType('r'))
    parser_build.add_argument('-g', '--gene-abundance-matrix-file', help="Gene abundance matrix file from ShotgunMG pipeline output", type=argparse.FileType('r'))
    parser_build.add_argument('-m', '--mapping-file', help="Mapping file used as input in the ShotgunMG pipeline output", type=argparse.FileType('r'))
    #parser_build.add_argument('-j', '--index-annotations-file', help="Link pickle file for linking gene annotation to row", type=argparse.FileType('w'))
    #parser_build.add_argument('-k', '--index-gene-abundance-file', help="Link pickle file for linking gene if to row", type=argparse.FileType('w'))
    
    parser_query = subparsers.add_parser('query')
    #parser_query.add_argument('-i', '--infile', help="Input file", type=argparse.FileType('r'))
    parser_query.add_argument('-t', '--type', help='Type of database to build (default: shotgunmg)', choices=["shotgunmg", "amplicons"], default="shotgunmg")
    parser_query.add_argument('-d', '--database-file', help="Database file", type=argparse.FileType('r'))
    parser_query.add_argument('-o', '--outfile', help="Output file, for <query>", type=argparse.FileType('r'))
    #parser_query.add_argument('-j', '--index-annotations-file', help="Link pickle file for linking gene annotation to row", type=argparse.FileType('r'))
    #parser_query.add_argument('-k', '--index-gene-abundance-file', help="Link pickle file for linking gene if to row", type=argparse.FileType('r'))

    
    args = parser.parse_args(arguments)
    #exit(0);   
    if args.database_file is None:
        raise ValueError('--database-file needed')
    
    if args.command == 'build':
        print("build")
        shotgun_mg = ShotgunMG(args.annotations_file.name, 
                               args.gene_abundance_matrix_file.name, 
                               args.mapping_file.name,
                               args.database_file.name,
                               write = 1,
                               dsetname = "Mock_community") 

        
        shotgun_mg.populate_annotations()
        shotgun_mg.populate_abundance()
        shotgun_mg.populate_mapping()

    elif args.command == 'query':
        print("query")
        shotgun_mg = ShotgunMG(h5db_file = args.database_file.name,
                               write = 0)

        #shotgun_mg.query_gene_id("gene_id_1140")
        #shotgun_mg.dump_mapping()
        #shotgun_mg.query_genus("Neisseria")
        #shotgun_mg.get_sample_ids_from_mapping("Treatment", "Staggered")

        # For instance. Pull our all genes that are 2x fold change between
        # Staggered vs Even for Pseudomonas genus. Display results as average of both conditions
        # or display counts per column. maybe do ascii plots as well.
        shotgun_mg.query(2, ["Treatment", "Staggered", "Even"], ["genus", "Pseudomonas"]) 
    

     
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


