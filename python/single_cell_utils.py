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
import scipy.io
import scrublet as scr
import numpy as np
from itertools import compress
import csv

cwd = os.getcwd()
sys.path.append(cwd + '/nrc')
from nrc.fastq_iterator import *
  
def main(arguments):
    
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(title='subcommands', description='valid subcommands', help='additional help', dest="command")
    parser_dmc1 = subparsers.add_parser('scrublet')
    parser_dmc1.add_argument('-i', '--infile-mtx', help=".mtx file (output of STAR)", type=argparse.FileType('r'))
    parser_dmc1.add_argument('-g', '--infile-genes', help="genes file (output of STAR)", type=argparse.FileType('r'))
    parser_dmc1.add_argument('-b', '--infile-barcodes', help="barcodes file (output of STAR)", type=argparse.FileType('r'))
    parser_dmc1.add_argument('-o', '--outdir', help="output directory", type=Path)
    parser_dmc1.add_argument('-v', '--min-gene-variability', type=int, default=78, help='min_gene_variability_pctl')
    parser_dmc1.add_argument('-m', '--min-counts', type=int, default=2, help='min_counts')
    parser_dmc1.add_argument('-c', '--min-cells', type=int, default=3, help='min_cells')
    parser_dmc1.add_argument('-p', '--n-prin-comps', type=int, default=30, help='n_prin_comps')
    parser_dmc1.add_argument('-e', '--expected-doublet-rate', type=float, default=0.08, help='Expected doublet rate')
    parser_dmc1.add_argument('--dry_run', default=False, action=argparse.BooleanOptionalAction, help='Dry run.')

    args = parser.parse_args(arguments)
    
    if args.command == 'scrublet':
        print("Running scrublet workflow", file=sys.stderr)
    
        infile_mtx = args.infile_mtx.name
        infile_genes = args.infile_genes.name
        outdir = args.outdir.name
        infile_barcodes = args.infile_barcodes.name

        print(infile_mtx, file=sys.stderr)
        print(infile_genes, file=sys.stderr)
        print(infile_barcodes, file=sys.stderr)
        print(outdir, file=sys.stderr)
        os.makedirs(outdir, exist_ok=True)

        base = os.path.splitext(os.path.basename(infile_mtx))[0]
        print("Basename : " + base, file=sys.stderr)

        counts_matrix = scipy.io.mmread(infile_mtx).T.tocsc()
        features = np.array(scr.load_genes(infile_genes, delimiter='\t', column=1))
        with open(args.infile_barcodes.name, 'r') as file:
            barcodes = [line.rstrip() for line in file]

        print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]), file=sys.stderr)
        print('Number of genes/features in gene list: {}'.format(len(features)), file=sys.stderr)
        print('Number of barcodes(i.e. cells) in barcodes list: {}'.format(len(barcodes)), file=sys.stderr)

        scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=args.expected_doublet_rate)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=args.min_counts, min_cells=args.min_cells, min_gene_variability_pctl=args.min_gene_variability, n_prin_comps=args.n_prin_comps)
        scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.08)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=78, n_prin_comps=30)
        predicted_singlets = [not elem for elem in predicted_doublets]

        """Then write files minus doublets"""
        counts_matrix_orig = scipy.io.mmread(infile_mtx).tocsc()
        mtx_doublets = counts_matrix_orig[:,predicted_doublets]
        barcodes_doublets = list(compress(barcodes, predicted_doublets))
        mtx_singlets = counts_matrix_orig[:,predicted_singlets]
        barcodes_singlets = list(compress(barcodes, predicted_singlets))
       
        print("Number of doublets: " + str(predicted_doublets.tolist().count(True)), file=sys.stderr)
        print("Number of singlets: " + str(predicted_doublets.tolist().count(False)), file=sys.stderr)
        scrub.plot_histogram()[0].savefig(os.path.join(outdir, base + "_histogram.png"))

        print('Running UMAP...', file=sys.stderr)
        scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
        scrub.plot_embedding('UMAP', order_points=True)[0].savefig(os.path.join(outdir, base + "_UMAP.png"))
        
        scipy.io.mmwrite("matrix_doublets.mtx", mtx_doublets)
        with open('barcodes_doublets.txt', 'w') as f:
            for barcode in barcodes_doublets:
                f.write(f"{barcode}\n")
       
        scipy.io.mmwrite("matrix_singlets.mtx", mtx_singlets)
        with open('barcodes_singlets.txt', 'w') as f:
            for barcode in barcodes_singlets:
                f.write(f"{barcode}\n")

        
     
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


