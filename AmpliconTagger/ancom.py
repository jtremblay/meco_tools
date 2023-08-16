#!/usr/bin/env python3
 
"""Implementation of ancomP (i.e. Python3 implementation of ANCOM).
returns OTUs differentially abundant between two conditions. Treatment vs control

Julien Tremblay - jtremblay514@gmail.com
"""

import sys
import os
import argparse
from ancomP.stats.ancom import ancom
import pandas as pd
import pathlib
import numpy as np
import csv
import multiprocessing
import logging
os.system("taskset -p 0xf %d" % os.getpid())
sys.stderr.write(str(os.system("taskset -p 0xf %d" % os.getpid())))

def main(arguments):
    pd.options.mode.chained_assignment = None
    
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--otu-table', help="OTU table file in tsv format - do not accept .biom files.", type=argparse.FileType('r'), required=True)
    parser.add_argument('-m', '--mapping-file', help="Mapping file in tsv format.", type=argparse.FileType('r'), required=True)
    parser.add_argument('-o', '--outdir', help="Outdir file containing differentially abundant OTUs, their logFC, CPM and taxonomy. One file per pairwise comparison.", type=str, required=True)
    parser.add_argument('-v', '--variable', help="Column in the mapping file in which all the possible pairwise comparisons will be performed.", type=str, default=1.5, required=True)
    parser.add_argument('-t', '--num-threads', help="Number of threads", type=int)
    
    args = parser.parse_args(arguments)
    otu_table_file = os.path.abspath(args.otu_table.name)
    mapping_file = os.path.abspath(args.mapping_file.name)
    outdir = os.path.abspath(args.outdir)
    variable = str(args.variable)
    num_threads = args.num_threads


    sys.stderr.write("[DEBUG]\n")
    sys.stderr.write("OTU table file: " + otu_table_file + "\n")
    sys.stderr.write("Mapping file: " + mapping_file + "\n")
    sys.stderr.write("outdir: " + outdir + "\n")
    sys.stderr.write("variable (i.e. column name): " + variable + "\n")
    sys.stderr.write("num_threads: " + str(num_threads) + "\n")
   
    # log
    fh_log = open(os.path.join(outdir, "ancom_log.txt"), "w")

    # for test:
    #otu_table_file = '~/Projects/Lallemand_AAD/16S/export/otu_tables/otu_table_filtered.tsv'
    #mapping_file = '~/Projects/Lallemand_AAD/16S/export/mapping_file.tsv'
    #outdir = '/home/jtrembla/Projects/Lallemand_AAD/16S/ancomp'
    #variable = 'Treatment'
    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)

    # Read OTU table
    otu_table = pd.read_csv(otu_table_file, sep='\t', skiprows=0, header=0, index_col=0)
    
    # Read mapping file and get all values for specified variable.
    mapping = pd.read_csv(mapping_file, sep='\t', skiprows=0, header=0, index_col=0)
    variables = mapping.loc[:, variable].unique().tolist()
    sys.stderr.write("Variables in selected variable: " + str(variables))


    # Build a list of task and execute them later.
    tasks = []
    for i in range(0,(len(variables)-1)):
        control = variables[i]
        sys.stderr.write("Current control : " + control + "\n")
        
        for j in range(i+1,len(variables)):
            treatment = variables[j]
            sys.stderr.write("Current treatment : " + treatment + "\n")
            # Select by row.
            mapping2 = mapping.loc[mapping[variable].isin([control,treatment])]

            samples_control = mapping2.loc[mapping[variable].isin([control])].index.tolist()
            samples_treatment = mapping2.loc[mapping[variable].isin([treatment])].index.tolist()
            
            samples_control = mapping2.loc[mapping[variable].isin([control])].index
            samples_treatment = mapping2.loc[mapping[variable].isin([treatment])].index
            # and double check that samples are actually in otu_table i.e. some samples can be left out after rarefaction.
            samples_control = samples_control[samples_control.isin(list(otu_table))].tolist()
            samples_treatment = samples_treatment[samples_treatment.isin(list(otu_table))].tolist()

             # Select by col.
            otu_table2 = otu_table.loc[:, samples_control + samples_treatment]
            otu_table2 = otu_table2 + 1
            column_names = list(otu_table2)
            design_control = [0] * len(samples_control)
            design_treatment = [1] * len(samples_treatment)
    

            tasks.append( (design_control, design_treatment,
                           samples_control, samples_treatment,
                           control, treatment,
                           otu_table2, otu_table, outdir, ) )
            
            #results = p.apply_async(do_work, [design_control, design_treatment,
            #                        samples_control, samples_treatment,
            #                        control, treatment,
            #                        otu_table2, otu_table, outdir])


    #p.close()
    #p.join()
    
    p = multiprocessing.Pool(num_threads)
    #os.system('taskset -cp 0-%d %s' % (num_threads, os.getpid()))

    for t in tasks:
        p.apply_async(do_work, t)
    p.close()
    p.join()

def do_work(design_control, design_treatment, samples_control, samples_treatment, control, treatment, otu_table2, otu_table, outdir):
    #sys.stderr.write("design_control: " + str(design_control) + "\n")
    #sys.stderr.write("design_treatment: " + str(design_treatment) + "\n")
    #sys.stderr.write("samples_control: " + str(samples_control) + "\n")
    #sys.stderr.write("samples_treatment: " + str(samples_treatment) + "\n")
    #sys.stderr.write("outdir: " + str(outdir) + "\n")
    print("outfile: " + os.path.join(outdir, treatment + '_vs_' + control + '.tsv') + "\n")

    grouping = pd.Series(design_control + design_treatment, 
                         index=samples_control + samples_treatment
    )
    for col in (samples_control + samples_treatment):
        otu_table2[col] = otu_table2[col].astype('float',copy=False)
        
    otu_table_transposed = otu_table2.T
    results = ancom(otu_table_transposed, grouping, significance_test='permutative-anova', permutations=100)
    
    otus = results.loc[results['reject'] == True].index.tolist()
    otu_table3 = otu_table2.loc[otu_table2.index.isin(otus)]
    otu_table3['Average_control'] = otu_table3.loc[:,samples_control].mean(axis=1)
    otu_table3['Average_treatment'] = otu_table3.loc[:,samples_treatment].mean(axis=1)
    otu_table3['Stdev_control'] = otu_table3.loc[:,samples_control].std(axis=1)
    otu_table3['Stdev_treatment'] = otu_table3.loc[:,samples_treatment].std(axis=1)
    otu_table3 = otu_table3.loc[:, ['Average_control','Average_treatment','Stdev_control','Stdev_treatment']]
    otu_table3['FC'] = otu_table3['Average_treatment'] / otu_table3['Average_control']
    otu_table3['logFC'] = np.log2(otu_table3['FC'])
    otu_table3 = otu_table3.apply(lambda x: pd.Series.round(x, 4))
    otu_table3 = pd.merge(otu_table3, otu_table.loc[:,['taxonomy']], how='left', left_index=True, right_index=True)
    pd.DataFrame(otu_table3).to_csv(os.path.join(outdir, treatment + '_vs_' + control + '.tsv'), sep='\t', quoting=csv.QUOTE_NONE)
    
    # Then do reverse: control becomes treatment and treatment becomes control
    otu_table3 = otu_table2.loc[otu_table2.index.isin(otus)]
    otu_table3['Average_control'] = otu_table3.loc[:,samples_treatment].mean(axis=1)
    otu_table3['Average_treatment'] = otu_table3.loc[:,samples_control].mean(axis=1)
    otu_table3['Stdev_control'] = otu_table3.loc[:,samples_treatment].std(axis=1)
    otu_table3['Stdev_treatment'] = otu_table3.loc[:,samples_control].std(axis=1)
    otu_table3 = otu_table3.loc[:, ['Average_control','Average_treatment','Stdev_control','Stdev_treatment']]
    otu_table3['FC'] = otu_table3['Average_treatment'] / otu_table3['Average_control']
    otu_table3['logFC'] = np.log2(otu_table3['FC'])
    otu_table3 = otu_table3.apply(lambda x: pd.Series.round(x, 4))
    otu_table3 = pd.merge(otu_table3, otu_table.loc[:,['taxonomy']], how='left', left_index=True, right_index=True)
    pd.DataFrame(otu_table3).to_csv(os.path.join(outdir, control + '_vs_' + treatment + '.tsv'), sep='\t', quoting=csv.QUOTE_NONE)

    #fh_log.write("Completed ancom " + str(os.path.join(outdir, control + '_vs_' + treatment + '.tsv')) + "\n")
    print("Completed ancom " + str(os.path.join(outdir, control + "_vs_" + treatment + ".tsv")) + "\n")
    sys.stdout.flush()


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

