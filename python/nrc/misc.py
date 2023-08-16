#!/usr/bin/env python

__author__ = "Julien Tremblay"
__copyright__ = "Copyright 2022, NRC tools" 
__credits__ = ["Julien Tremblay"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Julien Tremblay"
__email__ = ""
__status__ = "Release"


"""Contains code for the shotgunMG uilities project.
"""
import sys
import os
from tqdm import *
import h5py
import numpy as np 
#from numpy import array, concatenate, repeat, zeros, nan, asarray, float, where, isnan, float64, around
#from re import compile, sub
#from skbio.stats.ordination import pcoa
#from skbio.stats.distance import DistanceMatrix

class ShotgunMG:

    def __init__(self, annotations_file="", h5db_file="", number_of_genes=0, number_of_samples, dsetname=""):
        self._h5db_file = h5db_file
        self._annotations_file = annotations_file
        self._dsetname = dsetname
        
        with open(annotations_file) as f:
            # Skip header line
            first_line = f.readline()
            for line in f:
                #line_list = line.split('\t')
                number_of_genes = number_of_genes + 1
            f.close()
        self._number_of_genes = number_of_genes
        sys.stderr.write("Number of genes: " + str(self._number_of_genes) + "\n")
        
        with open(mapping_file) as f:
            # Skip header line
            first_line = f.readline()
            for line in f:
                #line_list = line.split('\t')
                number_of_samples = number_of_samples + 1
            f.close()
        self._number_of_samples = number_of_samples
        sys.stderr.write("Number of samples: " + str(self._number_of_samples) + "\n")

    @property
    def number_of_genes(self):
        return self._number_of_genes
    
    @property
    def annotations_file(self):
        return self._annotations_file

    @annotations_file.setter
    def annotations_file(self, value):
        self._annotations_file = value


    def populate_annotations(self):
        """Parses annotations.tsv output file from ShotgunMG pipeline and stores it into 
           a hdf5 database.

           Returns nothing at the moment.
        """
        def vprint(obj):
            sys.stderr.write("[DEBUG] " + str(obj) + "\n")
        
        #Internal methods for populate_annotations
        """Function to write annotations data into the hdf5 database.
        """
        def load_dict_to_attr(h5f, thisdict, grp_name) :

            if 'gene_id' not in thisdict:
                print('Dictionary missing gene_id key. Skipping function.')
                return
            
            for key, val in thisdict.items():
                h5f[grp_name].attrs[key] = val


        """Function to write annotations data into the hdf5 database - chunk mode.
        """
        def load_dict_to_attr_in_chunks(h5f, thisdict, grp_name) :

            if 'gene_id' not in thisdict:
                print('Dictionary missing gene_id key. Skipping function.')
                return
            
            for key, val in thisdict.items():
                h5f[grp_name].attrs[key] = val

            #return "Blabla..."

        """Show attributes of the hdf5 database.
        """
        def get_grp_attrs(name, node) :

            grp_dict = {}
            for k in node.attrs.keys():
                grp_dict[k]= node.attrs[k]

            print (grp_dict)
       
        h5db_file = self._h5db_file
        annotations_file = self._annotations_file

        # Initialize hdf5 file.
        h5w = h5py.File(h5db_file, 'w')
                    #grp.create_dataset("EmptyDataset", data=h5py.Empty("f"))
        #h5w.create_dataset(self._dsetname, data=h5py.Empty('S'), compression="gzip")

        dataset_name = "annotations"
        #if grp_name not in h5w:
        #    sys.stderr.write('Group: ' + grp_name + ' does not exist and will be created.\n')
        #    grp = h5w.create_group(grp_name)
        #    grp.create_dataset("annotations", dtype=h5py.special_dtype(vlen=str), compression="gzip")
        
        i = 1
        j = 1
        list_of_genes = []        #with tqdm(total=os.path.getsize(annotations_file)) as pbar:
        with tqdm(total=self._number_of_genes) as pbar:
            with open(annotations_file) as f:
                # Skip header line
                # lets try to populate the h5 file with 1000 rows (genes) at a time.
                first_line = f.readline()
                # create the empty dataset
                # Create the dataset at first
                for line in f:
                    line_list = line.split('\t')
                    curr_gene = dict(
                        contig_id            = line_list[0],
                        gene_id              = line_list[1],
                        product_name         = line_list[2],
                        kegg_entry           = line_list[3],
                        kegg_definition      = line_list[4],
                        kegg_module          = line_list[5],
                        kegg_module_desc     = line_list[6],
                        kegg_pathway         = line_list[7],
                        kegg_pathway_desc    = line_list[8],
                        pfam_access          = line_list[9],
                        pfam_product_name    = line_list[10],
                        pfam_desc            = line_list[11],
                        tigrfam_access       = line_list[12],
                        tigrfam_name         = line_list[13],
                        tigrfam_product_name = line_list[14],
                        tigrfam_desc         = line_list[15],
                        cog_access           = line_list[16],
                        cog_name             = line_list[17],
                        cog_function         = line_list[18],
                        cog_category         = line_list[19],
                        kog_access           = line_list[20],
                        kog_name             = line_list[21],
                        kog_function         = line_list[22],
                        kog_category         = line_list[23],
                        ublastp              = line_list[24],
                        ublastp_desc         = line_list[25],
                        tax_kingdom          = line_list[26],
                        tax_phylum           = line_list[27],
                        tax_class            = line_list[28],
                        tax_order            = line_list[29],
                        tax_family           = line_list[30],
                        tax_genus            = line_list[31],
                        tax_specie           = line_list[32]
                    )
                    pbar.update(int(i))
                    list_of_genes.append(curr_gene);

                    if(i % 1000 == 0):
                        vprint("Current gene number: " + str(i))
                        if(j == 1):
                            h5w.create_dataset("annotations", data=list_of_genes, maxshape=(None, 33))
                            vprint("Created dataset...")
                        else:
                            #load_dict_to_attr_in_chunks(h5w, curr_gene, grp_name)
                            h5w['annotations'].resize((h5w['annotations'].shape[0] + 1000), axis=0)
                            h5w['annotations'][-list_of_genes.shape[0]:] = list_of_genes


                        j = j + 1
                        list_of_genes = []
                    i = i + 1
            f.close()


