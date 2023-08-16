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
import csv
import pickle
from tqdm import *
import h5py
import tables as tb
import numpy as np 

def vprint(obj):
    sys.stderr.write("[DEBUG] " + str(obj) + "\n")

class ShotgunMG:

    def __init__(self, annotations_file="", abundance_file="", mapping_file="", h5db_file="", index_annotations_file="", index_gene_abundance_file="", write=0, number_of_genes=0, number_of_samples=0, dsetname=""):
        self._h5db_file = h5db_file
        self._annotations_file = annotations_file
        self._abundance_file = abundance_file
        self._dsetname = dsetname
        self._mapping_file = mapping_file
        self._index_annotations_file = index_annotations_file
        self._index_gene_abundance_file = index_gene_abundance_file
        self._index_annotations = None
        self._index_gene_abundance = None

        if(write == 1):
            #self._h5f_handle = h5py.File(h5db_file, 'w')
            self._h5f_handle = tb.open_file(h5db_file, mode="w", title="")
        elif(write == 0):
            #self._h5f_handle = h5py.File(h5db_file, 'r')
            self._h5f_handle = tb.open_file(h5db_file, mode="r")
        else:
            raise Exception("write argument has to be = 1 or = 0")

        if(annotations_file != ""):
            with open(annotations_file) as f:
                first_line = f.readline()
                for line in f:
                    number_of_genes = number_of_genes + 1
                f.close()
            self._number_of_genes = number_of_genes
            sys.stderr.write("Number of genes: " + str(self._number_of_genes) + "\n")

        if(mapping_file != ""):
            with open(mapping_file) as f:
                first_line = f.readline()
                for line in f:
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
    
    @property
    def abundance_file(self):
        return self._abundance_file

    def populate_annotations(self):
        """Parses annotations.tsv output file from ShotgunMG pipeline and stores it into 
           a hdf5 database.

           Returns nothing at the moment.
        """

        annotations_file = self._annotations_file
        h5w = self._h5f_handle

        dataset_name = "annotations"
        #dt = h5py.special_dtype(vlen=str)
        
        class DataDescr(tb.IsDescription):
            gene_id   = tb.StringCol(32, pos=0)
            contig_id = tb.StringCol(32, pos=1)
            KO        = tb.StringCol(32, pos=2)
            COG       = tb.StringCol(32, pos=3)
            PFAM      = tb.StringCol(32, pos=4)
            kingdom   = tb.StringCol(32, pos=5)
            phylum    = tb.StringCol(32, pos=6)
            Class     = tb.StringCol(32, pos=7)
            order     = tb.StringCol(32, pos=8)
            family    = tb.StringCol(32, pos=9)
            genus     = tb.StringCol(32, pos=10)
            species   = tb.StringCol(32, pos=11)

        with open(annotations_file) as f:
            # Skip header line
            first_line = f.readline()
            i = 1
            #group = h5w.create_group("/", 'gene_abundance', 'Blablabla')
            table = h5w.create_table("/", 'annotations', DataDescr, "Annotation table")
            my_pointer = table.row

            with tqdm(total=self._number_of_genes) as pbar:
                for line in f:
                    line = line.rstrip()
                    line_list = line.split('\t')
                    my_pointer['gene_id']   = line_list[1]
                    my_pointer['contig_id'] = line_list[0]
                    my_pointer['KO']        = line_list[3]
                    my_pointer['COG']       = line_list[9]
                    my_pointer['PFAM']      = line_list[16]
                    my_pointer['kingdom']   = line_list[26]
                    my_pointer['phylum']    = line_list[27]
                    my_pointer['Class']     = line_list[28]
                    my_pointer['order']     = line_list[29]
                    my_pointer['family']    = line_list[30]
                    my_pointer['genus']     = line_list[31]
                    my_pointer['species']   = line_list[32]
                        
                    my_pointer.append()
                    pbar.update(1)
                table.flush() 
            f.close()

    
    def populate_abundance(self):
        """Essentially stores a gene abundance matrix in the hdf5 file. PyTables version 
           Returns nothing at the moment.
        """
        
        abundance_file = self._abundance_file
        h5w = self._h5f_handle

        dataset_name = "gene_abundance"
        number_of_genes = self._number_of_genes
    
        # Create IsDescription structure only for gene id
        class DataDescr(tb.IsDescription):
            gene_id = tb.StringCol(32, pos=0)

        with open(abundance_file) as f:
            # Skip header line
            first_line = f.readline()
            first_line = first_line.rstrip()
            first_line_list = first_line.split('\t')
            first_line_list.pop(0)
            print(first_line_list)
            i = 1
            # Here we dynamically populate the abundance tables.
            for el in first_line_list:
                DataDescr.columns[el] = tb.Int32Col(pos=i)
            
            #h5w = tb.open_file("tutorial1.h5", mode="w", title="Test file")
            #group = h5w.create_group("/", 'gene_abundance', 'Blablabla')
            table = h5w.create_table("/", 'gene_abundance', DataDescr, "Gene abundance matrix")
            my_pointer = table.row

            with tqdm(total=self._number_of_genes) as pbar:
                for line in f:
                    line = line.rstrip()
                    line_list = line.split('\t')
                    row_id = line_list.pop(0)
                    my_pointer['gene_id'] = row_id
                        
                    for j in range(len(line_list)):
                        #vprint(str(j) + " " + row_id + " " + first_line_list[j] + " " + line_list[j])
                        my_pointer[first_line_list[j]] = int(line_list[j])
                    my_pointer.append()
                    pbar.update(1)
                table.flush()
        
        #print(rows)

    def populate_mapping(self):
        """Parses mapping_file that was used in input for the ShotgunMG pipeline and stores it into 
           a hdf5 database. ***Sample ID has to be the first column.***

           Returns nothing at the moment.
        """
        
        mapping_file = self._mapping_file
        h5w = self._h5f_handle
        dataset_name = "mapping"
        
        class DataDescr(tb.IsDescription):
            sample_id = tb.StringCol(256, pos=0)

        with open(mapping_file) as f:
            # Skip header line
            first_line = f.readline()
            first_line = first_line.rstrip()
            first_line_list = first_line.split('\t')
            first_line_list.pop(0)
            print(first_line_list)
            i = 1
            # Here we dynamically populate the abundance tables.
            for el in first_line_list:
                DataDescr.columns[el] = tb.StringCol(64, pos=i)
            
            #h5w = tb.open_file("tutorial1.h5", mode="w", title="Test file")
            #group = h5w.create_group("/", 'gene_abundance', 'Blablabla')
            table = h5w.create_table("/", 'mapping', DataDescr, "Metadata")
            my_pointer = table.row

            with tqdm(total=self._number_of_samples) as pbar:
                for line in f:
                    line = line.rstrip()
                    line_list = line.split('\t')
                    row_id = line_list.pop(0)
                    my_pointer['sample_id'] = row_id
                        
                    for j in range(len(line_list)):
                        my_pointer[first_line_list[j]] = str(line_list[j])
                    my_pointer.append()
                    pbar.update(1)
                table.flush()
        
    def query_gene_id(self, gene_id):
        """Query a gene_id and get its annotations and abundance values.
        """
        h5r = self._h5f_handle
        rows = h5r.root.gene_abundance.read_where('gene_id == b"' + gene_id + '"')
        print(rows)
    
    def query_genus(self, genus):
        """Query a genus and get its annotations and abundance values.
        """

        h5r = self._h5f_handle
        rows = h5r.root.annotations.read_where('genus == b"' + genus + '"')
        print(rows)
    
    def query(self, fc, variables, taxonomy):
        """General query.
        """

        variable = variables[0]
        treatment = variables[1]
        control = variable[2]
        tax_category = taxonomy[0]
        tax_value = taxonomy[1]

        h5r = self._h5f_handle
        #rows = h5r.root.annotations.read_where('genus == b"' + genus + '"')
        my_annotations_table = h5r.root.annotations

        condition = '(' + tax_category + ' == b"' + tax_value + '")'
        vprint(condition)
        #condition = '(name == b"Particle:      5") | (name == b"Particle:      7")'
        for record in my_annotations_table.where(condition):
            #print(record['gene_id'])
            print(type(record))


    def dump_mapping(self):
        """Dump mapping file.
        """
        h5r = self._h5f_handle
        table = h5r.root.mapping
        my_colnames = table.colnames
        my_dict = {}
        for i in range(len(my_colnames)):
            my_dict[str(i)] = my_colnames[i]

        print(my_dict)
        i = 0
        for r in table.iterrows():
            curr_header_name = my_dict[str(i)]
            print(r[curr_header_name], end='   ')
            
            i = i + 1
            if(i == len(my_colnames)):
               print("")
               i = 0

    def get_sample_ids_from_mapping(self, variable, value):
        """Dump mapping file.
        """
        h5r = self._h5f_handle
        table = h5r.root.mapping
        my_colnames = table.colnames
        my_dict = {}
        rows = h5r.root.mapping.read_where(variable + ' == b"' + value + '"')
        print("---")
        print(rows)



