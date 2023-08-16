#!/usr/bin/env python

__author__ = "Julien Tremblay"
__copyright__ = "Copyright 2023, Julien Tremblay" 
__credits__ = ["Julien Tremblay"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Julien Tremblay"
__email__ = "jtremblay514@gmail.com"
__status__ = "Release"


"""Contains general utility code in support of the scrna/dna pipeline.
"""
import sys
import os
from re import compile, sub
from itertools import islice

class Fastq(object):

    N = 4

    def __init__(self, fastq_file=""):
        self._fastq_file = fastq_file
        #self._fh = None
        self._fh = open(self._fastq_file)
        self._curr_entry = {}
    
    #def create_fh(self):
    def __iter__(self):
        return self

    def next(self):
        next_n_lines = list(islice(self._fh, Fastq.N))
        if not next_n_lines:
            self._fh.close()
            return False
        else:
            self._curr_entry = next_n_lines
            return True
    
    @property
    def seq(self):
        return (self._curr_entry[1].rstrip())
    
    #@property
    def seq_by_coords(self, x, y):
        return (self._curr_entry[1][int(x):int(y)])
   
    @property
    def qual(self):
        return self._curr_entry[3].rstrip()
    
    @property
    def header(self):
        return self._curr_entry[0].rstrip()
    
    @property
    def output(self):
        return self._curr_entry[0].rstrip() + "\n" + self._curr_entry[1].rstrip() + "\n" + self._curr_entry[2].rstrip() + "\n" + self._curr_entry[3].rstrip() + "\n"

