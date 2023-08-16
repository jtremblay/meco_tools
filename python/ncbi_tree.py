#!/usr/bin/env python
 
"""Tree parsing test scripts.

Julien Tremblay - jtremblay514@gmail.com
"""
 
import os
import sys
import argparse
import math
import PyQt4
from ete2 import PhyloTree
from ete2 import Tree, TreeStyle

def main(arguments):
 
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--infile', help="Input file", type=argparse.FileType('r'))
    parser.add_argument('-o', '--outfile', help="Output file", default=sys.stdout, type=argparse.FileType('w'))
     
    args = parser.parse_args(arguments)
    infile = args.infile
    infile = os.path.abspath(infile.name)
     
    print args

    #t = Tree(infile)
    #print t
    
    #t = Tree( '((H:1,I:1):0.5, A:1, (B:1,(C:1,D:1):0.5):0.5);')
    #print t
    #D = t.search_nodes(name="D")[0]
    #D = t&"D"
    #print D

    #node = D
    #path = []
    #while node.up:
    #    path.append(node.name)
    #    node = node.up

    #print path

    t = Tree( "((a,b),c);" )
    circular_style = TreeStyle()
    circular_style.mode = "c" # draw tree in circular mode
    circular_style.scale = 20
    t.render("mytree.pdf", w=183, units="mm", tree_style=circular_style)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


