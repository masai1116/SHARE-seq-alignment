#!/usr/bin/env python

# Author: Jason Buenrostro, Stanford University
# This will extract insert size distribution from input peaks

# The following program will extract insert information from peaks

##### IMPORT MODULES #####
# import necessary for python
import os
import sys
import numpy as np
from optparse import OptionParser
from bx import intervals
import scipy.io as sio
#from bx.intervals.intersection import Intersecter, Interval
# http://bcbio.wordpress.com/2009/04/29/finding-and-displaying-short-reads-clustered-in-the-genome/
# http://nullege.com/codes/search/bx.intervals.Intersecter?fulldoc=1

##### DEFINE FUNCTIONS #####

#### OPTIONS ####
# read options from command line
opts = OptionParser()
usage = "usage: %prog [options] [inputs]"
opts = OptionParser(usage=usage)
opts.add_option("-i", help="<Peaks to score> ATAC calls")
opts.add_option("-o", help="<Output>")
opts.add_option("--top", type="int", default=50000, help="<Number> number of summits to select")  # overrides above setting
options, arguments = opts.parse_args()

# return usage information if no argvs given
if len(sys.argv)==1:
    os.system(sys.argv[0]+" --help")
    sys.exit()

##### INPUTS AND OUTPUTS #####
# get intervals
print "Loading data..."
p1_ints = np.loadtxt(options.i,'str',delimiter="\t")
print "Calculating rank..."

# save as fast lookup
bedInts = {}
for j in range(0,len(p1_ints)):
    chr= p1_ints[j][0]; start = int(p1_ints[j][1]); end= int(p1_ints[j][2])
    if not bedInts.has_key(chr):
        bedInts[chr] = intervals.Intersecter()
    bedInts[chr].insert(start,end,np.append(p1_ints[j],j))

# pass filter idx
idx = np.array(np.ones(len(p1_ints)),dtype=int)

##### SCRIPT #####
# loop through files
print "Looping through intervals..."
for i in range(0,len(p1_ints)):
    # look for overlaps
    atacChr = p1_ints[i][0]
    atacStart = int(p1_ints[i][1])
    atacEnd = int(p1_ints[i][2])
    
    # if pass filter and only 1 match
    if idx[i] == 1:
        findData = np.array(bedInts[atacChr].find(atacStart,atacEnd))
        pos = np.array(np.array(findData)[:,-1],dtype='int')
        if len(findData[np.where(idx[pos])]) == 1:
            continue
    else:
        continue
    
    # if more than one
    pos = np.array(np.array(findData)[:,-1],dtype='int')
    while len(findData[np.where(idx[pos])]) > 1:
        a = np.array(findData[np.where(idx[pos])])
        a = a[np.array(a[:,-1],dtype=float).argsort()]
        idx[int(a[-1][-1])] = 0

# define output name
print "Printing top summits..."
np.savetxt(options.o,p1_ints[np.where(idx)][0:options.top],fmt='%s', delimiter='\t')

# sign out
print "Printed to: "+options.o
print "Created by Jason Buenrostro."
print "Completed."
