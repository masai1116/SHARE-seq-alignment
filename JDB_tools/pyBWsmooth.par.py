#!/usr/bin/env python

# Author: Jason Buenrostro, Stanford University
# This will make a bigWig of ATAC data

##### IMPORT MODULES #####
# import necessary for python
import os
import sys
import numpy as np
import subprocess
from optparse import OptionParser
import pysam
from multiprocessing import Pool
from bx.intervals.io import GenomicIntervalReader
from bx.bbi.bigwig_file import BigWigFile
# http://bcbio.wordpress.com/2009/04/29/finding-and-displaying-short-reads-clustered-in-the-genome
# http://nullege.com/codes/search/bx.intervals.Intersecter?fulldoc=1  

#### OPTIONS ####
# read options from command line
opts = OptionParser()
usage = "usage: %prog [options] [inputs]"
opts = OptionParser(usage=usage)
opts.add_option("-p",default=20,type='int', help="<Threads> Accepts an integer")
opts.add_option("-a", help="<bw> Accepts a bigwig file")
opts.add_option("-g", help="<Genome Size file>")
opts.add_option("-w", default=150,type='int', help="<Int> window size")
opts.add_option("-s", default=20,type='int', help="<Int> step size (span)")
opts.add_option("-c", default=0,type='int', help="chr chunk")
options, arguments = opts.parse_args()

# return usage information if no argvs given
if len(sys.argv)==1:
    os.system(sys.argv[0]+" --help")
    sys.exit()

##### DEFINE FUNCTIONS ##### 
def bwSmooth(c):
    # open bigwig
    t = open(options.a)
    bw = BigWigFile(t)
    # get data, pass if not available 
    chrN = c[0];sPos=int(c[1]);ePos=int(c[2])
    signal = bw.get_as_array(chrN,sPos,ePos)
    t.close()
    # smooth data
    if type(signal) == type(None): signal = np.zeros(ePos-sPos)
    #else:
    signal[np.isnan(signal)] = 0
    convM = np.convolve(signal,wSmooth,'same')
    # save
    sList = np.arange(sPos,ePos,step)
    eList = sList+step
    chrList = np.array([chrN]*len(sList))
    meanSig = convM[range(step/2,chunkSize+padLen,step)]
    # save out
    idx1 = meanSig>0; idx2 = eList<chrLen; idx = idx1*idx2
    idx[chunkSize/step:] = False
    pData = np.c_[chrList[idx],np.array(sList[idx],dtype=str),eList[idx],meanSig[idx]]
    return pData

##### INPUTS AND OUTPUTS #####
# get gSize file
gSizes = np.loadtxt(options.g,'str')
chunkSize = 2500000
padLen = 5000

# open out file
outName = options.a+'.smooth.'+str(options.c)+'.bed'
try: os.remove(outName)
except OSError: pass
outF = file(outName, 'a')

# smoothParems
wSize = options.w
wSmooth = np.ones(wSize)
step = options.s

#### SCRIPT #####
# split genome into chunks
#for i in range(0,len(gSizes)):
#for i in range(20,len(gSizes)):
for i in range(options.c-1,options.c): 
    # break chrs into pieces
    chrN = gSizes[i][0]
    chrLen = int(gSizes[i][1])
    sVals = np.arange(1,int(gSizes[i][1]),chunkSize)
    
    # organize
    chrName = [gSizes[i][0]]*len(sVals)
    chunks = np.c_[chrName,np.array(sVals,dtype=str),np.array(sVals+chunkSize+padLen,dtype=str)]
    
    # run on all pieces
    print "Smoothing: "+chrName[0]
    pool = Pool(processes=20)
    sigArray = pool.map(bwSmooth,chunks)
    # save to file
    print "Printing data for: "+chrName[0]
    np.savetxt(outF,np.vstack(sigArray), fmt='%s',delimiter='\t', newline='\n')

# close outF
outF.close()

# convert to bw
#print 'Converting to bigwig...'
#out = options.a.split('.dgf.bw')[0]+'.s'+str(step)+'.w'+str(wSize)+'sw.bw'
#tool='~wjg/Applications/UCSC_Tools/bedGraphToBigWig'
#os.system(tool+' '+'out.smooth.bed'+' '+options.g+' '+out)

# sign out
print "Created by Jason Buenrostro."
print "Completed"
