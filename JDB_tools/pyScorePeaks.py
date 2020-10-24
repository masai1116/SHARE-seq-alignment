#!/usr/bin/env python

# Author: Jason Buenrostro, Stanford University
# This will extract insert size distribution from input peaks

# The following program will extract insert information from peaks

##### IMPORT MODULES #####
# import necessary for python
import os
import sys
import pysam
import numpy as np
import matplotlib.pyplot as plt
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
opts.add_option("-a", help="<Peaks to score> ATAC calls")
opts.add_option("-b", help="<Directory> Directory of calls")
opts.add_option("--append", action="store_true", default=False, help="Append output to bed file")
opts.add_option("--PWM", action="store_true", default=False, help="Append output to bed file")
opts.add_option("--noBED", action="store_true", default=False, help="Append output to bed file")
opts.add_option("-o", help="<Output>")
options, arguments = opts.parse_args()

# return usage information if no argvs given
if len(sys.argv)==1:
    os.system(sys.argv[0]+" --help")
    sys.exit()

##### INPUTS AND OUTPUTS #####
# get intervals                                                                                                                                                                                                                                                                            
p1_ints = np.loadtxt(options.a,'str',delimiter="\t")
options.b = options.b+'/'

# remove unwanted chromosomes
#idx = []
#for i in range(0,len(p1_ints[:,0])):
#    if p1_ints[i,0].split('chr')[1].isdigit():
#        idx.append(i)
#    elif p1_ints[i,0].split('chr')[1] == 'X':
#        idx.append(i)
#p1_ints = p1_ints[idx,:]

# get list of dir contents
os.system("rm temp1 temp2 ReadFastqs")
cmd='ls -vdp  '+options.b+'* | grep -v README | grep -v "/$" > ReadChIP'
#cmd='ls -v  '+options.b+'* | grep -v README > ReadChIP'
os.system(cmd)
files = np.loadtxt('ReadChIP','str')

##### SCRIPT #####
# initialize out
try:iter = len(files)
except:files = [files]
outData = np.zeros([len(p1_ints),len(files)])
headers = []

# loop through files
for i in range(0,len(files)):
    # open file
    file = files[i].tolist()
    print 'Processing: ',file
    fName = file.split('/')[-1].split('.bed')[0].split('.narrowPeak')[0]
    #chipData = np.loadtxt(options.b+'/'+files[i],'str')
    chipData = np.loadtxt(file,'str') 
    headers.append([fName, len(chipData)])
        
    # save as fast lookup
    chip = {}
    for j in chipData:
        chr= j[0]; start = int(j[1]); end= int(j[2])
        # save a dict of chroms
        if not chip.has_key(chr):
            chip[chr] = intervals.Intersecter()
        # add data
        chip[chr].insert(start,end,j)
    
    # loop through peaks
    for j in range(0,len(p1_ints)):
        # get info
        atacChr = p1_ints[j][0]
        atacStart = int(p1_ints[j][1])
        atacEnd = int(p1_ints[j][2])
        
        # find peak
        try: findData = chip[atacChr].find(atacStart,atacEnd)
        except KeyError: out = []; findData=''
        if findData:
            if options.PWM: outData[j,i]=np.round(np.max(np.array(np.array(findData)[:,4],dtype='float')),3)
            else: outData[j,i]+=1
        else:
            pass


# define output name
print "Extracting signal complete, preparing to print..."
if not options.o:
    out = os.path.basename(options.a)+os.path.basename(options.b)+".mat"
else:
    out = options.o

# save data
if options.append == True:
    oldData = sio.loadmat(out)
    oldHeaders = []
    for i in oldData['dataHeaders']:
        oldHeaders.append([i[0][0], i[1][0][0]])
    # append to exisiting
    outData = np.hstack((oldData['peakAssoc'],outData))
    headers = np.vstack((oldHeaders,headers))

# save
if options.noBED:
    sio.savemat(out,{'peakAssoc':outData,'dataHeaders':np.array(headers, dtype=np.object)},oned_as='column',do_compression='True')
else:
    sio.savemat(out,{'peakAssoc':outData,'bed':np.array(p1_ints, dtype=np.object),'dataHeaders':np.array(headers, dtype=np.object)},oned_as='column',do_compression='True')

# sign out
print "Created by Jason Buenrostro."
print "Completed."
