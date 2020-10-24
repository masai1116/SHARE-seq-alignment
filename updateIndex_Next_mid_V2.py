# Python 3.6
'''
# Author: Sai Ma
# Date: 04/02/2018
# Objective: Splitting split ATAC projects as fastqs
'''
import argparse
import itertools
import gzip 
import os

def extractFastq(r1, r2, out, i1, i2):
    # print("Opening files")
    R1 = gzip.open(r1, 'r')
    R2 = gzip.open(r2, 'r')
    Index1 = gzip.open(i1, 'r')
    Index2 = gzip.open(i2, 'r')
    outR1 = gzip.open((out + ".R1.fastq.gz"), 'w')
    outR2 = gzip.open((out + ".R2.fastq.gz"), 'w')
    print(("Add index to read name of " + os.path.basename(r1)))
    i = 0
    totalReads = 0
    valid = False
    for line in zip(R1, R2, Index1, Index2):
        if (i == 0):
            Temp1 = line[0][0:-2]
            Temp2 = line[1][0:-2]
            totalReads += 1
        if (i == 1):
            id1 = line[2][0:-1]
            id2 = line[3]
            Temp3 = (Temp1 + id1 + "+" + id2)
            Temp4 = (Temp2 + id1 + "+" + id2)
        if (i == 1):
            outR1.write(Temp3)
            outR2.write(Temp4)
            outR1.write(line[0])
            outR2.write(line[1])
        if (i == 2):
            outR1.write(line[0])
            outR2.write(line[1])
        if (i == 3):
            outR1.write(line[0])
            outR2.write(line[1])
        i += 1
        if (i == 4):
            i = 0
    R1.close()
    R2.close()
    outR1.close()
    outR2.close()


def main():
    parser = argparse.ArgumentParser(
        description="Splits fastqs based on splitATAC project",
        epilog="Intended for splitATAC data")
    parser.add_argument(
        "-R1",
        metavar="Read 1",
        required=True,
        help="Path to the Read 1 fastq file")
    parser.add_argument(
        "-R2",
        metavar="Read 2",
        required=True,
        help="Path to the Read 2 fastq file")
    parser.add_argument(
        "--out",
        metavar="Output",
        required=True,
        help="Path to the output fastq files")
    parser.add_argument(
        "-Index1",
        metavar="Output",
        required=True,
        help="Path to the index1 files")
    parser.add_argument(
        "-Index2",
        metavar="Output",
        required=True,
        help="Path to the index2 files")
    args = parser.parse_args()
    extractFastq(args.R1, args.R2, args.out, args.Index1, args.Index2)


main()
