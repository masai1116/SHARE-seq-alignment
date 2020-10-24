# Python 3.6
'''
# Author: Sai Ma
# Date: 08/29/2018
# Objective: Splitting split ATAC projects as fastqs
'''
import argparse
import itertools

def barcodeSet(barcode):
    bases = "ATCGN"
    barcodeSet = set()
    barcodeSet.add(barcode)
    for i, c in enumerate(barcode):
        if c in bases:
            for base in bases:
                if c != base:
                    barcodeSet.add((barcode[:i] + base + barcode[i + 1:]))
    return barcodeSet

def extractFastq(r1, r2, out):
    R1 = open(r1, 'r')
    R2 = open(r2, 'r')
    outR1 = open((out + ".R1polyT"), 'w')
    outR2 = open((out + ".R2polyT"), 'w')
#    print(("Filter out mismatch poly dT of " + r1))
    i = 0
    Flag = 0
    barcode = 'TTTTTT'
    correctT = barcodeSet(barcode)
    for line in zip(R1, R2):
        if (i == 0):
            Temp1 = line[0]
            Temp2 = line[1]
        if (i == 1):
            if line[1][0:6] in correctT:
                Flag = 1
#                print(Flag)
                outR1.write(Temp1)
                outR2.write(Temp2)
                outR1.write(line[0])
                outR2.write(line[1])
        if (i == 2):
            if (Flag == 1):
                outR1.write(line[0])
                outR2.write(line[1])
        if (i == 3):
            if (Flag == 1):
                outR1.write(line[0])
                outR2.write(line[1])
        i += 1
        if (i == 4):
            i = 0
            Flag = 0
    R1.close()
    R2.close()
    outR1.close()
    outR2.close()


def main():
    parser = argparse.ArgumentParser(
        description="Match polyT on read2 for splitATAC project",
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
    args = parser.parse_args()

    extractFastq(args.R1, args.R2, args.out)

main()
