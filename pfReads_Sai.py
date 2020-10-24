# Python 3.5
'''
# Author: Sai Ma
# Date: 10/20/2018
# Objective: Generate a species mixing plot
'''

# Import necessary libraries
import pysam
import argparse

def speciesMixing(bamfile, outfile, p5, BS):
    bam = pysam.AlignmentFile(bamfile, "rb")
    out = open(outfile, 'w')
    barcodes = []
    m = int(p5)
    BarcodeSpace = int(BS) + 1
#    print("Barcoding space is ", BS)
    for x in range(1,BarcodeSpace):
        for y in range(1,BarcodeSpace):
            for z in range(1,BarcodeSpace):
                barcodes.append("R1." + "%02d" % (x,) + ",R2." + "%02d" % (y,) + ",R3." +"%02d" % (z,) + ",P1." +"%02d" % (m,))

    counts = dict()

    for i in barcodes:
        counts[i] = 0
    
    for read in bam.fetch():
        # Read data from bam alignment
        barcode = read.get_tag('RG')
        # Update counts
        if (barcode in barcodes):
            counts[barcode] = counts[barcode] + 1

        # Write counts to file
    out.write("R1,R2,R3,P5,fragments\n")
    for barcode, counts in counts.items():
        output = barcode + "," + str(counts) + "\n"
        out.write(output)

    bam.close()
    out.close()


def main():
    # Collect command line data
    parser = argparse.ArgumentParser(
        description="Generates a csv with barcodes & read counts",
        epilog="Intended for split-ATAC data")
    parser.add_argument(
        "--bam",
        metavar="Path to bam",
        required=True,
        help="Path to processed split-ATAC bamfile")
    parser.add_argument(
        "--out",
        metavar="Path to output csv")
    parser.add_argument(
        "--p5",
        metavar="P5 primer numner, e.g. 31")
    parser.add_argument(
        "--BS",
        metavar="Barcoding space, e.g. 24 or 96")
    args = parser.parse_args()
    # Compile data
    pysam.index(args.bam)
    speciesMixing(args.bam, args.out, args.p5, args.BS)


main()
