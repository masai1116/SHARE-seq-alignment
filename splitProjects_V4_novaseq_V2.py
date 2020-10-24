# Python 3.6
'''
# Author: Sai Ma
# Date: 10/20/2018
# Objective: Splitting split ATAC projects as fastqs
'''
import argparse
import itertools
import yaml
import gzip

P5 = {'TAGATCGC': 'P1.01',
      'CTCTCTAT': 'P1.02',
      'TATCCTCT': 'P1.03',
      'AGAGTAGA': 'P1.04',
      'GTAAGGAG': 'P1.05',
      'ACTGCATA': 'P1.06',
      'AAGGAGTA': 'P1.07',
      'CTAAGCCT': 'P1.08',
      'TGGAAATC': 'P1.09',
      'AACATGAT': 'P1.10',
      'TGATGAAA': 'P1.11',
      'GTCGGACT': 'P1.12',
      'TTTCTAGC': 'P1.13',
      'TAACCAAG': 'P1.14',
      'GTGTATCG': 'P1.15',
      'TCCATCAA': 'P1.16',
      'TTCGTGCA': 'P1.17',
      'AGGTTGCC': 'P1.18',
      'CCTTATGT': 'P1.19',
      'CAGCAACG': 'P1.20',
      'GGTTCAAT': 'P1.21',
      'ACATTCGT': 'P1.22',
      'GATTCCCA': 'P1.23',
      'CGGACTGC': 'P1.24',
      'AGCCGTTC': 'P1.25',
      'ATTGGGTC': 'P1.26',
      'TGCATACT': 'P1.27',
      'GGGCTTGG': 'P1.28',
      'GACGTGGC': 'P1.29',
      'GCAAATTT': 'P1.30',
      'GCAGCCTC': 'P1.31',
      'TCCGAGTT': 'P1.32',
      'GCATTAAG': 'P1.33',
      'ACGATAAC': 'P1.34',
      'CCTGCGGG': 'P1.35',
      'TGATTGTT': 'P1.36',
      'GGCACGGA': 'P1.37',
      'GATCATTC': 'P1.38',
      'ATGGTCAT': 'P1.39',
      'CGTACCAA': 'P1.40',
      'CCAGTTTA': 'P1.41',
      'ACCGGCCC': 'P1.42',
      'CTAGAAGT': 'P1.43',
      'CGCCAGAT': 'P1.44',
      'TCACATGG': 'P1.45',
      'GAACTCGA': 'P1.46',
      'CCACCGTT': 'P1.47',
      'TAAGTTAC': 'P1.48',
      'GAGACGTG': 'P1.49',
      'TTGCCTAA': 'P1.50',
      'TTAACTTG': 'P1.51',
      'CTTTAACA': 'P1.52',
      'CGTAGACC': 'P1.53',
      'TATTTGCG': 'P1.54',
      'ATCCAGGA': 'P1.55',
      'TGTTCCTG': 'P1.56',
      'ACGCGCAG': 'P1.57',
      'TCTGGCGA': 'P1.58',
      'AATCTACA': 'P1.59',
      'TACTGACC': 'P1.60',
      'CGATAGGG': 'P1.61',
      'ACTTAGAA': 'P1.62',
      'AGAGATCT': 'P1.63',
      'GGTGAAGG': 'P1.64',
      'ATCGAATG': 'P1.65',
      'TCAAGAGC': 'P1.66',
      'GCCCACGT': 'P1.67',
      'TGGGCGGT': 'P1.68',
      'CCCTTGGA': 'P1.69',
      'ATTACCGT': 'P1.70',
      'AGTCCGAG': 'P1.71',
      'ACTTGTTG': 'P1.72',
      'GTAATACA': 'P1.73',
      'GGCGTCTA': 'P1.74',
      'GCGCTGCT': 'P1.75',
      'GTGCCATT': 'P1.76',
      'TAGGTATG': 'P1.77',
      'AACACCTA': 'P1.78',
      'CTCCGAAC': 'P1.79',
      'CAACGGCA': 'P1.80',
      'CAATGTAG': 'P1.81',
      'GGCTACCC': 'P1.82',
      'AAAGTCCG': 'P1.83',
      'TTCCGCGG': 'P1.84',
      'AGGCACTT': 'P1.85',
      'CTTCAGTG': 'P1.86',
      'GCCGGTAG': 'P1.87',
      'TTCAATCC': 'P1.88',
      'CCACACAC': 'P1.89',
      'ATATTATC': 'P1.90',
      'CCGAAGCA': 'P1.91',
      'GTATCGGT': 'P1.92'}

Round1 = {'ATCACGTT': 'R1.01',
      'CGATGTTT': 'R1.02',
      'TTAGGCAT': 'R1.03',
      'TGACCACT': 'R1.04',
      'ACAGTGGT': 'R1.05',
      'GCCAATGT': 'R1.06',
      'CAGATCTG': 'R1.07',
      'ACTTGATG': 'R1.08',
      'GATCAGCG': 'R1.09',
      'TAGCTTGT': 'R1.10',
      'GGCTACAG': 'R1.11',
      'CTTGTACT': 'R1.12',
      'TGGTTGTT': 'R1.13',
      'TCTCGGTT': 'R1.14',
      'TAAGCGTT': 'R1.15',
      'TCCGTCTT': 'R1.16',
      'TGTACCTT': 'R1.17',
      'TTCTGTGT': 'R1.18',
      'TCTGCTGT': 'R1.19',
      'TTGGAGGT': 'R1.20',
      'TCGAGCGT': 'R1.21',
      'TGATACGT': 'R1.22',
      'TGCATAGT': 'R1.23',
      'TTGACTCT': 'R1.24',
      'TGCGATCT': 'R1.25',
      'TTCCTGCT': 'R1.26',
      'TAGTGACT': 'R1.27',
      'TACAGGAT': 'R1.28',
      'TCCTCAAT': 'R1.29',
      'TGTGGTTG': 'R1.30',
      'TACTAGTC': 'R1.31',
      'TTCCATTG': 'R1.32',
      'TCGAAGTG': 'R1.33',
      'TAACGCTG': 'R1.34',
      'TTGGTATG': 'R1.35',
      'TGAACTGG': 'R1.36',
      'TACTTCGG': 'R1.37',
      'TCTCACGG': 'R1.38',
      'TCAGGAGG': 'R1.39',
      'TAAGTTCG': 'R1.40',
      'TCCAGTCG': 'R1.41',
      'TGTATGCG': 'R1.42',
      'TCATTGAG': 'R1.43',
      'TGGCTCAG': 'R1.44',
      'TATGCCAG': 'R1.45',
      'TCAGATTC': 'R1.46',
      'TAGTCTTG': 'R1.47',
      'TTCAGCTC': 'R1.48',
      'TGTCTATC': 'R1.49',
      'TATGTGGC': 'R1.50',
      'TTACTCGC': 'R1.51',
      'TCGTTAGC': 'R1.52',
      'TACCGAGC': 'R1.53',
      'TGTTCTCC': 'R1.54',
      'TTCGCACC': 'R1.55',
      'TTGCGTAC': 'R1.56',
      'TCTACGAC': 'R1.57',
      'TGACAGAC': 'R1.58',
      'TAGAACAC': 'R1.59',
      'TCATCCTA': 'R1.60',
      'TGCTGATA': 'R1.61',
      'TAGACGGA': 'R1.62',
      'TGTGAAGA': 'R1.63',
      'TCTCTTCA': 'R1.64',
      'TTGTTCCA': 'R1.65',
      'TGAAGCCA': 'R1.66',
      'TACCACCA': 'R1.67',
      'TGCGTGAA': 'R1.68',
      'GGTGAGTT': 'R1.69',
      'GATCTCTT': 'R1.70',
      'GTGTCCTT': 'R1.71',
      'GACGGATT': 'R1.72',
      'GCAACATT': 'R1.73',
      'GGTCGTGT': 'R1.74',
      'GAATCTGT': 'R1.75',
      'GTACATCT': 'R1.76',
      'GAGGTGCT': 'R1.77',
      'GCATGGCT': 'R1.78',
      'GTTAGCCT': 'R1.79',
      'GTCGCTAT': 'R1.80',
      'GGAATGAT': 'R1.81',
      'GAGCCAAT': 'R1.82',
      'GCTCCTTG': 'R1.83',
      'GTAAGGTG': 'R1.84',
      'GAGGATGG': 'R1.85',
      'GTTGTCGG': 'R1.86',
      'GGATTAGG': 'R1.87',
      'GATAGAGG': 'R1.88',
      'GTGTGTCG': 'R1.89',
      'GCAATCCG': 'R1.90',
      'GACCTTAG': 'R1.91',
      'GCCTGTTC': 'R1.92',
      'GCACTGTC': 'R1.93',
      'GCTAACTC': 'R1.94',
      'GATTCATC': 'R1.95',
      'GTCTTGGC': 'R1.96'}

def barcodeSet(barcode):
    bases = "ATCGN"
    barcodeSet = set()
    barcodeSet.add(barcode)
    for i, c in enumerate(barcode):
        if c in bases:
            for base in bases:
                if c != base:
                    barcodeSet.add((barcode[:i] + base + barcode[i + 1:]))
                    barcodeSet.add((barcode[0:5]))
#                    barcodeSet.add((barcode[1:] + base))
#                    barcodeSet.add((base + barcode[:-1]))
    return barcodeSet
# allow 1 mismatch or 1 bp shift

def extractFastq(r1, r2, conf, project, out):
    R1 = gzip.open(r1, 'r')
    R2 = gzip.open(r2, 'r')
    inFile = open(conf, 'r')
    config = yaml.load(inFile, Loader=yaml.UnsafeLoader)
    projectNames = set()
    for key in config.keys():
        if ('Project' in key):
            projectNames.add(key)
    for proj in projectNames:
        metaData = config[proj]
        if (project == metaData['Name']):
            outR1 = gzip.open((out + ".R1.fq.gz"), 'w')
            outR2 = gzip.open((out + ".R2.fq.gz"), 'w')
            primer = dict()
            r1 = dict()
            primersub = metaData['Primer']
            r1sub = metaData['Round1']
            for barcode, name in P5.items():
                if (name in primersub):
                    barcodes = barcodeSet(barcode)
                    for bc in barcodes:
                        primer[bc] = name
            for barcode, name in Round1.items():
                if (name in r1sub):
                    barcodes = barcodeSet(barcode)
                    for bc in barcodes:
                        r1[bc] = name
            i = 0
            totalReads = 0
            sortedReads = 0
            barcodeMatch = 0
            for line in zip(R1, R2):
                if (i == 0):
                    index = line[0].find("+")
                    p5 = line[0][(index + 1):(index + 9)]
                    if (p5 in primer):
                        barcodeMatch += 1
                    elif (line[0][(index + 1):(index + 6)] in primer):
                        barcodeMatch += 1
                        p5 = line[0][(index + 1):(index + 6)]
                    indexR = line[0].find("1:N:0:")
                    round1 = line[0][(indexR + 21):(indexR + 29)]
                    if (round1 in r1):
                        barcodeMatch += 1
                    if (barcodeMatch == 2):
                        valid = True
                    else:
                        valid = False
                barcodeMatch = 0
                if (valid):
                    if (i == 0):
                        index1 = line[0].find(" ")
                        end1 = line[0].find("+")
                        index2 = line[1].find(" ")
                        end2 = line[1].find("+")
                        outR1.write(line[0][:index1] +
                                    " " +
                                    line[0][index1 + 1:end1] +
                                    "_" + p5 +
                                    "\n")
                        outR2.write(line[1][:index2] +
                                    " " +
                                    line[1][index2 + 1:end2] +
                                    "_" + p5 +
                                    "\n")
                    else:
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
        metavar="read 1",
        required=True,
        help="Path to the Read 1 fastq file")
    parser.add_argument(
        "-R2",
        metavar="Read 2",
        required=True,
        help="Path to the Read 2 fastq file")
    parser.add_argument(
        "-yaml",
        metavar="yaml",
        required=True,
        help="Path to yaml file contains barcodes for each project")
    parser.add_argument(
        "-Project",
        metavar="Project Name",
        required=True,
        help="Name of the project")
    parser.add_argument(
        "-out",
        metavar="Output",
        required=True,
        help="Path to the output fastq files")
    args = parser.parse_args()
    extractFastq(args.R1, args.R2, args.yaml, args.Project, args.out)


main()
