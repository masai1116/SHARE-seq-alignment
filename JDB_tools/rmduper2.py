#!/usr/bin/env python
import pysam
import argparse
import os

def rmduper2(bamfile, outfile):
    # Basic remove duplicate function
    bam = pysam.AlignmentFile(bamfile,"rb")

    keep = [0 for read in bam.fetch() if read.is_proper_pair]

    start = 0
    chrom = ''

    chrom_names = bam.references
    fragment_dict = {}

    forward_fragment_names = {}

    i = 0

    for read in bam.fetch():

        if not read.is_proper_pair:
            continue        

        if not read.is_reverse:
            
            l_pos = read.pos
            ilen = abs(read.template_length)
            r_pos = read.pos + ilen - 1
            chr_pos = chrom_names[read.reference_id]
            
            if l_pos > start or chr_pos != chrom:
            
                for key in fragment_dict.keys():

                    if len(fragment_dict[key]) == 1:
                        keep[fragment_dict[key][0][0]] = 1
                        forward_fragment_names[fragment_dict[key][0][1]] = 1
                    else:
                        for j in range(len(fragment_dict[key])):
                            forward_fragment_names[fragment_dict[key][j][1]] = 0
                
                chrom = chr_pos
                start = l_pos
                fragment_dict = {(l_pos,r_pos): [(i, read.query_name, read.mapping_quality)]}
                
            else:
                try:
                    fragment_dict[(l_pos,r_pos)].append((i, read.query_name, read.mapping_quality))
                except:
                    fragment_dict[(l_pos,r_pos)] = [(i, read.query_name, read.mapping_quality)]                           

        i += 1        
    
    for key in fragment_dict.keys():

        if len(fragment_dict[key]) == 1:
            keep[fragment_dict[key][0][0]] = 1
            forward_fragment_names[fragment_dict[key][0][1]] = 1
        else:
            for j in range(len(fragment_dict[key])):
                forward_fragment_names[fragment_dict[key][j][1]] = 0
    
    i = 0
    for read in bam.fetch():
        if not read.is_proper_pair:
            continue        
        if read.is_reverse:
            keep[i] = forward_fragment_names.pop(read.query_name)
        i += 1

    out = pysam.AlignmentFile(outfile,"wb",template=bam)

    i = 0
    for read in bam.fetch():
        if not read.is_proper_pair:
            continue        
        if keep[i] == 1:
            out.write(read)
        i += 1 

    bam.close()
    out.close()


def main():
    
    argparser = argparse.ArgumentParser(description = "%(prog)s -- Removes all copies of duplicated reads",
                                        epilog = "Designed for single-cell ATAC-seq")
    argparser.add_argument("--bam", metavar = "bamfile", required= True, help = "Bam file to be de-duplicated" )
    argparser.add_argument("--out", metavar = "output_name", help = "output name")
    args = argparser.parse_args() 


    if args.out is None:
        args.out = os.path.basename(args.bam)[:-4] + ".deduped.bam"

    rmduper2(args.bam, args.out)
    pysam.index(args.out)

if __name__ == "__main__":
    main()







