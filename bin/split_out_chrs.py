#!/usr/bin/env python3

import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""Split out haploid fastqs into separate files"""
    )
    parser.add_argument("-i", type=str, help="reference_fasta")
    parser.add_argument("--chrs", type=str, nargs='+', help="contigs to split off")

    args = parser.parse_args()

    chrs = ['others'] + args.chrs
    outfiles = {c:open(c+'.fa', 'w') for c in chrs}

    outfile = outfiles['others']
    with open(args.i, 'r') as ref:
        for line in ref:
            if '>' == line[0]:
                chrom = line[1:].split(' ')[0]
                if chrom in chrs:
                    outfile = outfiles[chrom]
            outfile.write(line)

    for outfile in outfiles.values():
        outfile.close()
