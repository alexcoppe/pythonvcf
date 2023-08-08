#!/usr/bin/env python

import argparse
import gzip
import os
import sys


def build_bed_from_dictionary(line):
    splitted_line = line.split("\t")
    exonStarts = splitted_line[9].split(",")
    exonEnds = splitted_line[10].split(",")
    if len(exonStarts) !=  len(exonEnds):
        return -1
    chrom = splitted_line[2]

    l = len(exonStarts)

    for i in range(l - 1):
        print("{}\t{}\t{}".format(chrom, exonStarts[i], exonEnds[i]))


def main():
    parser = argparse.ArgumentParser(description="Get the bed file from the ncbiRefSeqCurated table downloaded from UCSC Genome Browser")
    parser.add_argument('-f', '--file', action='store', type=str, help="The path to the ncbiRefSeqCurated tsv file", required=True)
    args = parser.parse_args()

    f = args.file

    f_exists = os.path.isfile(f)
    if f_exists == False:
        sys.stderr.write("File {} do not exists\n".format(f_exists))
        sys.exit(2)
    
    d = {}
    dictionary_of_lines = {}

    with (gzip.open if f.endswith(".gz") else open)(f) as file_content:
        for line in file_content:
            if line[0] == "#":
                continue
            line = line[:-1]
            splitted_line = line.split("\t")
            gene = splitted_line[12]
            exoncounts = int(splitted_line[8])
            if gene in d:
                if exoncounts > d[gene]:
                    d[gene] = exoncounts
                    dictionary_of_lines[gene] = line
            else:
                d[gene] = exoncounts
                dictionary_of_lines[gene] = line
    
        for k in d.keys():
            build_bed_from_dictionary(dictionary_of_lines[k])
 
if __name__ == "__main__":
     main()
