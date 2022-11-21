#!/usr/bin/env python

import argparse
import gzip
import os
import sys
import pythonvcf

def main():
    parser = argparse.ArgumentParser(description="Parse a the VCF file from Clinvar")
    parser.add_argument('-v', '--vcf', action='store', type=str, help="The vcf to be parsed", required=True)
    parser.add_argument('-i', '--index', action='store_true', help="Add an index column")
    args = parser.parse_args()

    vcf = args.vcf
    index = args.index
    i = 1

    with (gzip.open if vcf.endswith(".gz") else open)(vcf) as vcf_content:
        for line in vcf_content:
            if type(line) is str:
                pass
            else:
                line = line.decode('UTF-8')
            if line[0] != '#':
                if "\"" in line:
                    line = line.replace("\"", "")
                variant = pythonvcf.Variant_from_clinvar(line)
                if index == True:
                    print("{}\t{}\t{}".format(i, variant, variant.get_clinvar_stars()))
                else:
                    print("{}\t{}".format(variant, variant.get_clinvar_stars()))
                i += 1

if __name__ == "__main__":
    main()
