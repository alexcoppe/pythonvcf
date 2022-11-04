#!/usr/bin/env python

import argparse
import gzip
import os
import sys
import pythonvcf

def main():
    parser = argparse.ArgumentParser(description="Parse a the VCF file from Clinvar")
    parser.add_argument('-v', '--vcf', action='store', type=str, help="The vcf to be parsed", required=True)
    args = parser.parse_args()

    vcf = args.vcf
    with (gzip.open if vcf.endswith(".gz") else open)(vcf) as vcf_content:
        for line in vcf_content:
            if type(line) is str:
                pass
            else:
                line = line.decode('UTF-8')
            if line[0] != '#':
                variant = pythonvcf.Variant_from_clinvar(line)
                print(variant)

if __name__ == "__main__":
    main()
