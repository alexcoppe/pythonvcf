#!/usr/bin/env python

import argparse
import gzip
import os
import sys
import pythonvcf

def main():
    parser = argparse.ArgumentParser(description="Parse a the VCF file from Clinvar")
    parser.add_argument('-v', '--vcf', action='store', type=str, help="The vcf to be parsed", required=True)
    #parser.add_argument('-n', '--name', action='store', type=str, help="The sample name", required=False, default=None)
    #parser.add_argument('-f', '--first', action='store_true', help="Print only the first variant")
    #parser.add_argument('-t', '--no_header', action='store_true', help="Do not print the header")
    #parser.add_argument('-s', '--samples', action='store', type=int, help="The number of samples (1 by default)", default=1)
    args = parser.parse_args()

    vcf = args.vcf
    #name = args.name
    #first = args.first
    #no_header = args.no_header
    #number_of_samples = args.samples

    #samples_columns = create_samples_columns_names(number_of_samples)

    #vcf_exists = os.path.isfile(vcf)
    #if vcf_exists == False:
        #sys.stderr.write("File {} do not exists\n".format(vcf))
        #sys.exit(2)


    #trait = '.'

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
