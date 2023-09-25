#!/usr/bin/env python

import argparse
import gzip
import pythonvcf
import sys
import os.path

def main():
    parser = argparse.ArgumentParser(description="Filter a VCF file based on ClinVar gold stars")
    parser.add_argument('-v', '--vcf', action='store', type=str, help="The vcf to be filtered", required=True)
    parser.add_argument('-c', '--clinvar', action='store', type=str, help="The vcf from ClinVar", required=True)
    parser.add_argument('-s', '--stars', action='store', type=int, help="The minimum number of stars to keep", required=False, default=2)
    args = parser.parse_args()

    vcf = args.vcf
    clinvar_vcf = args.clinvar
    stars = args.stars

    clinvar_db = {}

    vcf_exists = os.path.isfile(vcf)
    if vcf_exists == False:
        sys.stderr.write("File {} do not exists\n".format(vcf))
        sys.exit(2)

    vcf_exists = os.path.isfile(clinvar_vcf)
    if vcf_exists == False:
        sys.stderr.write("File {} do not exists\n".format(vcf))
        sys.exit(2)

    with (gzip.open if clinvar_vcf.endswith(".gz") else open)(clinvar_vcf) as vcf_content:
        for line in vcf_content:
            if type(line) is str:
                pass
            else:
                line = line.decode('UTF-8')
            if line[0] != '#':
                variant = pythonvcf.Variant_from_clinvar(line)
                number_of_stars = variant.get_clinvar_stars()
                clinvar_db[','.join([variant.chromosome, str(variant.position), variant.reference, variant.alternative])] = number_of_stars

    with (gzip.open if vcf.endswith(".gz") else open)(vcf) as vcf_content:
        for line in vcf_content:
            if type(line) is str:
                pass
            else:
                line = line.decode('UTF-8')
            if line[0] != '#':
                variant = pythonvcf.Variant(line)
                if variant.chromosome.startswith("chr") or variant.chromosome.startswith("Chr"):
                    variant.chromosome = variant.chromosome[3:len(variant.chromosome)]
                key = ','.join([variant.chromosome, str(variant.position), variant.reference, variant.alternative])
                clinvar_variant_stars = clinvar_db.get(key)
                if clinvar_variant_stars == None:
                    pass
                else:
                    if clinvar_variant_stars >= stars:
                        splitted_line = line[:-1].split("\t")
                        info_field = splitted_line[7] + ";GOLD_STARS=" + str(clinvar_variant_stars)
                        splitted_line[7] = info_field
                        print("\t".join(splitted_line))
 
if __name__ == "__main__":
    main()
