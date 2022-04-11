#!/usr/bin/env python

import argparse
import gzip
import os
import sys
import pythonvcf


def build_set_from_bed(bed_file_name):
    file_exists = os.path.isfile(bed_file_name)
    if file_exists == False:
        sys.stderr.write("File {} do not exists\n".format(bed_file_name))
        sys.exit(2)

    genes_list = []
    with (gzip.open if bed_file_name.endswith(".gz") else open)(bed_file_name) as file_content:
        for line in file_content:
            try:
                gene_name = line[:-1].split()[3]
            except:
                sys.stderr.write("""Problem getting the gene name in this line:
                        {}""".format(line))
                sys.exit(2)
            genes_list.append(gene_name)
    genes_set = set(genes_list)
    return genes_set
            

def main():
    parser = argparse.ArgumentParser(description="Subset a VCF based on a list of genes")
    parser.add_argument('-v', '--vcf', action='store', type=str, help="The VCF file path to be filtered", required=True)
    parser.add_argument('-f', '--filtered_genes', action='store', type=str, help="The bed file path containing the genes to be filtered", required=True)
    args = parser.parse_args()

    vcf = args.vcf
    file_with_wanted_genes = args.filtered_genes

    set_of_wanted_genes_mutations = build_set_from_bed(file_with_wanted_genes)

    vcf_exists = os.path.isfile(vcf)
    if vcf_exists == False:
        sys.stderr.write("File {} do not exists\n".format(vcf))
        sys.exit(2)

    with (gzip.open if vcf.endswith(".gz") else open)(vcf) as vcf_content:
        for line in vcf_content:
            if type(line) is str:
                pass
            else:
                line = line.decode('UTF-8')
            if line[0] != '#':
                variant = pythonvcf.Variant(line)
                for transcript in variant.snpeff_transcipt_list:
                    gene = transcript.gene
                    if gene in set_of_wanted_genes_mutations:
                        if transcript.impact == "MODERATE" or transcript.impact == "HIGH":
                            print(line[:-1])
            else:
                print(line[:-1])


if __name__ == "__main__":
    main()
