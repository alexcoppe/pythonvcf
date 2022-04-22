#!/usr/bin/env python

import argparse
import gzip
import os
import sys
import pythonvcf

def build_set_from_txt(txt_file_name, l=[0]):
    file_exists = os.path.isfile(txt_file_name)
    if file_exists == False:
        sys.stderr.write("File {} do not exists\n".format(txt_file_name))
        sys.exit(2)
    genes_list = []
    with (gzip.open if txt_file_name.endswith(".gz") else open)(txt_file_name) as file_content:
        for line in file_content:
            try:
                splitted_line = line[:-1].split("\t")
            except:
                sys.stderr.write("""Problem getting the gene names in this line:
                        #{}""".format(line))
                sys.exit(2)
            for i in l:
                try:
                    genes_list.append(splitted_line[i])
                except IndexError:
                    sys.stderr.write("File {} has not the {} field\n".format(txt_file_name, i))
                    sys.exit(2)
            #genes_list.append(splitted_line[2])
    genes_set = set(genes_list)
    return genes_set
 



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
    parser.add_argument('-f', '--filtered_genes', action='store', type=str, help="The bed/txt file path containing the genes to be filtered (bed file or txt file with second and thirs columns with gene name)", required=True)
    parser.add_argument('-t', '--type_of_genes_file', action='store', type=str, help="The bed or the txt file with the secondo or third coloumns containing a gene name", required=False, default="bed")
    args = parser.parse_args()

    vcf = args.vcf
    file_with_wanted_genes = args.filtered_genes
    type_of_genes_file = args.type_of_genes_file

    if type_of_genes_file == "bed":
        set_of_wanted_genes_mutations = build_set_from_bed(file_with_wanted_genes)
    else:
        set_of_wanted_genes_mutations = build_set_from_txt(file_with_wanted_genes, [1,2])

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
