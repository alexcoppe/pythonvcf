#!/usr/bin/env python

import argparse
import gzip
import pythonvcf
import sys
import os.path


class Genome_range:
    def __init__(self, chromosome, start, end):
        self.chromosome = chromosome
        self.start = start
        self.end = end

    def __str__(self):
        return("{} {} {}".format(self.chromosome, self.start, self.end))

def get_ranges_from_bed(bed_path):
    l = []
    with (gzip.open if bed_path.endswith(".gz") else open)(bed_path) as bed_content:
        for line in bed_content:
            splitted_line = line.strip().split()
            if line.startswith('browser') or line.startswith('track') or line.startswith("#"):
                continue
            number_of_spaces = line.count(' ')
            number_of_tabs = line.count('\t')
            if number_of_tabs < 2 and number_of_spaces < 3:
                return -1
            chromosome = splitted_line[0]
            start = int(splitted_line[1])
            end = int(splitted_line[2])
            genome_range = Genome_range(chromosome, start, end)
            l.append(genome_range)
    return l

def is_variant_in_ranges(variant, bed_ranges):
    variant_chromosome = variant.chromosome
    if variant_chromosome.startswith("chr") or variant_chromosome.startswith("Chr"):
        variant_chromosome = variant_chromosome[3:]
    for bed_range in bed_ranges:
        if bed_range.chromosome ==  variant_chromosome:
            start = bed_range.start
            end = bed_range.end
            if variant.position > start and variant.position < end:
                return 1
    return 0

def main():
    parser = argparse.ArgumentParser(description="Filter a VCF using a BED file")
    parser.add_argument('-v', '--vcf', action='store', type=str, help="The vcf file path", required=True)
    parser.add_argument('-b', '--bed', action='store', type=str, help="The bed file path", required=True)
    args = parser.parse_args()

    vcf = args.vcf
    bed = args.bed

    if os.path.exists(vcf) == False:
        sys.exit("{} .VCF file do not exists".format(bed))

    if os.path.exists(bed) == False:
        sys.exit("{} .BED file do not exists".format(bed))
    
    bed_ranges = get_ranges_from_bed(bed)

    if bed_ranges == -1:
        sys.exit("{} is not a correct .BED file".format(bed))

    with (gzip.open if vcf.endswith(".gz") else open)(vcf) as vcf_content:
        for line in vcf_content:
            if type(line) is str:
                pass
            else:
                line = line.decode('UTF-8')
            if line[0] != '#':
                if line.count('\t') < 9:
                    sys.exit("{} is not a correct .VCF file".format(bed))
                if "\"" in line:
                    line = line.replace("\"", "")
                variant = pythonvcf.Variant(line)
                variant_in_range = is_variant_in_ranges(variant, bed_ranges)
                if variant_in_range:
                    print(line[:-1])

if __name__ == "__main__":
    main()


