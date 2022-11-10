#!/usr/bin/env python

import argparse
import gzip
import os
import sys
import pythonvcf
import re
import sqlite3

def create_connection(db_file):
    conn = None
    try:
        conn = sqlite3.connect(db_file)
    except Error as e:
        print(e)
    return conn

def update_variant_from_vcf(original_line_from_vcf, start, end,  new_clnsig):
    first_part_original_line_from_vcf = original_line_from_vcf[0:start]
    last_part_original_line_from_vcf = original_line_from_vcf[end:len(original_line_from_vcf)]
    print(original_line_from_vcf)
    print(first_part_original_line_from_vcf + "CLNSIG=" + new_clnsig + ";" + last_part_original_line_from_vcf)

def main():
    parser = argparse.ArgumentParser(description="Update CLINVAR data in a VCF")
    parser.add_argument('-v', '--vcf', action='store', type=str, help="The vcf to be parsed", required=True)
    parser.add_argument('-s', '--sqlitedb', action='store', type=str, help="The path of the SQLITE db to search updates", required=True)
    args = parser.parse_args()

    vcf = args.vcf
    sqlitedb_path = args.sqlitedb

    conn = create_connection(sqlitedb_path)

    with (gzip.open if vcf.endswith(".gz") else open)(vcf) as vcf_content:
        for line in vcf_content:
            if type(line) is str:
                pass
            else:
                line = line.decode('UTF-8')
            if line[0] != '#':
                match = re.search("CLNSIG=.+?;", line)
                if match:
                    #print(match)
                    variant = pythonvcf.Variant_with_genotype(line)
                    chromosome_in_variant = variant.chromosome
                    position_in_variant = variant.position
                    alternative_in_variant = variant.alternative
                    reference_in_variant = variant.reference
                    clnsig = line[int(match.start()):int(match.end())]
                    clnsig_value = clnsig[:-1].split("=")[1]
                    #print(clnsig_value)

                    cur = conn.cursor()
                    query = "SELECT * FROM CLINVAR WHERE CHR == {} AND POSITION == {} AND ALTERNATIVE == '{}' AND REFERENCE == '{}' LIMIT 1;".format(chromosome_in_variant[3:], position_in_variant, alternative_in_variant, reference_in_variant)
                    #print(query)
                    cur.execute(query)
                    rows = cur.fetchall()

                    if len(rows) == 0:
                        pass
                        #print(match)
                    else:
                        # line comes from the studied VCFs
                        # rows come from the Clinvar tab separated file
                        description_from_vcs = line[match.start():match.end()].split("=")[1][:-1]
                        #print(description_from_vcs)
                        descrition_from_clinvar = rows[0][6]
                        if description_from_vcs != descrition_from_clinvar:
                            update_variant_from_vcf(line, match.start(), match.end(), descrition_from_clinvar)
                            #print(line[match.start():match.end()])
                            #print(rows[0][6])
                            #print("\n")
                    
 
if __name__ == "__main__":
    main()
