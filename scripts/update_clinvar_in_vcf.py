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

def update_clinvar_column(original_line_from_vcf, field, clinvar_entry):
    clinvar_match = re.search("{}=.+?;".format(field), original_line_from_vcf)
    d = {"CLNREVSTAT": 7}
    #print(clinvar_entry)
    if clinvar_match:
        start = clinvar_match.start()
        end = clinvar_match.end()
        original_match  = original_line_from_vcf[int(clinvar_match.start()):int(clinvar_match.end())]
        new_match = clinvar_entry[d.get("CLNREVSTAT")]
        first_part_original_line_from_vcf = original_line_from_vcf[0:start]
        last_part_original_line_from_vcf = original_line_from_vcf[end:len(original_line_from_vcf)]
        new_line = first_part_original_line_from_vcf + field + "=" + new_match + ";" + last_part_original_line_from_vcf
        #print("original: ", original_line_from_vcf)
        #print("new: ", new_line)
        return new_line

def update_variant_from_vcf(original_line_from_vcf, start, end,  clinvar_entry):
    first_part_original_line_from_vcf = original_line_from_vcf[0:start]
    last_part_original_line_from_vcf = original_line_from_vcf[end:len(original_line_from_vcf)]
    new_line = first_part_original_line_from_vcf + "CLNSIG=" + clinvar_entry[6] + ";" + last_part_original_line_from_vcf

    return update_clinvar_column(new_line, "CLNREVSTAT", clinvar_entry)[:-1]

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
                    cur.execute(query)
                    rows = cur.fetchall()

                    if len(rows) == 0:
                        print(line[:-1])
                    else:
                        # line comes from the studied VCFs
                        # rows come from the Clinvar tab separated file
                        description_from_vcs = line[match.start():match.end()].split("=")[1][:-1]
                        descrition_from_clinvar = rows[0][6]
                        if description_from_vcs != descrition_from_clinvar:
                            new_variant = update_variant_from_vcf(line, match.start(), match.end(), rows[0])
                            print(new_variant)
            else:
                print(line[:-1])
                    
 
if __name__ == "__main__":
    main()
