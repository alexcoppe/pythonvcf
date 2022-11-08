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
                     variant = pythonvcf.Variant_with_genotype(line)
                     chromosome_in_variant = variant.chromosome
                     position_in_variant = variant.position
                     clnsig = line[int(match.start()):int(match.end())]
                     clnsig_value = clnsig[:-1].split("=")[1]

                     cur = conn.cursor()
                     query = "SELECT * FROM CLINVAR WHERE CHR == {} AND POSITION == {}  LIMIT 1;".format(chromosome_in_variant[3:], position_in_variant)
                     cur.execute(query)
                     rows = cur.fetchall()

                     for row in rows:
                         print(row)
 
if __name__ == "__main__":
    main()
