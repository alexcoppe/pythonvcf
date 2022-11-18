# update_clinvar_in_vcf

This script updates the clinvar fields inserted by ANNOVAR from a SQLite database build from the parse_clinvar_vcf.py script.
The steps to follow are these:

1 - Create the tab separated file from the clinvar downloaded VCF using parse_clinvar_vcf.py

```bash
parse_clinvar_vcf.py -v path_to_the_clinvar_vcf.vcf > tab_separated_file_from_clinvar_vcf.tab
```

2 - Create the SQLite database (just one table: CLINVAR)

```bash
sqlite3 clinvar.db

sqlite> CREATE TABLE CLINVAR (
   ...>     CHR VARCHAR(31) NOT NULL,
   ...>     POSITION INTEGER NOT NULL,
   ...>     ALTERNATIVE VARCHAR(300) NOT NULL,
   ...>     REFERENCE VARCHAR(300) NOT NULL,
   ...>     GENE VARCHAR(20) NOT NULL,
   ...>     CLNSIG VARCHAR(300) NOT NULL,
   ...>     CLNREVSTAT VARACHAR(300) NOT NULL,
   ...>     CLNDISDB VARACHAR(300) NOT NULL,
   ...>     CLNDN VARACHAR(300) NOT NULL
   ...> )
   ...> ;
sqlite> CREATE INDEX id ON CLINVAR (CHR, POSITION, REFERENCE, ALTERNATIVE);
```

3 - Populate the SQLite CLINVAR table

```bash
sqlite3 clinvar.db -separator '\t'  '.mode tabs' '.import clinvar_20221030.tsv CLINVAR' '.exit'
```

4 - Update a VCF file annotated by ANNOVAR 

```
update_clinvar_in_vcf.py -v /path_to_the_annovar_annotated_vcf.vcf -s /path_to_clinvar_sqlite_db.db > /path_to_the_newly_annovar_annotated_vcf.vcf
```
