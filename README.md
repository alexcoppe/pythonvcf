# pythonvcf :dna:

A Python package and a set of script for working with VCFs files.

# Available scripts :man_technologist:


| Script        | Description|
| ------------- |:-------------|
| parse_clinvar_vcf.py | Creates a tab separated file from the VCF used as input |
| update_clinvar_in_vcf.py | Updates a VCF file annotated with [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/) using a [SQLite](https://www.sqlite.org/index.html) made from the table produced by parse_clinvar_vcf.py script |

##  :snake: filter_vcf_by_clinvar_stars.py

A script that filters a VCF file based on ClinVar [gold stars](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/). It needs the VCF to filter, the VCF downloaded from [ClinVar](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/). 

Option | Long Option | What does it do
------------ | ------------- | -------------
-h | --help | Show help
-v | --vcf | The vcf to be filtered
-c | --clinvar | The vcf from ClinVar
-s | --stars | The minimum number of stars to keep

##### Example
```console
>>> filter_vcf_by_clinvar_stars.py -v your_vcf.vcf -c clinvar.vcf -s 3
```

##  :snake: filter_vcf_by_bed.py

This script filters a VCF by using a BED file. It needs the VCF to be filtered and the [.BED](https://en.wikipedia.org/wiki/BED_(file_format)) file with the rangers you want to filter.

Option | Long Option | What does it do
------------ | ------------- | -------------
-h | --help | Show help
-v | --vcf | The VCF to be filtered
-b | --bed | The BED used as filter

##### Example
```console
>>> filter_vcf_by_bed.py -v your_vcf.vcf -b your.bed
```


## Installation

pip install -e ./

## Uninstallation

- pip uninstall pythonvcf.py

- rm -rf scripts/pythonvcf.py.egg-info/
