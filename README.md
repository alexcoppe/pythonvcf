# pythonvcf

### Available scripts


| Script        | Description|
| ------------- |:-------------|
| parse_clinvar_vcf.py | Creates a tab separated file from the VCF used as input |
| update_clinvar_in_vcf.py | Updates a VCF file annotated with [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/) using a [SQLite](https://www.sqlite.org/index.html) made from the table produced by parse_clinvar_vcf.py script |


### Installation

pip install -e ./

### Uninstallation

- pip uninstall pythonvcf.py

- rm -rf scripts/pythonvcf.py.egg-info/
