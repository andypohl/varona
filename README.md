![Python Version from PEP 621 TOML](https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2Fandypohl%2Fvarona%2Fmain%2Fpyproject.toml&logo=python)
![main branch unittests](https://github.com/andypohl/varona/actions/workflows/unittest.yml/badge.svg)
[![codecov](https://codecov.io/gh/andypohl/varona/graph/badge.svg?token=Bdlakar3V6)](https://codecov.io/gh/andypohl/varona)

# varona
Variant Annotation library and command line tool to read and annotate VCF
file records, supplementing the information with data from the Ensembl VEP
API.  More information about the Varona is available in at the
[documentation site](https://varona.pages.dev/).

## Installation

It's recommended to install `varona` into a Python virtual environment.

```bash
pip install --extra-index-url https://pypi.pohl.io/simple/ varona
```

## Usage

```bash
varona input.vcf output.csv
```

Additionally, there are some options that can be used:

```bash
varona --help
## usage: varona [-h] [--log-level {debug,info,warning,error}]
##     [--assembly {GRCh37,GRCh38}]
##     [--maf {FR,BCFTOOLS,SAMPLES}] input_vcf output_csv
## 
## Annotate a VCF file.
## 
## positional arguments:
##   input_vcf             Path to the input VCF file
##   output_csv            Path to the output CSV file
## 
## options:
##   -h, --help            show this help message and exit
##   --log-level {debug,info,warning,error}
##                         Set the logging level (default: WARNING)
##   --assembly {GRCh37,GRCh38}
##                         genome assembly used in Ensembl VEP API (default: GRCh37)
##   --maf {FR,BCFTOOLS,SAMPLES}
##                         MAF calculation method (default: SAMPLES)
## 
```
