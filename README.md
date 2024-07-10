# varona
Variant Annotator.

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
