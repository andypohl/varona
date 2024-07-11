Usage
=====

Below is the usage information for the varona command:

.. code-block :: text

    usage: varona [-h] [--log-level {debug,info,warning,error}] 
        [--assembly {GRCh37,GRCh38}] 
        [--maf {FR,BCFTOOLS,SAMPLES}] input_vcf output_csv
    
    Annotate a VCF file.
    
    positional arguments:
      input_vcf             Path to the input VCF file
      output_csv            Path to the output CSV file
    
    options:
      -h, --help            show this help message and exit
      --log-level {debug,info,warning,error}
                            Set the logging level (default: WARNING)
      --assembly {GRCh37,GRCh38}
                            genome assembly used in Ensembl VEP API (default: GRCh37)
      --maf {FR,BCFTOOLS,SAMPLES}
                            MAF calculation method (default: SAMPLES)

The ``--maf`` option gives three options for calculating a minor allele
frequency (MAF) value to include in the output CSV file.  The BCFTOOLS option
is available if Varona can find an installation of `bcftools <https://samtools.github.io/bcftools/bcftools.html>`_
on the system.
