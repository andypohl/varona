Usage
=====

Below is the usage information for the varona command:

.. code-block :: text

  usage: varona [-h] [--log-level {debug,info,warning,error}] 
    [--assembly {GRCH37,GRCH38}]
    [--maf {FR,SAMPLES,BCFTOOLS}]
    [--no-vep] [--vep-data VEP_DATA] [--version]
  
  Annotate a VCF file.
  
  options:
    -h, --help            show this help message and exit
    --log-level {debug,info,warning,error}
                          Set the logging level (default: WARNING)
    --assembly {GRCH37,GRCH38}
                          genome assembly used in Ensembl VEP API (default: GRCh37)
    --maf {FR,SAMPLES,BCFTOOLS}
                          MAF calculation method (default: SAMPLES)
    --no-vep              Skip VEP API querying (no effect if --vep-data is provided)
    --vep-data VEP_DATA   Path to VEP output file
    --version             Show program's version number and exit

The ``--maf`` option gives three options for calculating a minor allele
frequency (MAF) value to include in the output CSV file.  The BCFTOOLS option
is available if Varona can find an installation of `bcftools <https://samtools.github.io/bcftools/bcftools.html>`_
on the system.
