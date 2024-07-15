Goals
=====

.. raw:: html

    <style type="text/css">
      span.bolditalic {
         font-weight: bold;
         font-style: italic;
      }
    </style>

.. role:: bolditalic
    :class: bolditalic

The goal of the Varona library is to provide a simple command-line utility,
and associated Python library to annotate variants in a VCF file. 

The CSV output tries to consolidate the following information:

**1. Depth of sequence coverage at the site of variation.**

    This is found within a variant's TC field within the INFO section. In the CSV, it is the "sequence_depth" column.

**2. Number of reads supporting the variant.**

    This is found within a variant's TR field within the INFO section. When there are
    more than one alternative allele present in the variant, the maximum of the values
    in the TR is taken. This is why the corresponding field in the CSV file is named
    "max_variant_reads".

**3. Percentage of reads supporting the variant versus those supporting reference reads.**

    This is calculated as "max_variant_reads" divided by "sequence_depth" to obtain
    a value between 0 and 1, then multiplied by 100 to obtain a percentage.

**4. Using the VEP hgvs API, get the gene of the variant, type of variation (substitution,insertion, CNV, etc.) and their effect (missense, silent, intergenic, etc.).  The API documentation is** `here <https://rest.ensembl.org/#VEP>`_.

    Columns for this data in the CSV are named: "gene_id" or "gene_name", "type", and 
    "effect", respectively. Sometimes the values are missing in the CSV file if there
    isn't an associated transcript.

**5. The minor allele frequency of the variant if available.**

    There are several ways to obtain this data and Varona's default is to calculate
    a MAF from the sample genotypes in the VCF file. This is the same strategy used
    by bcftools, and if bcftools is installed, Varona can also use "bcftools 
    +fill-tags" to add that field to the INFO section of each variant line.  
    Unfortunately, the bcftools supplied by Pysam does not have the plugins like 
    "fill-tags" compiled-in, so if this functionality is desired, bcftools needs to
    be installed separately.  The CSV file will have a column named "maf". A third
    way to obtain the MAF is to use the "FR" field in the INFO section in the
    same way as a "AF" fields are used. This is not the default behavior though,
    because in a small number of cases, the resulting MAF is not the same as the
    value produced by bcftools.

**6. Any additional information seeming relevant.**

    Varona keep positional and allele fields from the variant, and the "transcript_id" field, if available.
