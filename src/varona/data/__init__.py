"""Data module for Varona.

This module just has a .fai added from the 1000 genomes project.
It aids in the creation of fake VCF files for testing.

- ``verona.data.human_g1k_v37.fasta.fai``, downloaded from `1000 genomes FTP <http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/>`_.

Examples:
    .. code-block:: python

        import importlib.resources as pkg_resources
        from varona import data

        with pkg_resources.open_text(data, "human_g1k_v37.fasta.fai") as fai:
            for line in fai:
                print(line.strip())

"""
