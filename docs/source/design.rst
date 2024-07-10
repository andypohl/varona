Design
======

===
CLI
===

A command-line application, also named ``varona`` is included to execute
the main functionality of the library.  The command-line application has a few
options to select including how to calculate the minor allele frequency,
which genome assembly to use for the Ensembl API (`GRCh37 <http://grch37.rest.ensembl.org/>`_,
or `GRCh38 <http://rest.ensembl.org/>`_, with GRCh37 being the default),
and the verbosity of the logging output.


===
API
===

A major design goal with Varona is to have a flexible and extensible
programmatic interface. Because the output of Varona is a CSV file, the natural
progenitor of a CSV file is a DataFrame.  Varona will create two types of
DataFrames and join the two together: (a) a DataFrame based on the data embedded
in the VCF file only, and (b) a DataFrame based on the data obtained from the
Ensembl API.  Many times when a CSV file is created for someone, they will
add or remove columns to suit their needs. The exact data extracted from the
VCF records or the API can be adjusted by supplying alternative callback
functions, each of which returns a list of dictionaries each corresponding to
a DataFrame row. A different approach to callback functions could've been to
define a base class and let the programmer define subclasses but for this
version, the callbacks should provide a fair amount of flexibility.
