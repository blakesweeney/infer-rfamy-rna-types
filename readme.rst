==========================
Improving Rfam's RNA Types
==========================

This is some code to try to improve the ISNDC type that Rfam families get.

Logic
-----

The script in `bin/infer.py` will infer an INSDC RNA type for all Rfam
families. It has a series of different methods for selecting the name. In the
order they are attempted:

manual
  Check if there is a manually assigned RNA type for the family. Currently
  there are none.

name
  Check if the name contains a descriptive string like ``Y_RNA``. 

so-term
  Check if the SO terms assigned to the family have an RNA type assigned.

rna-type
  Use Rfam's RNA type to infer an INSDC RNA type.

so-search
  Use the given SO terms and look up the tree to find the nearest parent SO
  term(s) that have been assigned an RNA type.

The assignments to RNA types are found in ``data/manual-assignments.json``. Most
families have an inferred type (though ``other`` is a common one) with this
method. The last stage, searching the attached SO terms, is a minor
contribution to the inferred names, and can probably be ignored in general.

Also note that in theory it is possible (and correct) to infer more than one
ISNDC RNA type for some families. For example family []() is a snoRNA that can
act as an miRNA, so both types are correct. Here we use as the RNA type.

Results
-------

The inferred types are in ``data/infered-types.txt``. A few quick summaries:

.. code:: shell

  $ xsv select method data/infered-types.txt | sort | uniq -c
      479 fallbacks
        1 method
       36 name
     1056 rna-type
       15 so-search
     1002 so-term
  $ xsv select rna_types data/infered-types.txt | sort | uniq -c
      479 ""
        2 RNase_MRP_RNA
        5 RNase_P_RNA
        9 SRP_RNA
        4 Y_RNA
       33 antisense_RNA
       10 autocatalytically_spliced_intron
        5 hammerhead_ribozyme
      217 lncRNA
      484 other
      530 precursor_RNA
       15 rRNA
       12 ribozyme
        1 rna_types
       18 snRNA
      753 snoRNA
        2 tRNA
        4 telomerase_RNA
        5 tmRNA
        1 vault_RNA


The method called ``fallbacks`` means no INSDC type could be found.

Requirements
------------

- python3
- make
- wget
- gzip

You must also install the requirments in ``requirements.txt``.

.. code:: shell

  $ pip install -r requirements.txt

Usage
-----

To produce a csv of family to INSDC RNA types in ``data/infered-types.txt``.

.. code:: shell

  $ make
