maf-cut
=======

``maf-cut`` cuts out parts of MAF_ format alignments.  For example::

  maf-cut chr7:1234000-1235000 alignments.maf

This gets alignment parts within the range 1234000-1235000 of the
sequence named "chr7".

The ranges are zero-based.  For example, ``chr7:0-1000`` means the
first 1000 bases of chr7.

You can omit the sequence name::

  maf-cut 1234000-1235000 alignments.maf

Then, the range applies to the topmost sequence in each alignment.

.. _MAF: http://genome.ucsc.edu/FAQ/FAQformat.html#format5
