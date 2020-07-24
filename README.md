This folder contains:


* program cpgoe.pl

this program calculates the % of each base (A,C,G,T,N) and the ratios CpG (observed/expected) and GpC (observed/expected) for each sequence in a fasta file (ratios are calculated for a whole sequence).

It uses the transliteration operator (tr) in pattern matches for a rapid count of occurences of, e.g., "CG" sites.

This programs needs perl and BioPerl (using the Bio::Seq library)

run with: 

cpgoe.pl -f  file.fas
