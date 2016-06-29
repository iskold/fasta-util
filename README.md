# fasta-util
This repository contains some quick and dirty scripts I've found useful. I keep them here for backup purposes.
Contains the following:

* downsampler.py - (Pseudo)-Randomly downsamples fastq-files into a specific target read count. NOTE: Not normalization! If you want normalization, you need to find a normalizer, for example based on kmer frequencies, or map depth etc. Downsampling is useful if you care about read abundances and want to compare several samples that each have been sampled to different depths.
* shiftfasta.py - Shifts single-entry fasta files certain number of bp/aa up or down. Useful if you want to compare two sequences that do not start from the same gene.
* statidba.py - Lists the 10 largest fasta entries in a fasta file, or any number you want. Not memory-efficient (yet).
* reversecompliment.py - Reverse complements fasta entries in a nucleotide fasta file.
* firstfasta.py - extracts first fasta entry from a fasta file with multiple entries.
