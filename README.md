# ClinSeq
Automatic Annotation and deposition of UW Virology clinical sequencing data

This program takes a fasta file full of as many viruses as you would like to annotate and then annotates them and packages them nicely up for genbank submission.

# Installation
Installation should be very straight forward, just clone repository and change locataion of BLAST_DB
You'll also probobly want to just create a blank .csv file for metadata so the program will skip it

# Usage
Go into the python file and just set the "fasta_loc" to whatever fasta you want - you should make the names of the sequences (the things after >) what you would like the strain of the virus to be. Then hit run! Your .sqn files should each be in a folder named what you named your sequences. 

Right now we accomplish the annotation by searching blast for the best hit that is a complete genome and using a maaft alignment to generate annotations for the supplied sequence. We can handle some ribosomal slippage as well as RNA editing in HPIV/Measles/Mumps with support for more and more viruses as I end up having to deal with them. 
