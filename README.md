# ClinSeq
Automatic Annotation and deposition of UW Virology clinical sequencing data

Annotates clinical virus sequences and creates genbank approved files for submission. The only information you need to supply is a fasta with as many viruses as you'd like to annotate in it.

# Installation
Installation should be very straight forward.
1. Clone repository and change locataion of BLAST_DB

2. You'll also probobly want to just create a blank .csv file for metadata so the program will skip it - the variable for this is metadata_sheet_location

3. Put a template.sbt and comment.cmt file into the folder that this program will be running from - This is necessary for generating a useable .sqn file that you can send to genbank

# Usage
Go into the python file and just set the "fasta_loc" to whatever fasta you want - you should make the names of the sequences (the things after >) what you would like the strain of the virus to be. Then hit run! Your .sqn files should each be in a folder named what you named your sequences. 

Right now we accomplish the annotation by searching blast for the best hit that is a complete genome and using a maaft alignment to generate annotations for the supplied sequence. We can handle some ribosomal slippage as well as RNA editing in HPIV/Measles/Mumps with support for more and more viruses as I end up having to deal with them. 

# Future directions
Eventually I'll make this a command line tool but as I'm developing it just using it in my IDE is much better.
Looking into removing the need to generate fasta files from NGS reads
