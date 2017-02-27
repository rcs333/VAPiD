# ClinSeq
Automatic Annotation and deposition of UW Virology clinical sequencing data
This repository contains a few different programs

# Sorter
Annotates clinical virus sequences and creates genbank approved files for submission. The only information you need to supply is a fasta with as many viruses as you'd like to annotate in it. Requires BLAST+ and maaft 

# Virus_annotate
Searches a file with many fastas for ORFs and does gene prediction here - basically prokka in batch mode that automatically will blast and HHPRED every single predicted viral gene. Requires prokka (in -c mode) and BLAST+ and HHPRED

# Installation
Installation should be straightforward - It's been tested on unix and mac

1. Install all dependincies (Shoutout to the wonderful people who wrote these!)
  BLAST+ https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
  maaft http://mafft.cbrc.jp/alignment/software/
  prokka http://www.vicbioinformatics.com/software.prokka.shtml
  HHPRED https://github.com/soedinglab/hh-suite
  
1. Clone repository and change locataion of BLAST_DB (provided you have a local setup) 

2. You'll also probobly want to just create a blank .csv file for metadata so the program will skip it - the variable for this is metadata_sheet_location

3. Put a template.sbt and comment.cmt file into the folder that this program will be running from - This is necessary for generating a useable .sqn file that you can send to genbank

# Installation - Virus_annotate

1. You'll need to set the path to the HHPRED database as well as blast and prokka
2. Make two direcotries in the folder that you'll be running from called 'prokka_output' and 'temp_contig_dir'

# Usage - Sorter
Go into the python file and just set the "fasta_loc" to whatever fasta you want - you should make the names of the sequences (the things after >) what you would like the strain of the virus to be. Then hit run! Your .sqn files should each be in a folder named what you named your sequences. 

Right now we accomplish the annotation by searching blast for the best hit that is a complete genome and using a maaft alignment to generate annotations for the supplied sequence. We can handle some ribosomal slippage as well as RNA editing in HPIV/Measles/Mumps with support for more and more viruses as I end up having to deal with them. 

# Usage - Virus_annotate
To annotate a scaffold simply call 'python virus-annotate 'scaffold-path' 'minimum base cutoff' 'choose a name' 
You can also call it with the -h option to get a not very detailed description of these arguments.

If minimum base cutoff is passed a zero the program will take a pre assembled segmented genome and attempt to annotate and send it to tbl2asn

This generates a folder with the name 'choose a name'
Inside the folder will be all of the different contigs that were over the minimum base cutoff
Inside each of these folders will be the output of tbl2asn as well as an hhpred and blastp result for any ORF that prokka couldn't figure out

You would then have to go through and use your human brain to look at the results and figure if you want to rename anything.
After manually doing this you would have to run tbl2asn again yourself, although I'm going to at some point cook up a quick script to do this automatically.

# Future directions
Eventually I'll make sorter a command line tool but as I'm developing it just using it in my IDE is much better.
Eventually planning on removing the need to generate fasta files from NGS reads for sorter - looking at using centrifuge to find best reference sequence and then performing a directed alignment on that
