# ClinSeq
Automatic Annotation and deposition of UW Virology clinical sequencing data
This repository contains scripts for analysis of assembled fasta files of known human viruses. This program annotates these viruses and then packages them up for NCBI deposition. This repository is under active development and I also use the code here very frequently. 

If anything breaks, or issues come up, or even if you just need help installing it please feel free to open issues or email me at either rcs333@uw.edu or uwvirongs@gmail.com.

The only information you need to supply is a fasta with as many viruses as you'd like to annotate in it as well as an optional metadata sheet and a .sbt file for author lists in the .sqn and .gbf files that NCBI likes. 

# Installation
Installation should be straightforward - It's been tested on ubuntu and mac

1. Install all dependencies (Shoutout to the wonderful people who wrote these!)

  BLAST+ https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
  
  maaft http://mafft.cbrc.jp/alignment/software/ 
  
  edirect utilities https://www.ncbi.nlm.nih.gov/books/NBK179288/
  
  centrifuge https://github.com/infphilo/centrifuge
  
  tbl2asn https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/
  
2. Clone repository and change location of BLAST_DB (line 14 of ClinVirusSeq.py) (provided you have a local setup) if you don't then you can edit the blast command to be remote - this is almost always faster unless you're running from a pretty powerful workstation or server.

3. You'll also probably want to just create a blank .csv file for metadata so the program will skip it - the variable for this is metadata_sheet_location. Although if you have metadata (coverage is probably the most common) you could load it into a csv and point the program there and to the correct column.

4. I guess because all of the paths are hardcoded because it makes it easier for me to actually use this program to get it working yourself you'll need to find all the subprocess commands and edit the locations of the programs to where they are on your actual machine. For example I have blastn in ct-test/ncbi-blast-2.6.0+/bin/blastn (relative to the directory that you're running the script from) you'd need to change this path to where blast is installed on your personal system. (and do this for the other programs) Line 43, Lines 71-82, and Line 309 of ClinVirusSeq.py

4. Put a template.sbt and comment.cmt file into the folder that this program will be running from - This is necessary for generating a usable .sqn file that you can send to genbank


# Usage - ClinVirusSeq.py
Create your fasta file with all of the sequences that you would like to annotate - you should make the names of the sequences (the things after >) what you would like the strain of the virus to be. 

Then you would need to run the ClinVirusSeq.py script from the command line. i.e. cd to the directory that this is living it - if you cloned from github it'll be ClinVirusSeq/ 

You can run it with no arguments or the -h flag and it'll print out some usage information. But basically the way you run this is to type (without quotes) "ClinVirusSeq.py fasta_file metadata_info_sheet sbt_file_loc" where fasta_file is the location of all your sequences, metadata_info_sheet is the location of the semi-optional metadata sheet, and sbt_file_loc is the location of your .sbt file. I have it set to take a specified .sbt file so that you can create them with different author lists if you need to to make it easier.

Then hit enter! Your .sqn files should each be in a folder named what you named your strains. (so if the first line of your fasta file is >SC12309 then there will be a folder called SC12309 with all of the .gbf .tbl and ect files that you need and love)

The optional -r flag can be used as a quality of line functionality. When included it will just recompile the .sqn and .gbk files using tbl2asn with the default functionality. This is extremely useful for changing author lists or updating organisms names or collection dates after this program has been run. So you can launch this program on a huge number of viruses and then like manually edit the .tbl files or the source information and then re run with the -r flag.

# Implementation Details and Important Notes

Right now we accomplish the annotation by searching blast for the best hit that is a complete genome and using a maaft alignment to generate annotations for the supplied sequence. Basing the annotations off the best blast hit. We can handle some ribosomal slippage as well as RNA editing in HPIV/Measles/Mumps with support for more and more viruses as I end up having to deal with them. 

------------------------------------------------------------------------------------------------------------------------------------------
IMPORTANT SECURITY NOTE: DO NOT PUT THIS CODE ON A SERVER AND LET IT FACE UNKNOWN CLIENTS!!! 
This code uses direct subprocess commands with user supplied input so it would be trivial to insert bash commands into the .fasta file and completely destroy a server. I don't see any reason why anyone would do this but if you are putting it on a server make sure that only people you know won't destroy your server will use it. Also - I guess don't name your strains really stupid things like >rm -rf / or >:(){ :|: & };: I doubt this would work as I've written it so you'll be fine as long as you name your strains normal alphanumeric strings or normal human words. Just don't let potentially malicious users run this code.
------------------------------------------------------------------------------------------------------------------------------------------

# Future directions
Currently developing a way of going directly from NGS reads coming straight off the Illumina to reference, assemble off the reference and then annotate, although this is in no way working right now stay tuned! :)
I'll also be going over and doing some more QC of this code before I start really telling people about it. 
