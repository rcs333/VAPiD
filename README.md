# ClinSeq
Automatic Annotation and deposition of UW Virology clinical sequencing data
This repository contains scripts for analysis of assembled fasta files of known human viruses. This program annotates these viruses and then packages them up for NCBI deposition. This repository is under VERY active development and I also use the code here very frequently. Although now I'm starting to develop in the dev branch - I'll try to keep the main branch from breaking so at all times someone could come on here and download working code

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
  
2. Clone repository and change location of BLAST_DB (line 13 of ClinVirusSeq.py) (provided you have a local setup) if you don't then you can edit the blast command to be remote - this is almost always faster unless you're running from a pretty powerful workstation or server. Actually the current version uses remote by default so all you need is an internet connection and the blast utilities.

3. You'll also probably want to just create a blank .csv file for metadata so the program will skip it - the variable for this is metadata_sheet_location. Although if you have metadata (coverage is probably the most common) you could load it into a csv and point the program there and to the correct column. (line 212 for the column) Collection date is also pulled from this .csv file. You'd want the csv file to have in the first column the strain name and then collection date in the 6th column and coverage in the 5th (assuming we're calling the first column column 1)

4. I guess because all of the paths are hardcoded because it makes it easier for me to actually use this program to get it working yourself you'll need to find all the subprocess commands and edit the locations of the programs to where they are on your actual machine. For example I have blastn in ct-test/ncbi-blast-2.6.0+/bin/blastn (relative to the directory that you're running the script from) you'd need to change this path to where blast is installed on your personal system. (and do this for the other programs) Line 43, Lines 71-82, 85, and Line 321 of ClinVirusSeq.py  (Of course this wouldn't be a problem if you had installed these correctly on your machine - then you could just remove the direct paths and put in the program name)

4. Put a template.sbt file into the folder that this program will be running from - This is necessary for generating a usable .sqn file that you can send to genbank. You can generate these somewhere on the NCBI site - I made one once and haven't changed it since.

5. I guess you also need to make sure you have read/write permissions in the folder you run this from as this script utilizes subprocess and bash commands like mv and cp like way too much.

# Usage - ClinVirusSeq.py
Create your fasta file with all of the sequences that you would like to annotate - you should make the names of the sequences (the things after >) what you would like the strain of the virus (In your NCBI Genbank records) to be. 

Then you would need to run the ClinVirusSeq.py script from the command line. i.e. cd to the directory that this is living it - if you cloned from github it'll be ClinVirusSeq/ 

You can run it with no arguments or the -h flag and it'll print out some usage information. But basically the way you run this is to type (without quotes) "ClinVirusSeq.py fasta_file metadata_info_sheet sbt_file_loc" where fasta_file is the location of all your sequences, metadata_info_sheet is the location of the semi-optional metadata sheet, and sbt_file_loc is the location of your .sbt file. I have it set to take a specified .sbt file so that you can create them with different author lists if you need to to make it easier.

Running the script without the -r flag all of your records will be in individual folders named by "strain" (what came after > in your fasta file). As a side note - please don't put whitespace or backslashes in your strain names. I would recommend doing this if you know you'll need to run this several times due to bugs or other issues. For example if you're debugging otherwise just use -r the first time

If you're feeling bold just put the -r flag on the first time and all your .sqn files will get consolidated. Your output will be in a folder with the date YEAR_MONTH_DAY/virus_species/ The script creates consolidated .sqn files by virus species because that's what the lovely people at Genbank asked us to do. So you'll have a date folder for each fasta file - then another sub folder for each virus species and inside these folders will be species.sqn which should be what you wanna submit to NCBI. You can find individual records inside the subfolders (named by strains)  (so if the first line of your fasta file is >SC12309 then there will be a folder called SC12309 with all of the .gbf .tbl and ect files that you need and love) Note that this does not (yet) protect against the same virus blasting to separate records that other people (sometimes me honestly) named two different things. (This is actually a pretty hard problem to solve without hard coding every virus, which I may end up doing)

Also - any records that had stop codons will not be moved into the species folders - they'll be left out in your main directory for you to examine - most often this comes from errors with coverage but sometimes it's because I introduced a catastrophic bug.
tbl2asn will also spit out some errors and warnings onto the command line which are obviously pretty helpful. One known bug is that if you pass a sequence that doesn't contain every orf in its entirety you'll get start codons of 1< which obviously tbl2asn doesn't like

# Implementation Details and Important Notes

Right now we accomplish the annotation by searching blast for the best hit that is a complete genome and using a maaft alignment to generate annotations for the supplied sequence. Basing the annotations off the best blast hit. We can handle some ribosomal slippage as well as RNA editing in HPIV/Measles/Mumps/Sendai with support for more and more viruses as I end up having to deal with them. Coronavirus 229E and HPIV3 are annotated off hardcoded records due to variability in accuracy for these records in blast.

------------------------------------------------------------------------------------------------------------------------------------------
IMPORTANT SECURITY NOTE: DO NOT PUT THIS CODE ON A SERVER AND LET IT FACE UNKNOWN CLIENTS!!! 
This code uses direct subprocess commands with user supplied input so it would be trivial to insert bash commands into the .fasta file and completely destroy a server. I don't see any reason why anyone would do this but if you are putting it on a server make sure that only people you know won't destroy your server will use it. Also - I guess don't name your strains really stupid things like & rm -rf / or & :(){ :|: & };: I doubt this would work as I've written it so you'll be fine as long as you name your strains normal alphanumeric strings or normal human words. Just don't let potentially malicious users run this code.
------------------------------------------------------------------------------------------------------------------------------------------

# Future directions
Currently developing a way of going directly from NGS reads coming straight off the Illumina to reference, assemble off the reference and then annotate, although this is in no way working right now stay tuned! :) Also going to try to have some way of consolidating gene product annotations because NCBI was complaining about that and I agree it's kinda bad and right now I'm manually reviewing all my annotations, which is slow.
