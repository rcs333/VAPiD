import timeit
import subprocess
import glob
import sys
import argparse

start = timeit.default_timer()

# This program runs some shit and does some shit about clinical virus samples
# Gonna write more as I need too

# parser = argparse.ArgumentParser(description= 'Annotate a set of UW clinical viral samples, pulling virus information from prokka and blast')
# parser.add_argument('file_dir', help='Input file directory, all .fasta files will be processed and .seq and .gbf files will be produced in the format input_dir/output/FASTA_name')
# parser.add_argument('metadata_info_sheet_location', help='.csv file where all of the metadata is stored')
# parser.add_argument('sbt_file_loc', help='location of .sbt file for .gbf file creation')

# args = parser.parse_args()

# Here I assume that the .fasta file has multiple fastas as opposed to being given a directory, this is subject to later change
fasta_filename = '10fasta_UWViroClinSeq.fasta'
metadata_info_sheet = 'UWVIROCLINSEQ - SCCA.csv'
gff_file_loc = 'HPIV3_121416.gff'

# Takes the name of a clincical virus as specified on the metadata sheet and returns a list of the relevant metadata
def pull_metadata(virus_name):
    for line in open(metadata_info_sheet):
        if line.split(',')[1] == virus_name:
            # Parse and steal input
            # reutrn two strings, one for the cmt file and the other for the .fsa features


def parse_gff(gff_file_loc):
    # First two lines are garbarge
    # One line a sequence format: ##TYPE DNA virus_name
    # then sequences start:
    # FORMAT:
    # RNA NAME
    # SEQUENCE
    # end-
    # all of them, also in the same order as the first list
    # NAME GENEIOUS cds ## ## stupid shit then the names
    # all named, and also in order
    # Write this into lists
    # write the damn files right here
    # pull_metadata(name)
    # write the .tbl and .fsa right here

def write_output():
    # make a folder for each, name it the sample name
    # Go through and make .fsa and .tbl files out of our data

# TODO: generalize, but first I'mma run it with hard coded filepaths

def run_tbl():
    # run .tbl2asn on all of the folders and process the .sqn files for submission
    # Probobly entails throwing the .sbt file into each folder
    #


# Process the fasta_file




# Now we go through and actually work our magic on the viruses
for x in range(0,len(virus_name_list)):
    clin_data_list = pull_metadata(virus_name_list[x])
    # TODO: Modify fasta/cmt file
    # TODO: Run Prokka - with options stolen from sheet



