import timeit
import subprocess
from random import randint
import argparse

# This program takes in a gff file, the location of a .csv metadata file (IN A SPECIFIC FORMAT) and an .sbt file
# Containing author names, mainly, and produces folders with the .gbf and .sqn files for genbank submission as well as
# tbl2asn error checking outputs

start = timeit.default_timer()

parser = argparse.ArgumentParser(description='Package a set of UW clinical virus sequences for submission, pulling '
                                             'virus name information from blast and annotations are contained inside'
                                             'the .gff file passed to the script originally')
parser.add_argument('gff_file', help='Input file in .gff format, should contain annotations and genome information for'
                                     'each virus to be processed')
parser.add_argument('metadata_info_sheet', help='The metadata sheet that contains whatever we have on these samples')
parser.add_argument('sbt_file_loc', help='File path for the .sbt file that should contain author names mainly')

args = parser.parse_args()


# Takes the name of one of the virus samples and returns a list of all the metadata we have in the sheet
def pull_metadata(virus_name):
    for line in open(args.metadata_info_sheet):
        if line.split(',')[1] == virus_name:
            # Just returns all of the data in a list (this is a csv so split(',') just pulls it all
            return line.split(',')


# Takes the location of a .gff file with one or more virus genomes and annotations and two lists and a map the lists
# are of the viral strain names and the genomes in the same order so [x] in each will be the same, the map is a map of
# strain names to annotations
def parse_gff(gff_location):
    read = False
    virus_name_list = []
    dna_string = ''
    dna_list = []
    annotation_map = {}

    for line in open(gff_location):
        # if we're in the first half which has names and genomes then we need to grab all the names and all the genomes
        # while keeping them in order
        if line[0] == '#':
            if not read:
                virus_name_list.append(line.split()[1])
                read = True
            elif line[2:5] != 'end':
                # Grab the genetic information and force it to be DNA, the 2:-1 is so that we ignore the ## at the
                # start and the \n at the end of the line
                dna_string += line[2:-1].replace('U', 'T')
            else:
                dna_list.append(dna_string)
                dna_string = ''
                read = False

        # We're only interested in the lines that actually have useable annotation information on them
        elif line.split()[2].lower() == 'cds' or 'mat_peptide':
            sample_name = line.split()[0]
            # if we already have at least one annotation for this virus then add to the string the next seperated by *
            # this way we can use the split('*') command on the
            if sample_name in annotation_map:
                annotation_map[sample_name] += '*' + ' '.join(line.split()[2:])
            else:
                annotation_map[sample_name] = ' '.join(line.split()[2:])
    return
