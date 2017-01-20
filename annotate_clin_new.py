import timeit
import subprocess
from random import randint
import argparse
import os

"""
Ryan Shean
Clinical Virus Annotator
1/19/17

This program takes in a gff file, the location of a .csv metadata file (IN A SPECIFIC FORMAT) and an .sbt file
Containing author names, mainly, and produces folders with the .gbf and .sqn files for genbank submission as well as
tbl2asn error checking outputs
"""


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
    return virus_name_list, dna_list, annotation_map


# Grab the coverage data from the imported metadata, if there is no coverage information randomly generate it and add
# .1x so that we know that it was randomly generated
def pull_coverage(data_list):
    if data_list[20] != '':
        coverage = data_list[20]
    else:
        coverage = str(randint(20,100)) + '.1x'
    return coverage


# Take the name of a virus sample, and write the .cmt file for it using supplied coverage information
# NOTE: This NEEDS the virus sample directory to have been already created
def write_cmt(sample_name, coverage):
    cmt = open(sample_name + '/assembly.cmt', 'w')
    cmt.write('StructuredCommentPrefix\t##Assembly-Data-START##\n')
    cmt.write('Assembly Method\tGeneious v. 9.1\n')
    cmt.write('Coverage\t' + coverage + '\n')
    cmt.write('Sequencing Technology\tIllumina\n')
    cmt.write('StructuredCommentPrefix\t##Assembly-Data-END##\n')
    cmt.close()


# Takes a sample name and genome information and checks to see if the result file already exists, if it does then we
# just grab it again. Otherwise, blast the virus and pull the strain name out and return it.
def grab_name_of_virus(sample_name, genome_string):
    # Check to see if we've already blasted, there's no need to reblast if the file already exists
    if os.path.isfile(sample_name + '/' + sample_name + '.result'):
        # Write a temporary file that we'll blast to find the name of the virus
        file_to_blast = open(sample_name + '/temp.fasta', 'w')
        file_to_blast.write('>' + sample_name + '\n')
        file_to_blast.write(genome_string)
        file_to_blast.close()
        # Actually do the blasting and save the results in a file sample_name.result
        subprocess.call('blastn -query ' + sample_name + '/temp.fasta -db /Users/uwvirongs/Downloads/surpi-master/nt '
                        '-num_threads 8 -num_descriptions 0 -num_alignments 1 -word_size 28 | tee ' + sample_name + '/'
                        + sample_name + '.result', shell=True)
        # Remove the temporary file to clean up the directory just a little
        subprocess.call('rm ' + sample_name + '/temp.fasta', shell=True)

    # Open the file and grab the important bits then return them
    for line in open(sample_name + '/' + sample_name + '.result'):
        if line[0] == '>':
            strain = line.split('|')[4].split('strain')[0].split('isolate')[0].split(',')[0]
            break
    return strain


# Takes the strain name and the genome of a virus and generates the .fsa file
# NOTE: this can defiantly take more data, but we don't have any on the clinical sequencing database that I'd use right
# now so we don't even use it
def write_fsa(sample_name, virus_strain, virus_genome):
    fsa = open(sample_name + '/' + sample_name + '.fsa', 'w')
    fsa.write('>' + sample_name + ' [organism=' + virus_strain + '] [collection-date=2016] [country=USA] '
              '[moltype=genomic] [host=Human] [gcode=1] [molecule=RNA] [strain=' + sample_name + ']\n')
    fsa.write(virus_genome)
    fsa.close()


def write_tbl(sample_name, virus_annotation, genome):

    tbl = open(sample_name + '/' + sample_name + '.tbl', 'w')
    tbl.write('>Feature ' + sample_name)
    # Just putting a comment here for now
    product_map = {}

    # Go through the annotation string and make a map of products to type$start#end and if there are multiple areas
    # for the same product set them up like type$start#end^start#end
    for entry in virus_annotation.split('*'):
        annotation_type = entry.split()[0]
        start_base = entry.split()[1]
        end_base = entry.split()[2]
        product = entry.split('=')[1].split(';')[0]

        if product in product_map.keys():
            product_map[product] += '^' + start_base + '#' + end_base
        else:
            product_map[product] = annotation_type + '$' + start_base + '#' + end_base

    # go through the product map product by product and write to the .tbl file
    for k in product_map.keys():
        anno_type = product_map[k].split('$')[0].upper()
        times_string = product_map[k].split('$')[1]
        flag = ''
        # If we've got Ribosomal Slippage/RNA Editing as shown by '^'
        if '^' in times_string:
            start_1 = times_string.split('^')[0].split('#')[0]
            end_1 = times_string.split('^')[0].split('#')[1]
            start_2 = times_string.split('^')[1].split('#')[0]
            end_2 = times_string.split('^')[1].split('#')[1]
            tbl.write('\n' + start_1 + '\t' + end_1 + '\t' + anno_type + '\n')
            tbl.write(start_2 + '\t' + end_2 + '\n')
            tbl.write('\t\t\tproduct\t' + k + '\n')
            # TODO: Need someway of changing this to say RNA editing for PHIV(s) perhaps pass the strain and search
            tbl.write('\t\t\texception\tRibosomal slippage')

        else:
            start_int = times_string.split('#')[0]
            end_int = times_string.split('#')[1]
            if end_int == len(genome):
                if ((int(end_int) - int(start_int)) + 1 % 3) == 0 and genome[end_int - 3:end_int].upper() in 'TGA,TAA,TAG':
                    flag = ''
                else:
                    flag = '>'
            tbl.write('\n' + start_int + '\t' + flag + end_int + '\t' + anno_type + '\n')
            tbl.write('\t\t\tproduct\t' + k)
    tbl.close()


# Take a list of metadata, name of a virus, the genetic information, and the annotation and write up all the files
# Then crank them through tbl2asn
def write_a_virus(metadata_list, sample_name, virus_genome, virus_annotation):

    subprocess.call('mkdir ' + sample_name, shell=True)

    coverage = pull_coverage(metadata_list)

    # Create the .cmt and .sbt files in the directory for this virus
    write_cmt(sample_name, coverage)
    subprocess.call('cp ' + args.sbt_file_loc + sample_name + '/', shell=True)

    # Grab the strain name of the virus from blast
    virus_strain_name = grab_name_of_virus(sample_name, virus_genome)

    # Use the strain name to write the .fsa file
    write_fsa(sample_name, virus_strain_name, virus_genome)

    write_tbl(sample_name, virus_annotation, virus_genome)

    subprocess.call('tbl2asn -p ' + sample_name + '/ -t ' + sample_name + '/template.sbt -Y ' + sample_name +
                    '/assembly.cmt -V vb', shell=True)

if __name__ == '__main__':
    start = timeit.default_timer()
    # TODO: add a flag for redoing tbl2asn that ONLY does that - i.e. you could manually edit files then crank em out
    parser = argparse.ArgumentParser(description='Package a set of UW clinical virus sequences for submission, pulling '
                                                 'virus name information from blast and annotations are contained '
                                                 'inside the .gff file passed to the script originally')
    parser.add_argument('gff_file', help='Input file in .gff format, should contain annotations and genome information '
                                         'for each virus to be processed')
    parser.add_argument('metadata_info_sheet', help='The metadata sheet that contains whatever we have on these samples')
    parser.add_argument('sbt_file_loc', help='File path for the .sbt file that should contain author names mainly')

    args = parser.parse_args()

    sample_name_list, genome_list, annot_map = parse_gff(args.gff_file)

    for x in range(0, len(sample_name_list)):
        write_a_virus(pull_metadata(sample_name_list[x]), sample_name_list[x], genome_list[x],
                      annot_map[sample_name_list[x]])