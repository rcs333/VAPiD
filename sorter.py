# We gon make a small version of our blaster and sorter here that just takes a fasta, cause this is as far back as
# I really want to go for right now

import subprocess
import re
from random import randint
import argparse
import timeit
import os

# TODO: Clean up this whole top section
# Top line is ribosomal slippage
# bottom line is for rna editing, in case I want to handle these guys differently

#hard_list = ['human coronavirus 229E', 'Human coronavirus OC43', 'Human parainfluenza virus 3', 'astrovirus 1', 'west nile', 'HIV', 'HTLV', 'influenza A',
#             'HPIV', 'measels', 'mumps', 'ebola', 'Human parainfluenza virus 4', 'Human parainfluenza virus 2']
#non_specific_hard_list = ['coronavirus', 't-lymphotropic', 'immunodeficiency', 'influenza', 'measels', 'mumps', 'nipah', 'ebola']

# Set up the environment variablesd
fasta_file_loc = 'SC261.fasta'
BLAST_DB_LOCATION = '/Users/uwvirongs/Downloads/surpi-master/nt'

# grab the strain name so we can actually begin to generalize this
#TODO: when we up this to accept fasta with a lot of inputs we'll save the fasta inside the folder as well
#for line in open(fasta_file_loc):
#    if line[0] == '>':
#        strain = line[1:-1]

# make us a directory
#subprocess.call('mkdir -p ' + strain, shell=True)

# TODO: only commented out to save time each step
# blast us this fasta and save the top 25 hits into a file
#subprocess.call('blastn -query ' + fasta_file_loc + ' -db ' + BLAST_DB_LOCATION + ' -num_threads 8 -num_descriptions 0 '
#                '-num_alignments 25 -word_size 28 | tee ' + strain + '/' + strain + '.blastresults', shell=True)

# read through the top 25 hits and save the top hit that has the word complete in it - if we don't do this fail miserably


def read_fasta(fasta_file_loc):
    strain_list = []
    genome_list = []
    dna_string = ''
    for line in open(fasta_file_loc):
        if line[0] == '>':
            strain_list.append(line[1:-1])
            if dna_string != '':
                genome_list.append(dna_string)
                dna_string = ''
        else:
            dna_string += line[:-1]
    genome_list.append(dna_string)
    return strain_list, genome_list


def blast_n_shit(strain, our_fasta_loc):
    # If we've already done this one before skip the blasting step, should speed up error chekcing in the futuroe
    if not os.path.isfile(strain + '/' + strain +'.blastresults'):
        subprocess.call('blastn -query ' + our_fasta_loc + ' -db ' + BLAST_DB_LOCATION + ' -num_threads 8 -num_descriptions 0 '
                        '-num_alignments 25 -word_size 28 | tee ' + strain + '/' + strain + '.blastresults', shell=True)

    # read through the top 25 hits saved earlier and save the acsession number of the best hit that's complete
    read_next = False
    for line in open(strain + '/' + strain + '.blastresults'):
        # This line removed because blast is mostly screwy and I'm pretty sure we're garunteed complete genomes here
        if line[0] == '>': # and 'complete' in line:
            name_of_virus = line.split('|')[4].split('strain')[0].split('isolate')[0].strip()
            ref_seq_gb = line.split('|')[3]
            if 'complete' in line:
                break
            else:
                read_next = True
        elif read_next:
            if 'complete' in line:
                break
            else:
                read_next = False
    # This skips us the fact that silly genbank put a laboratory strain as the ref_seq and we get clinicals
    # Should become obsolete when we own the coronaviruses
    if 'CORONAVIRUS 229E' in name_of_virus.upper():
        ref_seq_gb = 'KY369913.1'

    #  here we take the blast results and save both the fasta and the gbk file for pulling of annotations
    cmd = '/Users/uwvirongs/edirect/esearch -db nucleotide -query ' + ref_seq_gb + ' | /Users/uwvirongs/edirect/efetch -format fasta > ' + strain + '/' + strain + '_ref.fasta'
    ps = subprocess.Popen(cmd, stdin=subprocess.PIPE, shell=True)
    ps.communicate()

    cmd = '/Users/uwvirongs/edirect/esearch -db nucleotide -query ' + ref_seq_gb + ' | /Users/uwvirongs/edirect/efetch -format gb > ' + strain + '/' + strain + '_ref.gbk'

    ds = subprocess.Popen(cmd, stdin=subprocess.PIPE, shell=True)
    ds.communicate()

    subprocess.call('cat ' + strain + '/' + strain + '_ref.fasta ' + strain + '/' + strain + '.fasta > ' + strain + '/' + strain + '.aligner', shell=True)
    subprocess.call('mafft --auto ' + strain + '/' + strain + '.aligner > ' + strain + '/' + strain + '.ali', shell=True)
    # returns name of virus - and has saved .ali file as well as a .gbk file

    # simply splits the aligned data into two different strings
    ali_list, ali_genomes = read_fasta(strain + '/' + strain + '.ali')
    ref_seq = ali_genomes[0]
    our_seq = ali_genomes[1]
    return name_of_virus, our_seq, ref_seq


# Takes the name of one of the virus samples and returns a list of all the metadata we have in the sheet
def pull_metadata(virus_name, meta_data_info_sheet_loc):
    for line in open(meta_data_info_sheet_loc):
        if line.split(',')[1] == virus_name:
            # Just returns all of the data in a list (this is a csv so split(',') just pulls it all
            return line.split(',')


# Takes in two alinged sequences, searches for a specific place in the unaligned reference sequence - marks an index
# Then goes to that exact index in our sequence and counts the relative location for our own sequence - should annotate
# correctly
def find_loc_in_our_seq_from_reference(start, our_seq, reference_seq):
    ref_index = 0
    # this converts our index into 0 based
    count = start
    while count > 0:
        nt = reference_seq[ref_index]
        if nt != '-':
            count -= 1
        ref_index += 1

    our_index = 0
    for x in range(0, ref_index):
        nt = our_seq[x]
        if nt != '-':
            our_index += 1

    return our_index


def build_num_arrays(our_seq, ref_seq):
    ref_count = 0
    our_count = 0
    ref_num_array = []
    our_num_array = []

    for x in range(0, len(ref_seq)):
        if ref_seq[x] != '-':
            ref_count += 1
            ref_num_array.append(ref_count)
        else:
            ref_num_array.append(-1)

        if our_seq[x] != '-':
            our_count += 1
            our_num_array.append(our_count)
        else:
            our_num_array.append(-1)

    return our_num_array, ref_num_array


def adjust(given_num, our_num_array, ref_num_array):

    # Go through our number array and search for the number of interest
    for x in range(0, len(our_num_array)):
        if ref_num_array[x] == given_num:
            print('what should be equal to ' + str(given_num) + ' is ' + str(ref_num_array[x]))
            print('and this is the index ' + str(x))
            print('and array at position ' + str(x) + ' is ' + str(our_num_array[x]))
            index = x
            break

    # now index is the absolute location of what we want
    return str(our_num_array[index])


def pull_correct_annotations(strain, our_seq, ref_seq):
    # Read the reference gff file and extract lists of all of the protein locations and annotations!
    gene_loc_list = []
    gene_product_list = []
    print(our_seq)
    print(ref_seq)
    for line in open(strain + '/' + strain + '_ref.gbk'):
        if ' CDS ' in line:
            #print(line)
            # this is now going to be a list of numbers, start-stop start-stop
            gene_loc_list.append(re.findall(r'\d+', line))
            allow_one = True
        if '/product="' in line and allow_one:
            allow_one = False
            #print(line)
            gene_product_list.append(line.split('=')[1][1:-2])
    #print('gene locs before adjustment')
    #print(gene_loc_list)
    our_seq_num_array, ref_seq_num_array = build_num_arrays(our_seq, ref_seq)
    print(our_seq_num_array)
    print(ref_seq_num_array)
    for x in range(0, len(gene_loc_list)):
        for y in range(0, len(gene_loc_list[x])):
            gene_loc_list[x][y] = adjust(int(gene_loc_list[x][y]), our_seq_num_array, ref_seq_num_array)
    #print('products and numbers after adjustment')
    #print(gene_product_list)
    #print(gene_loc_list)
    return gene_loc_list, gene_product_list


def write_fasta(strain, genome):
    w = open(strain + '/' + strain + '.fasta', 'w')
    w.write('>' + strain + '\n')
    w.write(genome)
    w.close()


# Grab the coverage data from the imported metadata, if there is no coverage information randomly generate it and add
# .1x so that we know that it was randomly generated
def pull_coverage(data_list):
    # This protects us from when the metadata list doesn't contain any info on the virus
    if data_list is not None and data_list[20] != '':
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


def write_tbl(strain, gene_product_list, gene_locations, genome):
    #print('now we in write_tbl')
    #print('gene products')
    #print(gene_product_list)
    #print('gene locations')
    #print(gene_locations)
    tbl = open(strain + '/' + strain + '.tbl', 'w')
    tbl.write('>Feature ' + strain)
    flag = ''
    # TODO: add the ribosomal slippage line here
    for x in range(0, len(gene_product_list)):
        product = gene_product_list[x]
        location_info = gene_locations[x]
        if len(location_info) == 4:
            start_1 = str(location_info[0])
            end_1 = str(location_info[1])
            start_2 = str(location_info[2])
            end_2 = str(location_info[3])
            tbl.write('\n' + start_1 + '\t' + end_1 + '\tCDS\n')
            tbl.write(start_2 + '\t' + end_2 + '\n')
            tbl.write('\t\t\tproduct\t' + product + '\n')
            tbl.write('\t\t\texception\tRibosomal Slippage\n')
        else:
            start = location_info[0]
            end = location_info[1]
            if end == len(genome):
                if ((end - start) + 1 % 3) == 0 and genome[end - 3:end].upper() in 'TGA,TAA,TAG,UGA,UAA,UAG':
                    flag = ''
                else:
                    flag = '>'
            tbl.write('\n' + str(start) + '\t' + flag + str(end) + '\tCDS\n')
            tbl.write('\t\t\tproduct\t' + product)
    tbl.close()

def rna_editing():
    x =0

def annotate_a_virus(strain, genome, metadata_location, sbt_loc):
    subprocess.call('mkdir -p ' + strain, shell=True)

    write_fasta(strain, genome)

    name_of_virus, our_seq, ref_seq = blast_n_shit(strain, strain + '/' + strain + '.fasta')

    gene_loc_list, gene_product_list = pull_correct_annotations(strain, our_seq, ref_seq)

    metadata_list = pull_metadata(strain, metadata_location)

    coverage = pull_coverage(metadata_list)

    write_cmt(strain, coverage)

    subprocess.call('cp ' + sbt_loc + ' ' + strain + '/', shell=True)

    write_fsa(strain, name_of_virus, genome)

    write_tbl(strain, gene_product_list, gene_loc_list, our_seq)

    if 'parainfluenza virus' in name_of_virus.lower():
        if '3' in name_of_virus:
            x = 0
            # go into phosphoprotein - find poly G track and 2 g's and +2 to stop codon and then manually translate it
            # Save translation so that we can overwrite the .sqn file
        elif '4b' in name_of_virus:
            x = 0
            # Same as 3 but do it with the D protein instead

    subprocess.call('tbl2asn -p ' + strain + '/ -t ' + strain + '/' + sbt_loc.split('/')[-1] +
                    ' -Y ' + strain + '/assembly.cmt -V vb', shell=True)


def write_fsa(strain, name_of_virus, virus_genome):
    fsa = open(strain + '/' + strain + '.fsa', 'w')
    fsa.write('>' + strain + ' [organism=' + name_of_virus + '] [collection-date=2016] [country=USA] '
              '[moltype=genomic] [host=Human] [gcode=1] [molecule=RNA] [strain=' + strain + ']\n')
    fsa.write(virus_genome)
    fsa.close()


if __name__ == '__main__':
    fasta_loc = 'HarvestSeqJ1-J2-J4-J5.fasta'
    start = timeit.default_timer()

    # TODO: add a flag for redoing tbl2asn that ONLY does that - i.e. you could manually edit files then crank em out
#    parser = argparse.ArgumentParser(description='Package a set of UW clinical virus sequences for submission, pulling '
#                                                 'virus name information from blast and annotations are contained '
#                                                 'inside the .gff file passed to the script originally')
 #   parser.add_argument('gff_file', help='Input file in .gff format, should contain annotations and genome information '
 #                                        'for each virus to be processed')
 #   parser.add_argument('metadata_info_sheet', help='The metadata sheet that contains whatever we have on these samples')
 #   parser.add_argument('sbt_file_loc', help='File path for the .sbt file that should contain author names mainly')
#
 #   args = parser.parse_args()

    virus_strain_list, virus_genome_list = read_fasta(fasta_loc)

    # TODO: change the name of the arguments
    # TODO: go look at my virus_annotate_code for how to read in fastas
    # TODO: then go ahead and change this to a looping generalizable program that'll actually function
    # TODO: test
    for x in range(0, len(virus_strain_list)):
        annotate_a_virus(virus_strain_list[x], virus_genome_list[x], 'UWVIROCLINSEQ.csv', 'template.sbt')

