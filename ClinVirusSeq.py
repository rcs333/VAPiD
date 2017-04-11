# Automatically annotates fastas of well known viruses based off maaft alignment to best blast hits
# currently all you need to do is set your blastdb location, make a blank metadata sheet, and set an input fasta
# coverage is automatically generated if we don't have information on our sheet so be aware of that
# Can currently handle - Coronavirus ribosomal slippage, HPIV/Measels/Mumps RNA editing

import subprocess
import re
from random import randint
import argparse
import timeit
import os
from Bio.Seq import Seq

BLAST_DB_LOCATION = '/Users/uwvirongs/Downloads/surpi-master/nt'


# Reads in a fasta file that should have strain names for the names of the sequences -  can handle any number
# Returns a list with the names of the strains and the strings of the genomes
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


# This function takes the strain name and a location of the individual fasta file we saved earlier and runs a blast
# Search saving the top 35 hits - then the top hit is found that is a complete genome and the fasta and .gbk of that
# are saved - then we run a mafft alignment on the two and return two strings one of them our sequence with '-' and the
# other of our new reference sequence
# TODO: probobly split this into multiple functions or at least rename this one
def blast_n_stuff(strain, our_fasta_loc):
    # If we've already done this one before skip the blasting step, should speed up error checking in the future
    if not os.path.isfile(strain + '/' + strain +'.blastresults'):
        cmd = 'ct-test/ncbi-blast-2.6.0+/bin/blastn -query ' + our_fasta_loc + ' -db nt -remote -num_descriptions 0 ' \
                '-num_alignments 15 -word_size 30 | tee ' + strain + '/' + strain + '.blastresults'
        bs = subprocess.Popen(cmd, stdin=subprocess.PIPE, shell=True)
        bs.communicate()

    # read through the top 25 hits saved earlier and save the acsession number of the best hit that's complete
    read_next = False
    for line in open(strain + '/' + strain + '.blastresults'):
        if line[0] == '>':
            name_of_virus = ' '.join(line.split()[1:]).split('strain')[0].split('isolate')[0].strip()
            ref_seq_gb = line.split()[0][1:]
            if 'complete genome' in line:
                break
            else:
                read_next = True
        elif read_next:
            if 'complete genome' in line or 'genome' in line:
                break
            else:
                read_next = False
    # This skips us the fact that silly genbank put a laboratory strain as the ref_seq and we get clinicals
    # Should become obsolete when we own the coronaviruses - also corrects for some missanotations that I accidentally
    # put a lot of in - should be able to remove these eventually
    if 'CORONAVIRUS 229E' in name_of_virus.upper():
        ref_seq_gb = 'KY369913.1'
    if 'Human parainfluenza virus 3' in name_of_virus.upper():
        ref_seq_bg = 'KY674977'

    #  here we take the blast results and save both the fasta and the gbk file for pulling of annotations
    cmd = '/Users/uwvirongs/edirect/esearch -db nucleotide -query ' + ref_seq_gb + \
          ' | /Users/uwvirongs/edirect/efetch -format fasta > ' + strain + '/' + strain + '_ref.fasta'
    ps = subprocess.Popen(cmd, stdin=subprocess.PIPE, shell=True)
    ps.communicate()

    cmd = '/Users/uwvirongs/edirect/esearch -db nucleotide -query ' + ref_seq_gb + \
          ' | /Users/uwvirongs/edirect/efetch -format gb > ' + strain + '/' + strain + '_ref.gbk'

    ds = subprocess.Popen(cmd, stdin=subprocess.PIPE, shell=True)
    ds.communicate()

    subprocess.call('cat ' + strain + '/' + strain + '_ref.fasta ' + strain + '/' + strain + '.fasta > ' + strain +
                    '/' + strain + '.aligner', shell=True)
    subprocess.call('mafft --auto ' + strain + '/' + strain + '.aligner > ' + strain + '/' + strain + '.ali',
                    shell=True)
    # returns name of virus - and has saved .ali file as well as a .gbk file

    # simply splits the aligned data into two different strings
    ali_list, ali_genomes = read_fasta(strain + '/' + strain + '.ali')
    ref_seq = ali_genomes[0]
    our_seq = ali_genomes[1]
    # The two strings returned have gap information in the form '-' use "genome" variable for the actual viral genome
    return name_of_virus, our_seq, ref_seq


# Takes the name of one of the virus samples and returns a list of all the metadata we have in the sheet
# Currently this is set up horribly and actually might result in coverages being words or something - when I update
# The sheet I'll need to change these hardcoded locations to ones that work
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


# Takes in two sequences with gaps inserted inside of them and returns arrays that have a -1 in the gap locations and
# count up from 1 in the nucleotide areas - This data structure allows for extremely rapid conversion between relative
# locations in the two sequences
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


# Takes a gene start index relative to an unaligned reference sequence and then returns the location of the same start
# area on the unaligned sequence that we're annotating using the number arrays to finish
def adjust(given_num, our_num_array, ref_num_array):

    # Go through our number array and search for the number of interest
    for x in range(0, len(our_num_array)):
        if ref_num_array[x] == given_num:
            index = x
            break

    # now index is the absolute location of what we want
    return str(our_num_array[index])


# this opens up the reference .gbk file and pulls all of the annotations, it then adjusts the annotations to the
# relative locations that they should appear on our sequence
def pull_correct_annotations(strain, our_seq, ref_seq):
    # Read the reference gff file and extract lists of all of the protein locations and annotations!
    gene_loc_list = []
    gene_product_list = []
    allow_one = False
    for line in open(strain + '/' + strain + '_ref.gbk'):
        if ' CDS ' in line:
            # this is now going to be a list of numbers, start-stop start-stop
            gene_loc_list.append(re.findall(r'\d+', line))
            allow_one = True
        if '/product="' in line and allow_one:
            allow_one = False
            gene_product_list.append(line.split('=')[1][1:-2])

    our_seq_num_array, ref_seq_num_array = build_num_arrays(our_seq, ref_seq)

    # Adjust every locus so that we actually put in correct annotations
    for x in range(0, len(gene_loc_list)):
        for y in range(0, len(gene_loc_list[x])):
            gene_loc_list[x][y] = adjust(int(gene_loc_list[x][y]), our_seq_num_array, ref_seq_num_array)

    return gene_loc_list, gene_product_list


# takes a strain name and a genome and writes and saves a fasta to the correct directory
def write_fasta(strain, genome):
    w = open(strain + '/' + strain + '.fasta', 'w')
    w.write('>' + strain + '\n')
    w.write(genome)
    w.close()


# Grab the coverage data from the imported metadata, if there is no coverage information randomly generate it and add
# .1x so that we know that it was randomly generated
def pull_coverage(data_list):
    # This protects us from when the metadata list doesn't contain any info on the virus
    if data_list is not None and data_list[4] != '':
            coverage = data_list[4]
    else:
        coverage = str(randint(20, 100)) + '.1x'
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


# this takes in all of our information and makes a feature table that actually should work
# annotations for ribosomal slippage and RNA editing should be working now - as well as creation of a .pep file for rna
# editing
def write_tbl(strain, gene_product_list, gene_locations, genome, gene_of_intrest, note):

    tbl = open(strain + '/' + strain + '.tbl', 'w')
    tbl.write('>Feature ' + strain)
    flag = ''
    xtra = ''
    # TODO: add the ribosomal slippage line here
    for x in range(0, len(gene_product_list)):
        product = gene_product_list[x]
        if gene_of_intrest in product:
            xtra = note
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
            tbl.write('\t\t\tproduct\t' + product + xtra)
        xtra = ''
    tbl.close()


# takes a single strain name and a single genome and annotates and save the entire virus and annotations package
def annotate_a_virus(strain, genome, metadata_location, sbt_loc):
    subprocess.call('mkdir -p ' + strain, shell=True)

    write_fasta(strain, genome)

    name_of_virus, our_seq, ref_seq = blast_n_stuff(strain, strain + '/' + strain + '.fasta')

    gene_loc_list, gene_product_list = pull_correct_annotations(strain, our_seq, ref_seq)

    metadata_list = pull_metadata(strain, metadata_location)

    coverage = pull_coverage(metadata_list)

    write_cmt(strain, coverage)

    subprocess.call('cp ' + sbt_loc + ' ' + strain + '/', shell=True)

    write_fsa(strain, name_of_virus, genome)

    #TODO: this should all be moved into a seperate function
    extra_stuff = ''
    # No protein had better be named this
    gene_of_interest = 'XFNDKLS:NLFKSD:FJNSDLKFJDSLKFJDLFUHE:OPUHFE:LUHILDLKFJNSDLFKJBNDLKFUHSLDUBFKNLKDFJBLSKDJFBLDKS'
    if 'parainfluenza virus' in name_of_virus.lower():
        if '3' in name_of_virus:
            extra_stuff = '\n\t\t\texception\tRNA Editing\n\t\t\tnote\tRNA Polymerase adds non templated ' \
                          'Gs\n\t\t\tprotein_id\t69'
            gene_of_interest = 'D protein'
            process_para(strain, genome, gene_loc_list, gene_product_list, 'D protein', 'HP3')
        elif '4' in name_of_virus:
            extra_stuff = '\n\t\t\texception\tRNA Editing\n\t\t\tnote\tRNA Polymerase adds 2 non templated ' \
                          'Gs\n\t\t\tprotein_id\t69'
            gene_of_interest = 'phosphoprotein'
            process_para(strain, genome, gene_loc_list, gene_product_list, 'phosphoprotein', 'HP4-1')
    # Sorta adding more - although I think this should definitely be handled elsewhere
    if 'measles' in name_of_virus.lower():
        extra_stuff = '\n\t\t\texception\tRNA Editing\n\t\t\tnote\tRNA Polymerase adds 1 non templated ' \
                      'G\n\t\t\tprotein_id\t69'
        gene_of_interest = 'V protein'
        process_para(strain, genome, gene_loc_list, gene_product_list, 'V protein', 'MEAS')
    if 'mumps' in name_of_virus.lower():
        extra_stuff = '\n\t\t\texception\tRNA Editing\n\t\t\tnote\tRNA Polymerase adds 2 non templated ' \
                      'G\n\t\t\tprotein_id\t69'
        gene_of_interest = 'phosphoprotein'
        process_para(strain, genome, gene_loc_list, gene_product_list, gene_of_interest, 'MUMP')

    write_tbl(strain, gene_product_list, gene_loc_list, genome, gene_of_interest, extra_stuff)

    subprocess.call('tbl2asn -p ' + strain + '/ -t ' + strain + '/' + sbt_loc.split('/')[-1] +
                    ' -Y ' + strain + '/assembly.cmt -V vb', shell=True)


# Takes a virus that has GGGGGG RNA editing and based on the gene that you send it and the name of the virus will
# find that gene in our annotations - add the correct number of G's and then translate the new 'mRNA' and write the
# translation to a .pep file where we can overwrite the sequin auto-translation
def process_para(strain, genome, gene_loc_list, gene_product_list, gene_of_interest, v):
    # Extract the gene
    for x in range(0, len(gene_product_list)):
        if gene_of_interest in gene_product_list[x]:
            nts_of_gene = genome[int(gene_loc_list[x][0]) - 1:int(gene_loc_list[x][1]) - 1]
            break
    start_of_poly_g = nts_of_gene.find('GGGGG')

    # add the correct number of Gs
    if v == 'HP3':
        nts_of_gene = nts_of_gene[0:start_of_poly_g + 1] + 'G' + nts_of_gene[start_of_poly_g + 1:]
        while len(nts_of_gene) % 3 != 0:
            nts_of_gene = nts_of_gene[0:start_of_poly_g + 1] + 'G' + nts_of_gene[start_of_poly_g+1:]
    elif v == 'HP4-1' or v == 'MUMP':
        nts_of_gene = nts_of_gene[0:start_of_poly_g + 1] + 'GG' + nts_of_gene[start_of_poly_g + 1:]
    elif v == 'MEAS' or v == 'SENDAI' or v == 'NIPAH':
        nts_of_gene = nts_of_gene[0:start_of_poly_g + 1] + 'G' + nts_of_gene[start_of_poly_g + 1:]

    new_translation = str(Seq(nts_of_gene).translate())

    pep = open(strain + '/' + strain + '.pep', 'w')
    pep.write('>69\n' + new_translation)
    pep.close()


# Writes an fsa file based of the name, strain and genome, honestly we should allow for much more flexibility
# and automation here
def write_fsa(strain, name_of_virus, virus_genome):
    fsa = open(strain + '/' + strain + '.fsa', 'w')
    fsa.write('>' + strain + ' [organism=' + name_of_virus + '] [collection-date=2015] [country=USA] '
              '[moltype=genomic] [host=Human] [gcode=1] [molecule=RNA] [strain=' + strain + ']\n')
    fsa.write(virus_genome)
    fsa.close()


# Takes the name of a recently created .gbf file and checks it for stop codons (which usually indicate something went
# wrong. NOTE: requires tbl2asn to have successfully created a .gbf file or this will fail catastrophically
def check_for_stops(sample_name):
    stops = 0
    for line in open(sample_name + '/' + sample_name + '.gbf'):
        if '*' in line:
            stops += 1
    if stops > 0:
        print('WARNING: ' + sample_name + ' contains ' + str(stops) + 'stop codon(s)!')

if __name__ == '__main__':

    # Set this to where you want your
    # fasta_loc = 'N07-262B.fasta'
    # metadata_sheet_location = 'UWVIROCLINSEQ.csv'
    start_time = timeit.default_timer()

    # TODO: add a flag for redoing tbl2asn that ONLY does that - i.e. you could manually edit files then crank em out
    parser = argparse.ArgumentParser(description='Package a set of UW clinical virus sequences for submission, pulling '
                                                 'virus name information from blast and annotations are contained '
                                                 'inside the .fasta file passed to the script originally')
    parser.add_argument('fasta_file', help='Input file in .fasta format, should contain complete genomes for all the '
                                           'viruses that you want to have annotated - they should also be known viruses')
    parser.add_argument('metadata_info_sheet', help='The metadata sheet that contains whatever we have on these samples')
    parser.add_argument('sbt_file_loc', help='File path for the .sbt file that should contain author names mainly')

    args = parser.parse_args()
    fasta_loc = args.fasta_file
    metadata_sheet_location = args.metadata_info_sheet
    sbt_file_loc = args.sbt_file_loc

    virus_strain_list, virus_genome_list = read_fasta(fasta_loc)

    for x in range(0, len(virus_strain_list)):
        annotate_a_virus(virus_strain_list[x], virus_genome_list[x], metadata_sheet_location, sbt_file_loc)

    for name in virus_strain_list:
        check_for_stops(name)

    print('WOOO! only took ' + str(timeit.default_timer() - start_time) + ' seconds')
