# Automatically annotates fastas of well known viruses based off maaft alignment to best blast hits
# currently all you need to do is set your blastdb location, make a blank metadata sheet, and set an input fasta
# Can currently handle - Coronavirus ribosomal slippage, HPIV/Measels/Mumps/Sendai RNA editing

import subprocess
import re
import argparse
import timeit
import os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Blast import NCBIWWW
import datetime
import platform
import shutil
import sys
from Bio import Entrez
Entrez.email = 'uwvirongs@gmail.com'


VERSION = 'v0.9'


# Reads in a fasta file that should have strain names for the names of the sequences -  can handle any number
# Returns a list with the names of the strains and the strings of the genomes
def read_fasta(fasta_file_loc):
    strain_list = []
    genome_list = []
    dna_string = ''
    for line in open(fasta_file_loc):
        if line[0] == '>':
            strain_list.append(line[1:].split()[0])
            if dna_string != '':
                genome_list.append(dna_string)
                dna_string = ''
        else:
            dna_string += line.strip()
    print(dna_string)
    x = 0
    while dna_string[x] == 'N' or dna_string[x] == '?':
        x += 1

    y = len(dna_string - 1)
    while dna_string[y] == 'N' or dna_string[y] == '?':
        y -= 1
    dna_string = dna_string[x:y]
    genome_list.append(dna_string)
    return strain_list, genome_list


# This function takes the strain name and a location of the individual fasta file we saved earlier and runs a blast
# Search saving the top 35 hits - then the top hit is found that is a complete genome and the fasta and .gbk of that
# are saved - then we run a mafft alignment on the two and return two strings one of them our sequence with '-' and the
# other of our new reference sequence
def blast_n_stuff(strain, our_fasta_loc):
    # If we've already done this one before skip the blasting step, should speed up error checking in the future
    if not os.path.isfile(strain + SLASH + strain + '.blastresults'):
        #cmd = 'blastn -query ' + our_fasta_loc + ' -db nt -remote -num_descriptions 0 ' \
        #        '-num_alignments 35 -word_size 28 | tee ' + strain + '/' + strain + '.blastresults'
        #bs = subprocess.Popen(cmd, stdin=subprocess.PIPE, shell=True)
        #bs.communicate()
        # This should cross platform the results and just save the blast results straight up in the same format as before
        record = open(our_fasta_loc).read()
        result_handle = NCBIWWW.qblast('blastn', 'nt', record, word_size=28, descriptions=0, alignments=35, format_type='Text')
        with open(strain + SLASH + strain + '.blastresults', 'w') as out_handle:
            out_handle.write(result_handle.read())
        result_handle.close()

    # read through the top hits saved earlier and save the accession number of the best hit that's complete
    read_next = False
    found = False
    for line in open(strain + SLASH + strain + '.blastresults'):
        if line[0] == '>':
            name_of_virus = ' '.join(line.split()[1:]).split('strain')[0].split('isolate')[0].strip()
            ref_seq_gb = line.split()[0][1:]

            # last part of these two logic checks is so we avoid the misassembled/mutated viruses
            # This is going to get really out of hand if we have to keep blacklisting records
            if 'complete genome' in line and ref_seq_gb.split('.')[0] not in 'KM551753 GQ153651':
                found = True
                break
            else:
                read_next = True
        elif read_next:
            if 'complete genome' in line or 'genome' in line and ref_seq_gb.split('.')[0] not in 'KM551753 GQ153651':
                found = True
                break
            else:
                read_next = False

    if not found:
        for line in open(strain + SLASH + strain + '.blastresults'):
            if line[0] == '>':
                name_of_virus = ' '.join(line.split()[1:]).split('strain')[0].split('isolate')[0].strip()
                ref_seq_gb = line.split()[0][1:]
                break


    # TODO: put a function in here that error checks and does blacklisting conversions as well as renaming
    # This skips us the fact that silly genbank put a laboratory strain as the ref_seq and we get clinicals
    # Should become obsolete when we own the coronaviruses - also corrects for some missanotations that I accidentally
    # put a lot of in - should be able to remove these eventually
    if args.r:
        ref_seq_gb = args.r
    if 'CORONAVIRUS 229E' in name_of_virus.upper():
        ref_seq_gb = 'KY369913.1'
    #if 'HUMAN PARAINFLUENZA VIRUS 3' in name_of_virus.upper():
    #    ref_seq_gb = 'KY674977'
    if 'HUMAN RESPIROVIRUS 3' in name_of_virus.upper():
        ref_seq_gb = 'KY369864'
    if 'HUMAN IMMUNODEFICIENCY VIRUS TYPE 1' in name_of_virus.upper():
        ref_seq_gb = 'L20587.1'
    if 'MEASLES' in name_of_virus.upper():
        ref_seq_gb = 'EU293548'
    #  here we take the blast results and save both the fasta and the gbk file for pulling of annotations
    #
    #cmd = 'esearch -db nucleotide -query ' + ref_seq_gb + \
    #      ' | efetch -format fasta > ' + strain + '/' + strain + '_ref.fasta'
    #ps = subprocess.Popen(cmd, stdin=subprocess.PIPE, shell=True)
    #ps.communicate()

    # This replacement code should work I think
    record = Entrez.read(Entrez.esearch(db='nucleotide', term=ref_seq_gb))
    h = Entrez.efetch(db='nucleotide', id=record["IdList"][0], rettype='fasta', retmode='text')
    d = open(strain + SLASH + strain + '_ref.fasta', 'w')
    d.write(h.read())
    d.close()

    #cmd = 'esearch -db nucleotide -query ' + ref_seq_gb + \
    #      ' | efetch -format gb > ' + strain + '/' + strain + '_ref.gbk'

    #ds = subprocess.Popen(cmd, stdin=subprocess.PIPE, shell=True)
    #ds.communicate()

    h2 = Entrez.efetch(db='nucleotide', id=record["IdList"][0], rettype='gb', retmode='text')
    e = open(strain + SLASH + strain + '_ref.gbk', 'w')
    e.write(h2.read())
    e.close()

    # cat x-platform code lol
    #subprocess.call('cat ' + strain + '/' + strain + '_ref.fasta ' + strain + '/' + strain + '.fasta > ' + strain +
    #                '/' + strain + '.aligner', shell=True)
    # this code should be completely unecessary
    g = open(strain + SLASH + strain + '.aligner', 'w')
    t1 = open(strain + SLASH + strain + '_ref.fasta')
    t2 = open(strain + SLASH + strain + '.fasta')
    for line in t1:
        g.write(line)
    for line in t2:
        g.write(line)
    g.close()


    #subprocess.call('mafft --auto ' + strain + '/' + strain + '.aligner > ' + strain + '/' + strain + '.ali',
    #                shell=True)
    # returns name of virus - and has saved .ali file as well as a .gbk file

    # simply splits the aligned data into two different strings
    #ali_list, ali_genomes = read_fasta(strain + SLASH + strain + '.ali')
    #ref_seq = ali_genomes[0]
    #our_seq = ali_genomes[1]
    # The two strings returned have gap information in the form '-' use "genome" variable for the actual viral genome
    from Bio import pairwise2
    seq1 = SeqIO.read(strain + SLASH + strain + '_ref.fasta', 'fasta')
    seq2 = SeqIO.read(strain + SLASH + strain + '.fasta', 'fasta')
    alignments = pairwise2.align.globalxx(seq1.seq, seq2.seq)
    ref_seq = alignments[0][0]
    our_seq = alignments[0][1]
    return name_of_virus, our_seq, ref_seq


# Takes the name of one of the virus samples and returns a list of all the metadata we have in the sheet
# Currently this is set up horribly and actually might result in coverages being words or something - when I update
# The sheet I'll need to change these hardcoded locations to ones that work
def pull_metadata(virus_name, meta_data_info_sheet_loc):
    for line in open(meta_data_info_sheet_loc):
        if line.split(',')[1] == virus_name:
            # Just returns all of the data in a list (this is a csv so split(',') just pulls it all
            return line.split(',')


# Takes in two aligned sequences, searches for a specific place in the unaligned reference sequence - marks an index
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
# locations in the two sequences although does assume that these genes are of uniform length
# NOTE: This means that when we have reads that like don't have the start codons of the first gene or something we'll
# get a -1 for the start location on our annotation
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
    for line in open(strain + SLASH + strain + '_ref.gbk'):
        if ' CDS ' in line:
            # this is now going to be a list of numbers, start-stop start-stop
            gene_loc_list.append(re.findall(r'\d+', line))
            allow_one = True
        if '/product="' in line and allow_one:
            allow_one = False
            # Inconsistent naming of protein products
            px = line.split('=')[1][1:-2]
            if px == 'phospho protein':
                px = 'phoshoprotein'
            gene_product_list.append(px)
    our_seq_num_array, ref_seq_num_array = build_num_arrays(our_seq, ref_seq)


    # Adjust every locus so that we actually put in correct annotations
    for entry in range(0, len(gene_loc_list)):
        for y in range(0, len(gene_loc_list[entry])):
            gene_loc_list[entry][y] = adjust(int(gene_loc_list[entry][y]), our_seq_num_array, ref_seq_num_array)

    return gene_loc_list, gene_product_list


# takes a strain name and a genome and writes and saves a fasta to the correct directory
def write_fasta(strain, genome):
    w = open(strain + SLASH + strain + '.fasta', 'w')
    w.write('>' + strain + '\n')
    w.write(genome)
    w.close()


# Grab the coverage data from the imported metadata, if there's no information coverage is set to zero and axed from the
# .cmt file in the write_cmt() function
def pull_coverage(data_list):
    # This protects us from when the metadata list doesn't contain any info on the virus or there's no coverage
    # also means you can just pass a blank csv file if you don't have coverages
    coverage = ''
    if data_list is not None and data_list[4] != '':
        coverage = data_list[4]

    return coverage


# grab the collection date from the imported metadata - if nothing is put in the metadata sheet collection date will not
# be included in the meta data
def pull_col_date(data_list):
    col_date = ''
    if data_list is not None and data_list[5] != '':
        # formatted a bit differently because this goes in the fsa file header not the comment file
        col_date = ' [collection-date=' + data_list[5] + ']'
    return col_date


# Take the name of a virus sample, and write the .cmt file for it using supplied coverage information
# This is something that only we'll know for sure and don't want to reuire it, - so I'll keep it untill release
# TODO: remove the assembly data and method and technology before release
def write_cmt(sample_name, coverage):
    cmt = open(sample_name + SLASH +'assembly.cmt', 'w')
    cmt.write('##Assembly-Data-START##\n')
    cmt.write('Assembly Method\tGeneious v. 9.1\n')
    if coverage != '':
        cmt.write('Coverage\t' + coverage + '\n')
    cmt.write('Sequencing Technology\tIllumina\n')
    cmt.write('##Assembly-Data-END##\n')
    cmt.close()


# this takes in all of our information and makes a feature table that actually should work
# annotations for ribosomal slippage and RNA editing should be working now - as well as creation of a .pep file for rna
# editing
def write_tbl(strain, gene_product_list, gene_locations, genome, gene_of_intrest, note):

    tbl = open(strain + SLASH + strain + '.tbl', 'w')
    tbl.write('>Feature ' + strain)

    for x in range(0, len(gene_product_list)):
        flag = ''
        xtra = ''
        sflag = ''
        product = gene_product_list[x]
        if gene_of_intrest in product:
            xtra = note
        #TODO: this is inelegant, need to make it HIV specific i think
        if 'Pol polyprotein' == product or 'Pol' == product:
            sflag = '<'
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
            start = int(location_info[0])
            end = int(location_info[1])
            if end >= len(genome) and genome[end - 3:end].upper() not in 'TGA,TAA,TAG,UGA,UAA,UAG':
                flag = '>'

            it_count = 0
            modifid_orf = False
            while genome[end - 3:end].upper() not in 'TGA,TAA,TAG,UGA,UAA,UAG' and end < len(genome) - 3 and it_count <= 3:
                print('Modifying ORF length for ' + str(product))
                modifid_orf = True
                end += 3
                it_count += 1
            if modifid_orf and genome[end - 3:end].upper() not in 'TGA,TAA,TAG,UGA,UAA,UAG':
                end -= it_count * 3

            # This should now correctly annotate assemblies that come in with the very edges chopped off
            if int(start) < 1:
                sflag = '<'
                start = (int(end) % 3) + 1

            if int(end) < 1:
                end = len(genome)
                flag = '>'
            tbl.write('\n' + sflag + str(start) + '\t' + flag + str(end) + '\tCDS\n')
            tbl.write('\t\t\tproduct\t' + product + xtra)
        # TODO: examine the effects of this line is it necessary - is it deleting important info? I thnk it does nothing
        xtra = ''
    tbl.write('\n')
    tbl.close()


# takes a single strain name and a single genome and annotates and save the entire virus and annotations package
# returns the "species" of the virus for consolidated .sqn packaging
def annotate_a_virus(strain, genome, metadata, coverage, sbt_loc):
    # subprocess.call('mkdir -p ' + strain, shell=True) below code should replicate mkdir -p
    if not os.path.exists(strain):
        os.makedirs(strain)

    write_fasta(strain, genome)

    name_of_virus, our_seq, ref_seq = blast_n_stuff(strain, strain + SLASH + strain + '.fasta')

    gene_loc_list, gene_product_list = pull_correct_annotations(strain, our_seq, ref_seq)

    #metadata_list = pull_metadata(strain, metadata_location)

    #coverage = pull_coverage(metadata_list)

    #col_date = pull_col_date(metadata_list)

    write_cmt(strain, coverage)

    # subprocess.call('cp ' + sbt_loc + ' ' + strain + '/', shell=True)
    # shutil.copyfile(sbt_loc, strain + SLASH + strain + '.sbt')

    write_fsa(strain, name_of_virus, genome, metadata)

    extra_stuff = ''
    # No protein had better be named this
    gene_of_interest = 'XFNDKLS:NLFKSD:FJNSDLKFJDSLKFJDLFUHE:OPUHFE:LUHILDLKFJNSDLFKJBNDLKFUHSLDUBFKNLKDFJBLSKDJFBLDKS'
    if 'respirovirus' in name_of_virus.lower() or 'parainfluenza virus 3' in name_of_virus.lower():
        if '3' in name_of_virus:
            extra_stuff = '\n\t\t\texception\tRNA Editing\n\t\t\tnote\tRNA Polymerase adds non templated ' \
                          'Gs\n\t\t\tprotein_id\tn_' + strain
            gene_of_interest = 'D protein'
            process_para(strain, genome, gene_loc_list, gene_product_list, 'D protein', 'HP3')
        elif '1' in name_of_virus:
            extra_stuff = 'WEGOTAPARA1'
            gene_of_interest ='C\' protein'
    if 'parainfluenza virus 4a' in name_of_virus.lower():
        extra_stuff = '\n\t\t\texception\tRNA Editing\n\t\t\tnote\tRNA Polymerase adds 2 non templated ' \
                      'G\n\t\t\tprotein_id\tn_' + strain
        gene_of_interest = 'phosphoprotein'
        process_para(strain, genome, gene_loc_list, gene_product_list, 'phoshoprotein', 'HPIV4a')
    # Sorta adding more - although I think this should definitely be handled elsewhere
    if 'measles' in name_of_virus.lower():
        extra_stuff = '\n\t\t\texception\tRNA Editing\n\t\t\tnote\tRNA Polymerase adds 1 non templated ' \
                      'G\n\t\t\tprotein_id\tn_' + strain
        gene_of_interest = 'V protein'
        process_para(strain, genome, gene_loc_list, gene_product_list, 'V protein', 'MEAS')
    if 'mumps' in name_of_virus.lower():
        extra_stuff = '\n\t\t\texception\tRNA Editing\n\t\t\tnote\tRNA Polymerase adds 2 non templated ' \
                      'G\n\t\t\tprotein_id\tn_' + strain
        gene_of_interest = 'phosphoprotein'
        process_para(strain, genome, gene_loc_list, gene_product_list, gene_of_interest, 'MUMP')
    if 'rubulavirus 4' in name_of_virus:
        extra_stuff = '\n\t\t\texception\tRNA Editing\n\t\t\tnote\tRNA Polymerase adds 2 non templated ' \
                      'Gs\n\t\t\tprotein_id\tn_' + strain
        gene_of_interest = 'phosphoprotein'
        process_para(strain, genome, gene_loc_list, gene_product_list, 'phoshoprotein', 'HP4-1')
    #if 'respirovirus' in name_of_virus.lower():
    #    extra_stuff = '\n\t\texception\tRNA Editing\n\t\t\t\tnote\tRNA Polymerase adds 2 non templated ' \
    #                  'G\n\t\t\tprotein_id\tn_' + strain
    #    gene_of_interest = 'D protein'
    #    process_para(strain, genome, gene_loc_list, gene_product_list, gene_of_interest, 'SENDAI')
    if 'metapneumovirus' in name_of_virus.lower():
        put_start = int(gene_loc_list[7][0])
        print('start of gene is ' + str(put_start))
        orf = genome[put_start - 1:put_start + 4000]
        print('orf length is ' + str(len(orf)))
        print('genome length is ' + str(len(genome)))
        orf_trans = str(Seq(orf).translate())
        print('orf translation is ' + orf_trans)
        put_end = (orf_trans.find('*') * 3)
        print('putative end is ' + str(put_end))
        gene_loc_list[7][1] = put_start + put_end

    write_tbl(strain, gene_product_list, gene_loc_list, genome, gene_of_interest, extra_stuff)
    #     if SLASH == '/':
    #     subprocess.call('tbl2asn -p ' + strain + '/ -t ' + strain + '/' + sbt_loc.split('/')[-1] +
    #                    ' -Y ' + strain + '/assembly.cmt -V vb', shell=True)
    # else:
    # replaced it all with an x-platform call

    cmd = 'tbl2asn -p ' + strain + SLASH + ' -t ' + sbt_loc + ' -Y ' + strain + SLASH + 'assembly.cmt -V vb'
    subprocess.call(cmd, shell=True)
    return name_of_virus


# this is now a dummy function for if I ever have to do c term scanning for stop codons again
# TODO: get some test fastas and get this working
def fix_orf(start, len):
    return start


# This function takes a nucleotide sequence that has non-templated G's inserted an unknown number of times and trims
# Them from the end to read into a correct frame  - translates the sequence and picks the one with the least number of
# stop codons returns the raw sequence
def pick_correct_frame(one, two):

    # we already added the G's to the start - and generally we hit the stop codon exactly where the annotation says we
    # will so this just avoids some weirdness with passing sequences of length not divisible by three to seq.translate()
    while len(one) % 3 != 0:
        one = one[:-1]
    while len(two) % 3 != 0:
        two = two[:-1]
    one_trans = str(Seq(one).translate())
    two_trans = str(Seq(two).translate())
    one_count = one_trans.count('*')
    two_count = two_trans.count('*')
    # Troubleshooting code that I'm keeping in for when we add more viruses that have non templated G's
    print('adding two Gs gives ' + str(two_count) + ' stop codon(s)')
    print(two_trans)
    print('adding one G gives ' + str(one_count) + ' stop codon(s)')
    print(one_trans)
    if one_count < two_count:
        print('chose one')
        return one
    else:
        print('chose two')
        return two


# Takes a virus that has GGGGGG RNA editing and based on the gene that you send it and the name of the virus will
# find that gene in our annotations - add the correct number of G's and then translate the new 'mRNA' and write the
# translation to a .pep file where we can overwrite the sequin auto-translation
def process_para(strain, genome, gene_loc_list, gene_product_list, gene_of_interest, v):
    # Extract the gene protected because everything we throw in here are guaranteed to have the gene of interest
    for g in range(0, len(gene_product_list)):
        print(gene_of_interest)
        print(gene_product_list[g])
        # flipping this covers whack spacing in protein products
        if gene_of_interest in gene_product_list[g]:
            nts_of_gene = genome[int(gene_loc_list[g][0]) - 1:int(gene_loc_list[g][1]) - 1]
            break

    # Since we now have a priori knowledge of where this RNA editing takes place we can safely search within a substring
    # start_of_poly_g = nts_of_gene.find('GGGGG')

    # add the correct number of Gs
    if v == 'HP3':
        start_of_poly_g = nts_of_gene.find('GGGGG', 700, 740)
        nts_of_gene_1 = nts_of_gene[0:start_of_poly_g + 1] + 'G' + nts_of_gene[start_of_poly_g + 1:]
        nts_of_gene_2 = nts_of_gene[0:start_of_poly_g + 1] + 'GG' + nts_of_gene[start_of_poly_g + 1:]
        nts_of_gene = pick_correct_frame(nts_of_gene_1, nts_of_gene_2)
    # despite viral zone saying SENDAI adds 1 G adding two removes stop codon's - remains tbd if variable
    elif v == 'SENDAI':
        start_of_poly_g = nts_of_gene.find('AAAAGGG')
        nts_of_gene_1 = nts_of_gene[0:start_of_poly_g + 1] + 'G' + nts_of_gene[start_of_poly_g + 1:]
        nts_of_gene_2 = nts_of_gene[0:start_of_poly_g + 1] + 'GG' + nts_of_gene[start_of_poly_g + 1:]
        nts_of_gene = pick_correct_frame(nts_of_gene_1, nts_of_gene_2)
    elif v == 'HP4-1':
        start_of_poly_g = nts_of_gene.find('AAGAGG', 435, 460)
        nts_of_gene = nts_of_gene[0:start_of_poly_g + 1] + 'GG' + nts_of_gene[start_of_poly_g + 1:]
    elif v == 'MUMP':
        start_of_poly_g = nts_of_gene.find('AAGAGG', 445, 465)
        nts_of_gene = nts_of_gene[0:start_of_poly_g + 1] + 'GG' + nts_of_gene[start_of_poly_g + 1:]
    elif v == 'MEAS':
        start_of_poly_g = nts_of_gene.find('AAAAAGG', 674, 695)
        nts_of_gene = nts_of_gene[0:start_of_poly_g + 1] + 'G' + nts_of_gene[start_of_poly_g + 1:]
    elif v == 'NIPAH':
        start_of_poly_g = nts_of_gene.find('AAAAAAGG', 705, 725)
        nts_of_gene = nts_of_gene[0:start_of_poly_g + 1] + 'G' + nts_of_gene[start_of_poly_g + 1:]
    elif v ==  'HPIV4a':
        start_of_poly_g = nts_of_gene.find('AAGAGG', 439, 460)
        nts_of_gene = nts_of_gene[0:start_of_poly_g + 1] + 'GG' + nts_of_gene[start_of_poly_g + 1:]


    new_translation = str(Seq(nts_of_gene).translate())

    pep = open(strain + SLASH + strain + '.pep', 'w')
    pep.write('>n_' + strain + '\n' + new_translation)
    pep.write('\n')
    pep.close()


# Writes an fsa file based of the name, strain and genome, honestly we should allow for much more flexibility
# and automation here
def write_fsa(strain, name_of_virus, virus_genome, col_date):
    fsa = open(strain + SLASH + strain + '.fsa', 'w')
    # new metadata input strat
    fsa.write('>' + strain + ' [organism=' + name_of_virus + ']' + '[moltype=genomic] [host=Human] [gcode=1] [molecule=RNA]' + metadata + '\n')
    fsa.write(virus_genome)
    fsa.write('\n')
    fsa.close()


# TODO: document how to use this functinoality, provide a example document from which to fill shit out and also write up a tutorial for using the program and actually producing something of value
def do_meta_data(strain, sheet_exists):
    meta_data = ''
    first = True
    s = ''
    coverage = ''
    if sheet_exists:
        for line in open(metadata_sheet_location):
            if first:
                names = line.split(',')
                first = False
            elif line.split(',')[0] == strain:
                for dex in range(0, len(names)):
                    if names[dex].strip() == 'coverage':
                        coverage = line.split(',')[dex].strip()
                    else:
                        s = s + ' [' + names[dex].strip() + '=' + line.split(',')[dex].strip() + ']'
                break

    if s == '':
        print('metadata not found in provided .csv or .csv not created -  time for manual entry for sequence - ' + strain)
        col = ' [collection-date=' + raw_input('Enter collection date in the format (23-Mar-2005, Mar-2005, or 2005): ').strip() + ']'
        con = ' [country=' + raw_input('Enter country sample was collected in (example - USA): ').strip() + ']'
        st = ' [strain=' + raw_input('Enter strain name - if unknown just put ' + strain + ': ').strip() + ']'
        cov = raw_input('Enter coverage as a number example 42.3, if unknown just leave this blank and hit enter: ')
        coverage = cov
        meta_data = col + con + st
    else:
        meta_data = s
    return meta_data, coverage


# Takes the name of a recently created .gbf file and checks it for stop codons (which usually indicate something went
# wrong. NOTE: requires tbl2asn to have successfully created a .gbf file or this will fail catastrophically
# Now also returns if stop codons are in it or not so they'll be omitted during the packaging phase
def check_for_stops(sample_name):
    stops = 0
    for line in open(sample_name + SLASH + sample_name + '.gbf'):
        if '*' in line:
            stops += 1
    if stops > 0:
        print('WARNING: ' + sample_name + ' contains ' + str(stops) + 'stop codon(s)!')
        return True
    else:
        return False


def check_os():
    if platform.system() == 'Linux' or platform.system() == 'Darwin':
        return '/'
    else:
        return '\\'


if __name__ == '__main__':

    start_time = timeit.default_timer()
    SLASH = check_os()
    parser = argparse.ArgumentParser(description='Version ' + VERSION + '\nPackage a set of UW clinical virus sequences'
                                                 ' for submission, pulling virus name information from blast and '
                                                 'annotations are contained inside the .fasta file passed to the script'
                                                 ' originally')
    parser.add_argument('fasta_file', help='Input file in .fasta format, should contain complete genomes for all the '
                                           'viruses that you want to have annotated - they should be known viruses')
    #parser.add_argument('metadata_info_sheet', help='The metadata sheet that contains whatever we have for the samples')
    parser.add_argument('author_template_file_loc', help='File path for the NCBI provided sequence author template file (should have a .sbt extension')
    parser.add_argument('--metadata_loc', help='If you\'ve input the metadata in the provided csv specify the location with this optional argument, otherwise all metadata will be manually prompted for')
    parser.add_argument('--r', help='If you want to specifiy a specific NCBI reference, post the accession number here - '
                                    'must be the exact accession number - note feature only works with one virus at a time ')
    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    fasta_loc = args.fasta_file

    sbt_file_loc = args.author_template_file_loc

    virus_strain_list, virus_genome_list = read_fasta(fasta_loc)

    strain2species = {}
    strain2stops = {}

    meta_list = []
    coverage_list = []
    if args.metadata_loc:
        metadata_sheet_location = args.metadata_loc
        for x in range(0, len(virus_strain_list)):
            metadata, coverage = do_meta_data(virus_strain_list[x], True)
            meta_list.append(metadata)
            coverage_list.append(coverage)
    else:
        for x in range(0, len(virus_strain_list)):
            metadata, coverage = do_meta_data(virus_strain_list[x], False)
            meta_list.append(metadata)
            coverage_list.append(coverage)

    for x in range(0, len(virus_strain_list)):
        strain2species[virus_strain_list[x]] = annotate_a_virus(virus_strain_list[x], virus_genome_list[x],
                                                                meta_list[x], coverage_list[x], sbt_file_loc)
        # now we've got a map of [strain] -> name of virus (with whitespace)

    for name in virus_strain_list:
        # now we've got a map of [strain] -> boolean value if there are stops or not
        strain2stops[name] = check_for_stops(name)

    # now we only consolidate the sequin files if the flag is passed, which will save time if I have to blast stuff
    # multiple times for troubleshooting
    # TODO: make this x-platform
    # if args.r:
    # =args.r
    #     strain2species_nostops = {}
    #     for item in strain2species.keys():
    #         if not strain2stops[item]:
    #             strain2species_nostops[item] = strain2species[item]
    #     # now we've got map of strain 'folder names' to virus names with no stops
    #     virus_species_list = []
    #     for item in strain2species_nostops.values():
    #         if '_'.join(item.split()) not in virus_species_list:
    #             virus_species_list.append('_'.join(item.split()))
    #     # now we've got a list of all the viruses in our fasta file
    #     now = datetime.datetime.now()
    #     date = now.strftime("%Y_%m-%d")
    #     # subprocess.call('mkdir -p ' + date, shell=True) Changed to os.
    #     os.makedirs(date)
    #     for item in virus_species_list:
    #         # subprocess.call('mkdir -p ' + date + '/' + item, shell=True)
    #         os.makedirs(date + SLASH + item)
    #     # now we've got all the folders for the viruses to go into
    #     for strain in strain2species_nostops.keys():
    #         # TODO: factor out this redundant loops and logic and data structures, however it *does* work
    #         species = '_'.join(strain2species_nostops[strain].split())
    #         # cmd = 'mv ' + strain + '/ ' + date + '/' + species
    #         # subprocess.call(cmd, shell=True)
    #         shutil.move(strain + SLASH, date + SLASH + species)
    #
    #
    #     # now we gotta extract the files and tbl2asn em'
    #     # This is absolutely disgusting code and I really really need to factor this out into it's own method
    #     for item in virus_species_list:
    #         subprocess.call('cat ' + date + '/' + item + '/*/*.tbl > ' + date + '/' + item + '/' + item + '.tbl',
    #                         shell=True)
    #         subprocess.call('cat ' + date + '/' + item + '/*/*.fsa > ' + date + '/' + item + '/' + item + '.fsa',
    #                         shell=True)
    #         subprocess.call('cat ' + date + '/' + item + '/*/*.pep > ' + date + '/' + item + '/' + item + '.pep',
    #                         shell=True)
    #         subprocess.call('cat ' + date + '/' + item + '/*/*.cmt > ' + date + '/' + item + '/' + item + '.cmt',
    #                         shell=True)
    #         subprocess.call('tbl2asn -p ' + date + '/' + item + ' -t ' + sbt_file_loc + ' -Y ' + date + '/' + item +
    #                         '.cmt -a s -V vb', shell=True)

    print('Done, did  ' + str(len(virus_strain_list)) + ' viruses in ' + str(timeit.default_timer() - start_time) +
          ' seconds')
