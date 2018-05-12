# VAPiD v1.2

# VAPiD is an extremely lightweight virus genome annotator that takes any number of viral genomes and annotates them
# producing files suitable for NCBI submission

import subprocess
import re
import argparse
import timeit
import os
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
import platform
import sys
from Bio import Entrez
import time
Entrez.email = 'uwvirongs@gmail.com'


VERSION = 'v1.0'


# Reads in a fasta file that should have strain names for the names of the sequences -  can handle any number of
# sequences. Also strips leading and trailing Ns or ?s from the provided sequence. Returns two lists with the names of
# the strains in the first one and the genomes as strings in the second list, also changes U's to T's
def read_fasta(fasta_file_loc):
    strain_list = []
    genome_list = []
    dna_string = ''

    for line in open(fasta_file_loc):
        if line[0] == '>':
            strain_list.append(line[1:].split()[0])
            if dna_string != '':
                # strip leading and trailing Ns or ?s because there's no reason to submit them
                xip = 0
                while dna_string[xip] == 'N' or dna_string[xip] == '?':
                    xip += 1

                y = len(dna_string)
                while dna_string[y-1] == 'N' or dna_string[y-1] == '?':
                    y -= 1

                dna_string = dna_string[xip:y]

                genome_list.append(dna_string)
                dna_string = ''
        else:
            dna_string += line.strip()
    # Just to make sure all our sequences are on the same page
    genome_list.append(dna_string.replace('U', 'T'))
    return strain_list, genome_list


# This function takes the strain name and a location of the individual fasta file we saved earlier and runs a blast
# Search saving the top 35 hits - then the top hit is found that is a complete genome and the fasta and .gbk of that
# are saved - then we run alignment on the two and return two strings one of them our sequence with '-' and the
# other of our new reference sequence
def blast_n_stuff(strain, our_fasta_loc):
    # If we've already done this one before skip the blasting step, speeds up error checking for development
    # if user provided own reference use that one - also use our specifically chosen reference for some viruses
    if args.r:
        ding = args.r
        ref_seq_gb = ding

    elif args.db:
        local_database_location = args.db
        print('Searching local blast database at ' + local_database_location)
        # may need to tweak the output method - need to test first
        local_blast_cmd = 'blastn -db ' + local_database_location + ' -query ' + our_fasta_loc + \
                          ' -num_alignments 1 -word_size 28 -outfmt 6 -out ' + strain + SLASH + strain \
                          + '.blastresults'
        subprocess.call(local_blast_cmd, shell=True)

        # pull first accession number
        for line in open(strain + SLASH + strain + '.blastresults'):
            ref_seq_gb = line.split('|')[3]
            break

    elif args.online:
        if not os.path.isfile(strain + SLASH + strain + '.blastresults'):
            print('Searching NCBI for the best reference sequence (may take longer for multiple requests due to NCBI '
                      'throttling)')

            record = open(our_fasta_loc).read()

            result_handle = NCBIWWW.qblast('blastn', 'nt', record, word_size=28, descriptions=0, alignments=15, entrez_query = 'txid10239[ORGN]',
                                            format_type='Text')
            with open(strain + SLASH + strain + '.blastresults', 'w') as out_handle:
                out_handle.write(result_handle.read())
            result_handle.close()

            # read through the top hits saved earlier and save the accession number of the best hit that's complete
            read_next = False
            found = False
            for line in open(strain + SLASH + strain + '.blastresults'):
                if line[0] == '>':
                    name_of_virus = ' '.join(line.split()[1:]).split('strain')[0].split('isolate')[0].strip()
                    name_of_virus = name_of_virus.split('/')[0]
                    ref_seq_gb = line.split()[0][1:]

                    # last part of these two logic checks is so we avoid the misassembled/mutated viruses
                    # This is going to get really out of hand if we have to keep blacklisting records
                    if 'complete' in line and ref_seq_gb.split('.')[0] not in 'KM551753 GQ153651 L08816 HIVANT70C L20587':
                        found = True
                        break
                    else:
                        read_next = True
                elif read_next:
                    if 'complete genome' in line and ref_seq_gb.split('.')[0] not in 'KM551753 GQ153651 L08816 HIVANT70C L20587':
                        found = True
                        break
                    else:
                        read_next = False

            # if we don't find any complete genomes just pull the top hit from blast and go from there
            if not found:
                for line in open(strain + SLASH + strain + '.blastresults'):
                    if line[0] == '>':
                        name_of_virus = ' '.join(line.split()[1:]).split('strain')[0].split('isolate')[0].strip()
                        ref_seq_gb = line.split()[0][1:]
                        break
    # default case -- use either of the provided reference databases that we will include
    else:
        if os.path.isfile('all_virus.fasta'):
            local_database_location = 'all_virus.fasta'
        elif os.path.isfile('virus_compressed.fasta'):
            local_database_location = 'virus_compressed.fasta'
        else:
            print('No local blast database found in this folder! Please install from the github releases page! Or use vapid with --online')
            print('Exiting...')
            exit(0)
        print('Searching local blast database at ' + local_database_location)

        local_blast_cmd = 'blastn -db ' + local_database_location + ' -query ' + our_fasta_loc + \
                          ' -num_alignments 1 -word_size 28 -outfmt 6 -out ' + strain + SLASH + strain \
                          + '.blastresults'

        subprocess.call(local_blast_cmd, shell=True)

        # pull first accession number
        for line in open(strain + SLASH + strain + '.blastresults'):
            ref_seq_gb = line.split('|')[3]
            break

    # Download the reference fasta file from Entrez
    record = Entrez.read(Entrez.esearch(db='nucleotide', term=ref_seq_gb))
    # Download .gbk from Entrez, we'll pull annotations from this file later
    h2 = Entrez.efetch(db='nucleotide', id=record["IdList"][0], rettype='gb', retmode='text')
    e = open(strain + SLASH + strain + '_ref.gbk', 'w')
    e.write(h2.read())
    e.close()
    time.sleep(1)

    if not args.online:
        for line in open(strain + SLASH + strain + '_ref.gbk'):
            if 'DEFINITION' in line:
                name_of_virus = ' '.join(line.split()[1:]).split('strain')[0].split('isolate')[0].strip()

    print(ref_seq_gb + ' was the selected reference')
    print(name_of_virus + ' was the parsed name of the virus')

    h = Entrez.efetch(db='nucleotide', id=record["IdList"][0], rettype='fasta', retmode='text')
    d = open(strain + SLASH + strain + '_ref.fasta', 'w')
    d.write(h.read())
    d.close()

    # mafft rules and the developer of mafft is awesome
    z = open(strain + SLASH + strain + '_aligner.fasta', 'w')
    fe = open(our_fasta_loc)
    for line in fe:
        z.write(line)
    fe.close()
    ge = open(strain + SLASH + strain + '_ref.fasta')
    z.write('\n')
    for line in ge:
        z.write(line)
    ge.close()
    z.close()

    # Windows
    if SLASH == '\\':
        s = 'mafft-win\\mafft.bat --quiet ' + strain + SLASH + strain + '_aligner.fasta > ' + strain + SLASH + strain + '.ali'
        subprocess.call(s, shell=True)
    else:
        try:
            subprocess.call('mafft --quiet ' + strain + SLASH + strain + '_aligner.fasta > ' + strain + SLASH + strain + '.ali',
                    shell=True)
        except:
            print('Running on a non windows system, which means you need to install mafft and put it on the sys path '
                  'yourself.\nI suggest using brew or apt')
    ali_list, ali_genomes = read_fasta(strain + SLASH + strain + '.ali')

    # SWAPPED THESE DURING DEBUGGING
    ref_seq = ali_genomes[1]
    our_seq = ali_genomes[0]

    return name_of_virus, our_seq, ref_seq


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
def adjust(given_num, our_num_array, ref_num_array, genome):

    # Handles gene lengths that go off the end of the genome
    if given_num == len(genome):
        return len(genome)

    # Go through our number array and search for the number of interest
    if our_num_array[given_num] == '-1':

        in_dex = given_num
        while our_num_array[in_dex != '-1']:
            in_dex += 1
            break
        return str(our_num_array[in_dex])

    else:
        found = False
        for x in range(0, len(our_num_array)):
            if ref_num_array[x] == given_num:
                index = x
                found = True
                break

    # now index is the absolute location of what we want
    if found:
        return str(our_num_array[index])
    else:
        return str(len(genome))


# this opens up the reference .gbk file and pulls all of the annotations, it then adjusts the annotations to the
# relative locations that they should appear on our sequence
def pull_correct_annotations(strain, our_seq, ref_seq, genome):

    # Read the reference gbk file and extract lists of all of the protein locations and annotations!
    gene_loc_list = []
    gene_product_list = []
    allow_one = False

    for line in open(strain + SLASH + strain + '_ref.gbk'):
        if ' CDS ' in line:
            # this is now going to be a list of numbers, start-stop start-stop
            # this line simply makes sure we read in reversed start-stops in the true reversed direction
            if 'complement' in line:
                whack = re.findall(r'\d+', line)
                whack.reverse()
                gene_loc_list.append(whack)
            else:
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
            gene_loc_list[entry][y] = adjust(int(gene_loc_list[entry][y]), our_seq_num_array, ref_seq_num_array, genome)
    return gene_loc_list, gene_product_list


# takes a strain name and a genome and writes and saves a fasta to the correct directory
def write_fasta(strain, genome):
    w = open(strain + SLASH + strain + '.fasta', 'w')
    w.write('>' + strain + '\n')
    w.write(genome)
    w.close()


# Take the name of a virus sample, and write the .cmt file for it using supplied coverage information
# NOTE: only writes coverage length - so now if we want to say our sequencing platform we have to edit this code
def write_cmt(sample_name, coverage):
    cmt = open(sample_name + SLASH +'assembly.cmt', 'w')
    cmt.write('##Assembly-Data-START##\n')
    if coverage != '':
        cmt.write('Coverage\t' + coverage + '\n')
    cmt.write('##Assembly-Data-END##\n')
    cmt.close()


# this takes in all of our information and makes a feature table that actually should work
# annotations for ribosomal slippage and RNA editing should be working now - as well as creation of a .pep file for rna
# editing
def write_tbl(strain, gene_product_list, gene_locations, genome, gene_of_intrest, note, name_o_vir):

    tbl = open(strain + SLASH + strain + '.tbl', 'w')
    tbl.write('>Feature ' + strain)

    for x in range(0, len(gene_product_list)):
        print(gene_product_list[x] + ' ' + str(gene_locations[x]))
        flag = ''
        xtra = ''
        sflag = ''
        product = gene_product_list[x]

        if gene_of_intrest in product:
            xtra = note

        if 'HIV' in name_o_vir and ('Pol polyprotein' == product or 'Pol' == product):
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
            if 'HEPATITIS B' not in name_o_vir and 'BK polyamavirus' not in name_o_vir:
                tbl.write('\t\t\texception\tRibosomal Slippage\n')

        else:
            start = int(location_info[0])
            end = int(location_info[1])

            if end >= len(genome) and genome[end - 3:end].upper() not in 'TGA,TAA,TAG,UGA,UAA,UAG':
                flag = '>'

            it_count = 0
            modifid_orf = False
            #won't execute this block of code for complemented genes...
            if end > start:
                while genome[end - 3:end].upper() not in 'TGA,TAA,TAG,UGA,UAA,UAG' and end < len(genome) - 3 and it_count <= 3:
                    print('Modifying ORF length for ' + str(product))
                    modifid_orf = True
                    end += 3
                    it_count += 1
                if modifid_orf and genome[end - 3:end].upper() not in 'TGA,TAA,TAG,UGA,UAA,UAG':
                    end -= it_count * 3

            # This should now correctly annotate assemblies that come in with the very edges chopped off
            pie = ''
            die = ''
            if int(start) < 1:
                sflag = '<'
                pie = str((int(end) % 3) + 1)
                start = '1'

            if int(end) < 1:
                end = len(genome)
                flag = '>'

            if 'HPIV-1' in name_o_vir or 'human parainfluenza virus 1' in name_o_vir.lower():
                if 'C\'' in product or 'Y2' in product:
                    die = '\n\t\t\ttransl_except\t(pos:' + str(start) + '..' + str(int(start) + 2) + ',aa:Met)'
            tbl.write('\n' + sflag + str(start) + '\t' + flag + str(end) + '\tCDS\n')
            tbl.write('\t\t\tproduct\t' + product + xtra)
            if pie != '':
                tbl.write('\n\t\t\tcodon_start\t' + pie)
            if die != '':
                tbl.write(die)

    tbl.write('\n')
    tbl.close()


# takes a single strain name and a single genome and annotates and save the entire virus and annotations package
# returns the "species" of the virus for consolidated .sqn packaging
def annotate_a_virus(strain, genome, metadata, coverage, sbt_loc):

    if not os.path.exists(strain):
        os.makedirs(strain)

    write_fasta(strain, genome)

    name_of_virus, our_seq, ref_seq = blast_n_stuff(strain, strain + SLASH + strain + '.fasta')

    gene_loc_list, gene_product_list = pull_correct_annotations(strain, our_seq, ref_seq, genome)

    write_cmt(strain, coverage)

    write_fsa(strain, name_of_virus, genome, metadata)

    extra_stuff = ''

    # prime gene of interest so unless we're in one of the specific cases nothing will trigger
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

    if 'parainfluenza virus 4' in name_of_virus.lower():
        extra_stuff = '\n\t\t\texception\tRNA Editing\n\t\t\tnote\tRNA Polymerase adds 2 non templated ' \
                      'G\n\t\t\tprotein_id\tn_' + strain
        gene_of_interest = 'phosphoprotein'
        if 'P' in gene_product_list:
            gene_of_interest = 'P'
        elif 'P protein' in gene_product_list:
            gene_of_interest = 'P protein'
        process_para(strain, genome, gene_loc_list, gene_product_list, gene_of_interest, 'HPIV4a')

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

    write_tbl(strain, gene_product_list, gene_loc_list, genome, gene_of_interest, extra_stuff, name_of_virus)

    cmd = 'tbl2asn -p ' + strain + SLASH + ' -t ' + sbt_loc + ' -Y ' + strain + SLASH + 'assembly.cmt -V vb'
    try:
        subprocess.call(cmd, shell=True)
    except:
        print('tbl2asn not installed, go to https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/ and download the appropriate version')
    return name_of_virus


# This function takes a nucleotide sequence that has non-templated G's inserted an unknown number of times and trims
# Them from the end to read into a correct frame  - translates the sequence and picks the one with the least number of
# stop codons returns the raw sequence
def pick_correct_frame(one, two):
    print('Picking correct reading frame for RNA editing:')
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
    print('Looks like this virus has RNA editing, fixing it now')
    found_ = False
    for g in range(0, len(gene_product_list)):
        # flipping this covers whack spacing in protein products
        if gene_of_interest in gene_product_list[g]:
            nts_of_gene = genome[int(gene_loc_list[g][0]) - 1:int(gene_loc_list[g][1]) - 1]
            found_ = True
            break

    if found_:

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

        elif v == 'HPIV4a':
            start_of_poly_g = nts_of_gene.find('AAGAGG', 439, 460)
            nts_of_gene = nts_of_gene[0:start_of_poly_g + 1] + 'GG' + nts_of_gene[start_of_poly_g + 1:]

        new_translation = str(Seq(nts_of_gene).translate())

        pep = open(strain + SLASH + strain + '.pep', 'w')
        pep.write('>n_' + strain + '\n' + new_translation)
        pep.write('\n')
        pep.close()


# Writes an fsa file based of the name, strain and genome, honestly we should allow for much more flexibility
# and automation here
def write_fsa(strain, name_of_virus, virus_genome, metadata):
    fsa = open(strain + SLASH + strain + '.fsa', 'w')
    fsa.write('>' + strain + ' [organism=' + name_of_virus + ']' + '[moltype=genomic] [host=Human] [gcode=1] '
              '[molecule=RNA]' + metadata + '\n')
    fsa.write(virus_genome)
    fsa.write('\n')
    fsa.close()


# Build the metadata for every virus that's been submitted
def do_meta_data(strain, sheet_exists):
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
        print('metadata not found in provided .csv or .csv not created -  time for minimal manual entry for sequence - ' + strain)
        col = ' [collection-date=' + raw_input('Enter collection date in the format (23-Mar-2005, Mar-2005, or 2005): ').strip() + ']'
        con = ' [country=' + raw_input('Enter country sample was collected in (example - USA): ').strip() + ']'
        st = ' [strain=' + raw_input('Enter strain name - if unknown just put ' + strain + ': ').strip() + ']'
        cov = raw_input('Enter coverage as a number (example 42.3), if unknown just leave this blank and hit enter: ')
        meta_data = col + con + st
        coverage = cov
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
        print('WARNING: ' + sample_name + ' contains ' + str(stops) + ' stop codon(s)!')
        return True
    else:
        return False


# quick check to make sure slashes go the right way on both Windows and Mac/Linux
def check_os():
    if platform.system() == 'Linux' or platform.system() == 'Darwin':
        return '/'
    else:
        return '\\'


if __name__ == '__main__':

    start_time = timeit.default_timer()
    SLASH = check_os()
    parser = argparse.ArgumentParser(description='Version ' + VERSION + '\nPrepares FASTA file for NCBI Genbank submission '
                                                 'through local or online blastn-based annotation of viral sequences. '
                                                 'In default mode, VAPiD searches this folder for our viral databases.')
    parser.add_argument('fasta_file', help='Input file in .fasta format containing complete or near complete '
                                           'genomes for all the viruses that you want to have annotated')
    parser.add_argument('author_template_file_loc', help='File path for the NCBI-provided sequence author template file'
                        ' (should have a .sbt extension)\n https://submit.ncbi.nlm.nih.gov/genbank/template/submission/')
    parser.add_argument('--metadata_loc', help='If you\'ve input the metadata in the provided csv, specify the location '
                        'with this optional argument. Otherwise all metadata will be manually prompted for.')
    parser.add_argument('--r', help='If you want to specify a specific NCBI reference, put the accession number here '
                        '- must be the exact accession number - note: feature forces all sequences in FASTA to be this viral species.')
    parser.add_argument('--db', help='specify the local blast database name.  You MUST have blast+ with blastn'
                                    'installed correctly on your system path for this to work.')
    parser.add_argument('--online',action='store_true', help='Force VAPiD to blast against online database.  This is good for machines that don\'t '
                                         'have blast+ installed or if the virus is really strange.'
                                         'Warning: this can be EXTREMELY slow, up to ~5-25 minutes a virus')
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
    time = str(timeit.default_timer() - start_time)
    newtime = time.split('.')[0] + '.' + time.split('.')[1][0:1]
    print('Done, did  ' + str(len(virus_strain_list)) + ' viruses in ' + newtime + ' seconds')
