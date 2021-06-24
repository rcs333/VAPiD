# VAPiD is an extremely lightweight virus genome annotator that takes any number of viral genomes and annotates them
# producing files suitable for NCBI submission

# Vapid Version
VERSION = 'v1.6.7'

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
from Bio import SeqIO
import time
import shutil
import re
Entrez.email = 'uwvirongs@gmail.com'


# Reads in a fasta file that should have strain names for the names of the sequences -  can handle any number of
# sequences. Also strips leading and trailing Ns or ?s from the provided sequence. Also takes an optional boolean 
# value to determine if names should be allowed to have system protected characters like / in them. Returns three lists with the names of
# the strains in the first one and the genomes as strings in the second list, the third list is only used when the slashes argument is true
# also changes U's to T's
def read_fasta(fasta_file_loc, slashes=False):
    strain_list = []
    genome_list = []
    full_name_list = []
    dna_string = ''

    for line in open(fasta_file_loc):
        if line[0] == '>':
            if slashes:
                full_name_list.append(line[1:])
            else:
                full_name_list.append('')

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
    genome_list.append(dna_string)
    return strain_list, genome_list, full_name_list


# Spell checking functionality provided by Entrez
# takes an input string and returns the Entrez corrected string as long as it exists
def spell_check(query_string):
    # new entrez rules limit requests to no more than 3 a second, this will ensure we don't send more than two
    time.sleep(0.5)
    handle = Entrez.espell(term=query_string)
    record = Entrez.read(handle)
    corrected_query = record["CorrectedQuery"]
    # Entrez returns blank strings for numerals or things that are spelled correctly
    # Since this is based on NCBI's spell checking protein names are included and correct
    # However this won't correct SUPER messed up words or made up words
    if corrected_query != '':
        print('Checking spelling on ' + query_string)
        print(query_string + ' was corrected to: ' + corrected_query)
        return corrected_query
    else:
        return query_string


# This function takes the strain name and a location of the individual fasta file we saved earlier and runs a blast
# Search saving the top 15 hits - then the top hit is found that is a complete genome and the fasta and .gbk of that
# are saved - then we run alignment on the two and return two strings one of them our sequence with '-' and the
# other of our new reference sequence. If a sequence is submitted that needs to be reverse complemented to align
# to a reference the sequence will stay this way 
def blast_n_stuff(strain, our_fasta_loc):

    # if user provided own reference use that one - also use our specifically chosen reference for some viruses
    if args.r:
        ding = args.r
        ref_seq_gb = ding

    # if the user provided a database to use print location and use database
    elif args.db:
        local_database_location = args.db
        print('Searching local blast database at ' + local_database_location)
        # blastn with word size of 28 because we're searching off a provided reference we're just going to pull the top
        local_blast_cmd = 'blastn -db ' + local_database_location + ' -query ' + our_fasta_loc + \
                          ' -num_alignments 1 -word_size 28 -outfmt 6 -out ' + strain + SLASH + strain \
                          + '.blastresults'
        subprocess.call(local_blast_cmd, shell=True)

        # pull first accession number from our reference database
        for line in open(strain + SLASH + strain + '.blastresults'):
            ref_seq_gb = line.split('|')[3]
            break

    # online search
    elif args.online:
        print('Searching NCBI for the best reference sequence (may take longer for multiple requests due to NCBI '
              'throttling)')

        record = open(our_fasta_loc).read()

        # used to have entrez_query = 'txid10239[ORGN]' as an argument to restrict searches to taxonomic group 'Viruses'
        result_handle = NCBIWWW.qblast('blastn', 'nt', record, word_size=28, descriptions=0, alignments=15,
                                            format_type='Text')
        with open(strain + SLASH + strain + '.blastresults', 'w') as out_handle:
            out_handle.write(result_handle.read())
        result_handle.close()

        # read through the top hits saved earlier and save the accession number of the best hit that's complete
        read_next = False
        found = False
        for line in open(strain + SLASH + strain + '.blastresults'):
            if line[0] == '>':
                name_of_virus = ' '.join(line.split()[1:]).split('strain')[0].split('isolate')[0].split('complete')[0].split('partial')[0].split('genomic')[0].split('from')[0].strip()
                name_of_virus = name_of_virus.split('/')[0]
                ref_seq_gb = line.split()[0][1:]

                # last part of these two logic checks is so we avoid the misassembled/mutated viruses
                # This is going to get really out of hand if we have to keep blacklisting records
                # TODO: pull these out of the database, regenerate and reupload 
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
                    name_of_virus = ' '.join(line.split()[1:]).split('strain')[0].split('isolate')[0].split('complete')[0].split('partial')[0].split('genomic')[0].split('from')[0].strip()
                    ref_seq_gb = line.split()[0][1:]
                    break
    # default case -- use either of the provided reference databases that we will include
    else:
        # all virus is the preferred database but we'll switch to compressed if that's what the user downloaded
        # this list IS ordered by how much I recommend using these databases
        if os.path.isfile('all_virus.fasta.nin'):
            local_database_location = 'all_virus.fasta'
        elif os.path.isfile('virus_compressed.fasta.nin'):
            local_database_location = 'virus_compressed.fasta'
        elif os.path.isfile('ref_seq_vir.nin'):
            local_database_location = 'ref_seq_vir'
        # print a helpful error message and exit
        else:
            print('No local blast database found in this folder! Please install from the github releases page! '
                '(https://github.com/rcs333/VAPiD/releases) Or use vapid with --online (not reccomended)')
            print('Exiting...')
            exit(0)
        print('Searching local blast database at ' + local_database_location)

        # we're only going to save one because these our pretty decent reference databases
        local_blast_cmd = 'blastn -db ' + local_database_location + ' -query ' + our_fasta_loc + \
                          ' -num_alignments 1 -word_size 28 -outfmt 6 -out ' + strain + SLASH + strain \
                          + '.blastresults'

        subprocess.call(local_blast_cmd, shell=True)

        # pull first accession number
        for line in open(strain + SLASH + strain + '.blastresults'):
            ref_seq_gb = line.split('|')[3]
            break

    # Download the reference fasta file from 
    record = Entrez.read(Entrez.esearch(db='nucleotide', term=ref_seq_gb))
    # Download .gbk from Entrez, we'll pull annotations from this file later
    h2 = Entrez.efetch(db='nucleotide', id=record["IdList"][0], rettype='gb', retmode='text')
    e = open(strain + SLASH + strain + '_ref.gbk', 'w')
    e.write(h2.read())
    e.close()

    # if we've specified our own file remove the one we just found and delete it 
    if args.f:
        os.remove(strain + SLASH + strain + '_ref.gbk')
        shutil.copyfile(args.f, strain + SLASH + strain + '_ref.gbk')

    # NCBI online tools don't want more than like 1 request every .2 seconds or something so we just sleep for a second here
    time.sleep(1)

    # because the difference in how this stuff gets saved we have to pull online differently
    if not args.online:
        for line in open(strain + SLASH + strain + '_ref.gbk'):
            if 'DEFINITION' in line:
                # this forces removal of 'complete/partial annotations because those get added by genbank and there is no reason to include them'
                # we also want to strip off specific strain and isolate names in order to tend towards being more general 
                name_of_virus = ' '.join(line.split()[1:]).split('strain')[0].split('isolate')[0].split('complete')[0].split('partial')[0].split('genomic')[0].split('from')[0].strip()
    # let the user know 
    print(ref_seq_gb + ' was the selected reference')
    print(name_of_virus + ' was the parsed name of the virus')


    if args.f:
        SeqIO.convert(args.f, "genbank", strain  + SLASH + strain + "_ref.fasta", "fasta")
    else:
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
    print('Aligning reference and query...')
    # Windows
    if SLASH == '\\':
        # since we include the windows installation of maaft with vapid we can hard code the path
        s = 'mafft-win\\mafft.bat --adjustdirection --quiet ' + strain + SLASH + strain + '_aligner.fasta > ' + strain + SLASH + strain + '.ali'
        subprocess.call(s, shell=True)
    else:
        try:
            subprocess.call('mafft --adjustdirection --quiet ' + strain + SLASH + strain + '_aligner.fasta > ' + strain + SLASH + strain + '.ali 2>/dev/null',
                    shell=True)
        # print a helpful error message and exit
        except:
            print('Running on a non windows system, which means you need to install mafft and put it on the sys path '
                  'yourself.\nI suggest using brew or apt')
            exit(0)
    
    ali_list, ali_genomes, dumy_var_never_used = read_fasta(strain + SLASH + strain + '.ali')

    need_to_rc = False 
    # this will create a weird 
    if '_R_' in ali_list[1]:
        if ali_list[1][0:3] == '_R_':
            # we need to RC input query and redo 
            need_to_rc = True
    print('Done alignment')
    # this is the reverse of what I expect but it works
    ref_seq = ali_genomes[1]
    our_seq = ali_genomes[0]

    # now also returning the accession of the reference for use in the .cmt file as well as a bool for if we need to rerun this block 
    return name_of_virus, our_seq, ref_seq, ref_seq_gb, need_to_rc


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
    found = False
    # Handles gene lengths that go off the end of the genome
    # 1.6.4 - this is obsolete and a bad implementation, the block at the end of this function takes care of this
    # better, I'm leaving this in comments for a while just in case 
    #if given_num >= len(genome):
    #    return len(genome)

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
        if our_num_array[index] >= len(genome):
        # this is the new handling of when genes run off the end of the submitted genome 
            return len(genome)
        return str(our_num_array[index])
    else:
        return str(len(genome))


# this opens up the reference .gbk file and pulls all of the annotations, it then adjusts the annotations to the
# relative locations that they should appear on our sequence
def pull_correct_annotations(strain, our_seq, ref_seq, genome):
    # Read the reference gbk file and extract lists of all of the protein locations and annotations!

    # now we're doing this at the top because we'recalling this earlier
    our_seq_num_array, ref_seq_num_array = build_num_arrays(our_seq, ref_seq)

    gene_loc_list = []
    gene_product_list = []
    allow_one = False

    all_loc_list = []
    all_product_list = []
    name_of_the_feature_list = []
    # Experimental code for transferring 'gene' annotations from NCBI reference sequence
    if args.all:
        for line in open(strain + SLASH + strain + '_ref.gbk'):
            if ('..' in line) and ( ('gene' in line) or  ('mat_peptide' in line) or ('UTR' in line) or ('repeat_region' in line)  ):
                name_of_the_feature_list.append(line.split()[0])

                if 'complement' in line:
                    whack = re.findall(r'\d+', line.split()[1])
                    whack.reverse()
                    all_loc_list.append(whack)
                else:
                    all_loc_list.append(re.findall(r'\d+', line.split()[1]))
                allow_one = True
                if 'UTR' in line:
                    allow_one = False
                    all_product_list.append('')
            elif allow_one : 
                if '/product' in line:
                    allow_one = False
                    px_all = line.split('=')[1][1:-2]
                    all_product_list.append(px_all)
                elif '/gene' in line:
                    allow_one = False
                    px_all = line.split('=')[1][1:-2]
                    all_product_list.append(px_all)
                elif 'UTR' in name_of_the_feature_list[-1]:
                    px_all = line.split('=')[1][1:-2]
                    all_product_list.append(px_all)
                    allow_one = False
                elif '/rpt_type' in line:
                    px_all = line.split('=')[1][0:-1]
                    all_product_list.append(px_all)
                    allow_one = False
                elif '/db_xref' in line:
                    name_of_the_feature_list.pop()
                    all_loc_list.pop()
                    allow_one = False

        # adjust gene list
        #print(all_loc_list)
        #print(all_product_list)
        #print(name_of_the_feature_list)
        for entry in range(0, len(all_loc_list)):
            for y in range(0, len(all_loc_list[entry])):
                all_loc_list[entry][y] = adjust(int(all_loc_list[entry][y]), our_seq_num_array, ref_seq_num_array, genome)
    #print("DONE WITH THE ALL STUFF")
    allow_one = False
    for line in open(strain + SLASH + strain + '_ref.gbk'):
        if ' CDS ' in line and '..' in line:
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
            # for some weird reason - this is the way we only go in when the flag is passed 
            # TODO: need to make sure that this is activating properly 
            if args.spell_check:
                new_list = []
                px_word_list = px.split()
                for word in px_word_list:
                    if '1' or '2' or '3' or '4' or '5' or '6' or '7' or '8' or '9' or '0' not in word: 
                        new_list.append(spell_check(word))

                px = ' '.join(new_list)
            if px == 'phospho protein':
                px = 'phoshoprotein'
            gene_product_list.append(px)

    # Adjust every locus so that we actually put in correct annotations
    for entry in range(0, len(gene_loc_list)):
        for y in range(0, len(gene_loc_list[entry])):
            gene_loc_list[entry][y] = adjust(int(gene_loc_list[entry][y]), our_seq_num_array, ref_seq_num_array, genome)
    return gene_loc_list, gene_product_list, all_loc_list, all_product_list, name_of_the_feature_list


# takes a strain name and a genome and writes and saves a fasta to the correct directory
def write_fasta(strain, genome):
    w = open(strain + SLASH + strain + '.fasta', 'w')
    w.write('>' + strain + '\n')
    w.write(genome)
    w.close()


# Take the name of a virus sample, and write the .cmt file for it using supplied coverage information
# NOTE: only writes coverage length - so now if we want to say our sequencing platform we have to edit this code
# Now also writes in the comment the reference that this subission was annotated off - this should provide some more
# accountability
def write_cmt(sample_name, coverage, ref_gb, did_we_rc):
    cmt = open(sample_name + SLASH + 'assembly.cmt', 'w')
    cmt.write('##Assembly-Data-START##\n')
    if coverage != '':
        cmt.write('Coverage\t' + coverage + '\n')
    if did_we_rc:
        cmt.write('Original input sequence was reverse complemented by MAFFT during the alignment phase')
    cmt.write('Created with VAPiD' + VERSION + ' Reference annotations were pulled from ' + ref_gb + '\n')
    cmt.write('##Assembly-Data-END##\n')
    cmt.close()


# this takes in all of our information and makes a feature table that contains correct annotations for for ribosomal slippage and RNA editing
# - as well as creation of a .pep file for rna editing -- Now we also pass two possibly empty lists to write tbl so we can write gene annotations
def write_tbl(strain, gene_product_list, gene_locations, genome, gene_of_intrest, note, name_o_vir, all_loc_list, all_product_list, full_name, name_of_the_feature_list):
    # covers the nipah situation where there's RNA editing on more than 1 protein - if this happens for more viruses I'll need to code a more
    # robust sollution, but for now this works 
    if 'nipah' in name_o_vir.lower():
        pep = open(strain + SLASH + strain + '.pep', 'w')

    tbl = open(strain + SLASH + strain + '.tbl', 'w')
    tbl.write('>Feature ' + full_name)

    # This block should write all gene annotations to tbl file as long as we got passed genes, and the only way that will ever happen is if the
    # User put the -all flag
    if len(all_product_list) > 0:
        for x in range(0, len(all_product_list)):
            print(all_product_list[x] + str(all_loc_list[x]))

            e_flag = ''
            s_flag = ''
            s_all = all_loc_list[x][0]
            e_all = all_loc_list[x][1]
            p_all = all_product_list[x]
            if int(e_all) >= len(genome):
                e_flag = ''
            if int(s_all) < 1:
                s_flag = '<'
                s_all = '1'

            if int(e_all) < 1:
                e_all = len(genome)
                e_flag = '>'
            if p_all == 'inverted terminal repeat':
                if int(s_all) <  (len(genome) / 2):
                    s_all = 1
                else:
                    e_all = len(genome)
                tbl.write('\n' + s_flag + str(s_all) + '\t' + e_flag + str(e_all) + '\trepeat_region\n')
                tbl.write('\t\t\tnote\t' + p_all + '\n')
                tbl.write('\t\t\trpt_type\tinverted')
            else:
                tbl.write('\n' + s_flag + str(s_all) + '\t' + e_flag + str(e_all) + '\t' +name_of_the_feature_list[x] + '\n')
                if 'UTR' in name_of_the_feature_list[x] or 'gene' in name_of_the_feature_list[x]:
                    feat_des = 'gene'
                elif 'mat_peptide' in name_of_the_feature_list[x] or 'CDS' in name_of_the_feature_list:
                    feat_des = 'product'
                else:
                    feat_des =  'rpt_type'
                tbl.write('\t\t\t' + feat_des + '\t' + p_all)
    
    for x in range(0, len(gene_product_list)):
        print(gene_product_list[x] + ' ' + str(gene_locations[x]))
        flag = ''
        xtra = ''
        sflag = ''
        product = gene_product_list[x]

        if gene_of_intrest in product:
            xtra = note

        if 'nipah' in name_o_vir.lower():
            nts_of_gene = genome[int(gene_locations[x][0]) - 1:int(gene_locations[x][1]) - 1]
            
            if product.lower() == 'v protein':
                xtra = '\n\t\t\texception\tRNA Editing\n\t\t\tnote\tRNA Polymerase adds non templated ' \
                      'G\n\t\t\tprotein_id\tn_1' + strain
                start_of_poly_g = nts_of_gene.find('AAAAAGG')
                nts_of_gene = nts_of_gene[0:start_of_poly_g + 1] + 'G' + nts_of_gene[start_of_poly_g + 1:]
                new_translation = str(Seq(nts_of_gene).translate())
                pep.write('>n_1' + strain + '\n' + new_translation)
                pep.write('\n')
            elif product.lower() == 'w protein':
                xtra = '\n\t\t\texception\tRNA Editing\n\t\t\tnote\tRNA Polymerase adds non templated ' \
                      'G\n\t\t\tprotein_id\tn_2' + strain
                start_of_poly_g = nts_of_gene.find('AAAAAGG')
                nts_of_gene = nts_of_gene[0:start_of_poly_g + 1] + 'GG' + nts_of_gene[start_of_poly_g + 1:]
                new_translation = str(Seq(nts_of_gene).translate())
                pep.write('>n_2' + strain + '\n' + new_translation)
                pep.write('\n')

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

            it_count = 0
            modifid_orf = False
            # won't execute this block of code for complemented genes
            #print(genome[end - 3:end].upper())
            # this makes sure that our end is in frame, and will adjust hopefully this doesn't break everything 
            # added a check to only do this in the absence of RNA editing or ribosomal sliippage 
            if xtra == '' and end != len(genome):
                if ((end - start) + 1) % 3 != 0:
                    end +=1
                if ((end - start) + 1) % 3 != 0:
                    end +=1

            
            if end > start and 'IIIA' not in product.upper():
                if (genome[end - 3:end].upper() not in 'TGA,TAA,TAG,UGA,UAA,UAG') and (end < len(genome) - 3) and not re.search('[MRWSYKVHDBN]',genome[end - 3:end].upper()):
                    if re.search('[MRWSYKVHDBN]',genome[end - 3:end].upper()):
                        print('Ambiguous base detected in a putative stop codon, this can cause problems with VAPiD annotations')
                    print('Modifying ORF length for ' + str(product))
                    end = find_end_stop(genome, start, end)
            # This should now correctly annotate assemblies that come in with the very edges chopped off
            #print(genome[end - 3:end].upper())
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
    if 'nipah' in name_o_vir.lower():
        pep.close()



# Takes a nucleotide sequence and a start and end position [1 indexed] and search for a stop codon from the start
# to the end + 60 so every codon in the provided gene and then 3 after it. Return the first stop codon found or if no
# stop codon is found return the original end value and print a warning
def find_end_stop(genome, start, end):
    # save the provided end
    start -= 1 
    old_end = end
    end = start + 3
    # Search for stop codons in DNA and RNA space until 3 codons after the provided end.
    # Turns out 3 codons isn't enough 
    while genome[end -3:end].upper() not in 'TGA,TAA,TAG,UGA,UAA,UAG' and end <= (old_end + 60):
        end += 3
    if end == old_end + 60:
        print('WARNING no stop codon found, examine reference and original sequence')
        return old_end
    else:
        return end


# takes a single strain name and a single genome and annotates and save the entire virus and annotations package
# returns the "species" of the virus for consolidated .sqn packaging
def annotate_a_virus(strain, genome, metadata, coverage, sbt_loc, full_name, nuc_a_type):
    did_we_reverse_complement = False 
    if not os.path.exists(strain):
        os.makedirs(strain)

    if '_R_' in strain:
        if strain[0:3] == '_R_':
            print('WARNING: ' + strain + ' has _R_ as the first part of the sequence charachters YOU HAVE TO CHANGE THIS')
    write_fasta(strain, genome)

    name_of_virus, our_seq, ref_seq, ref_accession, need_to_rc = blast_n_stuff(strain, strain + SLASH + strain + '.fasta')

    if need_to_rc:
        print('Input sequence needed to be reverse complemented to align properly.')
        new_seq = Seq(genome)
        # reverse complement input sequence and overwrite variable
        genome = str(new_seq.reverse_complement())
        did_we_reverse_complement = True
        # overwrite our fasta and this happens before we call write fsa 
        write_fasta(strain, genome)

        name_of_virus, our_seq, ref_seq, ref_accession, need_to_rc = blast_n_stuff(strain, strain + SLASH + strain + '.fasta')

    gene_loc_list, gene_product_list, all_loc_list, all_product_list, name_of_the_feature_list = pull_correct_annotations(strain, our_seq, ref_seq, genome)

    write_cmt(strain, coverage,ref_accession, did_we_reverse_complement)

    write_fsa(strain, name_of_virus, genome, metadata, full_name, nuc_a_type)
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
        #print('start of gene is ' + str(put_start))
        orf = genome[put_start - 1:put_start + 4000]
        #print('orf length is ' + str(len(orf)))
        #print('genome length is ' + str(len(genome)))
        orf_trans = str(Seq(orf).translate())
        #print('orf translation is ' + orf_trans)
        if orf_trans.find('*') != -1:
            put_end = (orf_trans.find('*') * 3)
            print('putative end is ' + str(put_end))
            gene_loc_list[7][1] = put_start + put_end

    if 'parainfluenza virus 2' in name_of_virus.lower() or 'rubulavirus 2' in name_of_virus.lower():
        #print('Custom code for HPIV2 runnning')
        extra_stuff = '\n\t\t\texception\tRNA Editing\n\t\t\tnote\tRNA Polymerase adds 2 non templated ' \
                      'G\n\t\t\tprotein_id\tn_' + strain
        gene_of_interest = 'P protein'
        process_para(strain, genome, gene_loc_list, gene_product_list, gene_of_interest, 'HPIV2')

    


    write_tbl(strain, gene_product_list, gene_loc_list, genome, gene_of_interest, extra_stuff, name_of_virus, all_loc_list, all_product_list, full_name, name_of_the_feature_list)

    cmd = 'tbl2asn -p ' + strain + SLASH + ' -t ' + sbt_loc + ' -Y ' + strain + SLASH + 'assembly.cmt -V vb '
    try:
        subprocess.call(cmd, shell=True)
    except:
        print('tbl2asn not installed, go to https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/ and download the appropriate version')
    print('Done with: ' + strain)
    print('')
    print('')
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
    #print('adding two Gs gives ' + str(two_count) + ' stop codon(s)')
    #print(two_trans)
    #print('adding one G gives ' + str(one_count) + ' stop codon(s)')
    #print(one_trans)
    if one_count < two_count:
        #print('chose one')
        return one
    else:
        #print('chose two')
        return two


# Takes a virus that has GGGGGG RNA editing and based on the gene that you send it and the name of the virus will
# find that gene in our annotations - add the correct number of G's and then translate the new 'mRNA' and write the
# translation to a .pep file where we can overwrite the sequin auto-translation
def process_para(strain, genome, gene_loc_list, gene_product_list, gene_of_interest, v):
    # Extract the gene protected because everything we throw in here are guaranteed to have the gene of interest
    print('Looks like this virus has RNA editing, fixing it now')
    #print(v)

    found_ = False
    #print('gene of interest = '+ gene_of_interest)
    for g in range(0, len(gene_product_list)):
        # flipping this covers whack spacing in protein products
        #print('product = '+ gene_product_list[g])
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

        elif v == 'HPIV2':
            #print('HPIV2 is getting here correctly')
            start_of_poly_g = nts_of_gene.find('AAGAGG', 450, 490)
            nts_of_gene_1 = nts_of_gene[0:start_of_poly_g + 1] + 'G' + nts_of_gene[start_of_poly_g + 1:]
            nts_of_gene_2 = nts_of_gene[0:start_of_poly_g + 1] + 'GG' + nts_of_gene[start_of_poly_g + 1:]
            nts_of_gene = pick_correct_frame(nts_of_gene_1, nts_of_gene_2)

        new_translation = str(Seq(nts_of_gene).translate())

        pep = open(strain + SLASH + strain + '.pep', 'w')
        pep.write('>n_' + strain + '\n' + new_translation)
        pep.write('\n')
        pep.close()


# Writes an fsa file based of the name, strain and genome, honestly we should allow for much more flexibility
# and automation here
def write_fsa(strain, name_of_virus, virus_genome, metadata, full_name, nucleic_acid_type):
    fsa = open(strain + SLASH + strain + '.fsa', 'w')
    fsa.write('>' + full_name.strip() + ' [organism=' + name_of_virus + ']' + '[moltype=genomic] [host=Human] [gcode=1] '
              '[molecule=' + nucleic_acid_type + ']' + metadata + '\n')
    fsa.write(virus_genome)
    fsa.write('\n')
    fsa.close()


# Build the metadata for every virus that's been submitted
def do_meta_data(strain, sheet_exists, full_name):
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
                    elif names[dex].strip() == 'full_name':
                        if line.split(',')[dex].strip() != '':
                            full_name = line.split(',')[dex].strip()

                    else:
                        s = s + ' [' + names[dex].strip() + '=' + line.split(',')[dex].strip() + ']'
                break

    if s == '':
        print('metadata not found in provided .csv or .csv not created -  time for minimal manual entry for sequence - ' + strain)
        col = ' [collection-date=' + input('Enter collection date in the format (23-Mar-2005, Mar-2005, or 2005): ').strip() + ']'
        con = ' [country=' + input('Enter country sample was collected in (example - USA): ').strip() + ']'
        st = ' [strain=' + input('Enter strain name - if unknown just put ' + strain + ': ').strip() + ']'
        cov = input('Enter coverage as a number (example 42.3), if unknown just leave this blank and hit enter: ')
        meta_data = col + con + st
        coverage = cov
        # Here's one line of code to unilaterally standardize defualt naming scheme
        if full_name == '':
            full_name = strain + ' (' + con.split('=')[1][:-1] + '/' + col.split('=')[1][:-1] + ')'

    else:
        meta_data = s
        if full_name == '':
            full_name = strain
            print('Automatic strain naming failed but submission will proceed without metadata appended to the fasta header.')
    return meta_data, coverage, full_name


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
    parser.add_argument('--f', help='specify a custom gbf file that you would like to annotate off of')
    parser.add_argument('--db', help='specify the local blast database name.  You MUST have blast+ with blastn'
                                    'installed correctly on your system path for this to work.')
    parser.add_argument('--online',action='store_true', help='Force VAPiD to blast against online database.  This is good for machines that don\'t '
                                         'have blast+ installed or if the virus is really strange.'
                                         'Warning: this can be EXTREMELY slow, up to ~5-25 minutes a virus')
    parser.add_argument('--spell_check', action='store_true', help='Turn on spellchecking for protein annoations ')
    parser.add_argument('--all', action='store_true', help='Use this flag to transfer ALL annotations from reference, this is largely untested')
    parser.add_argument('--slashes', action='store_true', help='Use this flag to allow any characters in the name of your virus - This allows '
                        'you to submit with a fasta file formated like >Sample1 (Human/USA/2016/A) Complete CDS'
                        ' make sure that your metadata file only contains the first part of your name \'Sample1\' in the example above. '
                        'You can also submit names with slashes by specifying in the metadata sheet under the header full_name, if you do that '
                        'you do not need to use this flag')
    parser.add_argument('--dna', action='store_true', help='Make all files annotated by this run be marked as DNA instead of the default (RNA)')

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        sys.exit(0)

    fasta_loc = args.fasta_file
    if args.dna:
        nuc_acid_type = 'DNA'
    else:
        nuc_acid_type = 'RNA'
        
    sbt_file_loc = args.author_template_file_loc

    virus_strain_list, virus_genome_list, full_name_list = read_fasta(fasta_loc, args.slashes)

    strain2species = {}
    strain2stops = {}

    meta_list = []
    coverage_list = []
    if args.metadata_loc:
        metadata_sheet_location = args.metadata_loc
        for x in range(0, len(virus_strain_list)):
            metadata, coverage, full_name_list[x] = do_meta_data(virus_strain_list[x], True, full_name_list[x])
            meta_list.append(metadata)
            coverage_list.append(coverage)
    else:
        for x in range(0, len(virus_strain_list)):
            metadata, coverage, full_name_list[x] = do_meta_data(virus_strain_list[x], False, full_name_list[x])
            meta_list.append(metadata)
            coverage_list.append(coverage)

    for x in range(0, len(virus_strain_list)):
        strain2species[virus_strain_list[x]] = annotate_a_virus(virus_strain_list[x], virus_genome_list[x],
                                                                meta_list[x], coverage_list[x], sbt_file_loc, full_name_list[x],nuc_acid_type)
        # now we've got a map of [strain] -> name of virus (with whitespace)

    for name in virus_strain_list:
        # now we've got a map of [strain] -> boolean value if there are stops or not
        strain2stops[name] = check_for_stops(name)
        if len(name) > 23:
            print('WARNING: ' + name + ' is over 23 characters, which means that your gbf file will be corrupted')
    time = str(timeit.default_timer() - start_time)
    newtime = time.split('.')[0] + '.' + time.split('.')[1][0:1]
    print('Done, did  ' + str(len(virus_strain_list)) + ' viruses in ' + newtime + ' seconds')
