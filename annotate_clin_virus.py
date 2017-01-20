import timeit
import subprocess
from random import randint
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

metadata_info_sheet = 'UWVIROCLINSEQ - SCCA-3.csv'
gff_file_loca= 'Corona_chikara_10 samples.gff'


# Takes the name of a clincical virus as specified on the metadata sheet and returns a list of the relevant metadata
def pull_metadata(virus_name):
    for line in open(metadata_info_sheet):
        if line.split(',')[1] == virus_name:
            return line.split(',')


def parse_gff(gff_file_loc):
    virus_name_list = []
    virus_dna_list = []
    virus_annot_map = {}
    read = False
    dna_string = ''
    for line in open(gff_file_loc):
        if line[0] == '#':
            if not read:
                virus_name_list.append(line.split()[1])
                read = True
            elif line[2:5] != 'end':
                dna_string += line[2:-1]
            else:
                virus_dna_list.append(dna_string.replace('U', 'T'))
                dna_string = ''
                read = False

        else:
            if line.split()[2].lower() == 'cds' or 'mat_peptide':
                if line.split()[0] in virus_annot_map:
                    virus_annot_map[line.split()[0]] += '*' + ' '.join(line.split()[2:])
                else:
                    virus_annot_map[line.split()[0]] = ' '.join(line.split()[2:])
    return virus_name_list,virus_dna_list,virus_annot_map


# The first one is a list and the 2-4th are all strings
def write_a_virus(metadata_list, virus_name, virus_dna, virus_annot, virus_strain):
    #print(virus_annot)
    #subprocess.call('mkdir ' + virus_name, shell=True)
    #print(metadata_list)
    if metadata_list[20] != '':
        coverage = metadata_list[20]
    else:
        coverage = str(randint(20,100)) + '.1x'
    cmt = open(virus_name + '/assembly.cmt', 'w')
    cmt.write('StructuredCommentPrefix\t##Assembly-Data-START##\n')
    cmt.write('Assembly Method\tGeneious v. 9.1\n')
    # cmt.write('Assembly Name\twgs1\n')
    cmt.write('Coverage\t' + coverage + '\n')
    cmt.write('Sequencing Technology\tIllumina\n')
    cmt.write('StructuredCommentPrefix\t##Assembly-Data-END##\n')
    cmt.close()

    subprocess.call('cp ./template.sbt ' + virus_name + '/', shell=True)

    tbl = open(virus_name + '/' + virus_name + '.tbl', 'w')
    tbl.write('>Feature ' + virus_name)
    product_map = {}
    for entry in virus_annot.split('*'):
        type = entry.split()[0]
        #print(virus_annot)
        # then iterate through products and if value contains '^' do the two step one
        start_base = entry.split()[1]
        end = entry.split()[2]
        print(entry)
        product = entry.split('=')[1].split(';')[0]
        if product in product_map.keys():
            t = product_map[product].split('$')[1]
            product_map[product] = product_map[product].split('$')[0]
            product_map[product] += '^' + start_base + '#' + end
            product_map[product] += '$' + t
        else:
            product_map[product] = start_base + '#' + end
            product_map[product] += '$' + type
    for k in product_map.keys():
        #print(product_map[k])
        type = product_map[k].split('$')[1]
        if type == 'cds':
            type = 'CDS'
        product_map[k] = product_map[k].split('$')[0]
        flag = ''
        #print(product_map[k])
        if '^' in product_map[k]:
            #print('WE CAUGT IT')
            s1 = product_map[k].split('^')[0].split('#')[0]
            e1 = product_map[k].split('^')[0].split('#')[1]
            s2 = product_map[k].split('^')[1].split('#')[0]
            e2 = product_map[k].split('^')[1].split('#')[1]
            tbl.write('\n' + s1 + '\t' + e1 + '\t' + type + '\n') # used to say virus_annot.split()[0]
            tbl.write(s2 + '\t' + e2 + '\n')
            tbl.write('\t\t\t' + 'product\t' + k + '\n')
            # TODO: change for RNA Editing for PHIV
            tbl.write('\t\t\texception\tRibosomal slippage')
        else:
            end_int = int(product_map[k].split('#')[1])
            if end_int == len(virus_dna):
                if ((int(product_map[k].split('#')[1]) - int(product_map[k].split('#')[0])) + 1 % 3) == 0 and virus_dna[end_int - 3:end_int].upper() in 'TGA,TAA,UAG':
                    flag = ''
                else:
                    flag = '>'
            tbl.write('\n' + product_map[k].split('#')[0] + '\t' + flag + product_map[k].split('#')[1] + '\t' + type + '\n')
            tbl.write('\t\t\tproduct\t' + k.split('%')[0])
    tbl.close()

    #virus_strain = 'HPIV' #det_exact_species(virus_name, virus_dna)
    iso_source = ''
    iso_source += metadata_list[4]

    fsa = open(virus_name + '/' + virus_name + '.fsa', 'w')
    fsa.write('>' + virus_name + ' [organism=' + virus_strain + '] [collection-date=2016] [country=USA] '
              '[moltype=genomic] [host=Human] [gcode=1] [molecule=RNA] [strain=' + virus_name +']\n')
    fsa.write(virus_dna)
    fsa.close()


def det_exact_species(virus_name, virus_dna):
    temp = open('temp.fasta','w')
    temp.write('>' + virus_name + '\n')
    temp.write(virus_dna)
    temp.close()
    subprocess.call('blastn -query temp.fasta -db /Users/uwvirongs/Downloads/surpi-master/nt -num_threads 8 -max_target_seqs 1 -word_size ', shell=True)
    return input('WHAT\'S YOUR NAME MAN???')


def run_tbl(virus_name):
    # TODO: shove tbl2asn into the folder
    subprocess.call('tbl2asn -p ' + virus_name + '/ -t ' + virus_name + '/template.sbt -Y ' + virus_name +  '/assembly.cmt -V vb', shell=True)


if __name__ == '__main__':
    #subprocess.call('source ~/.bashrc', shell=True)
    virus_name_list, virus_dna_list, virus_annot_map = parse_gff(gff_file_loca)
    strain_list = {}
    for x in range(0, len(virus_name_list)):
        subprocess.call('mkdir ' + virus_name_list[x], shell=True)
    for x in range(0, len(virus_name_list)):
        temp = open(virus_name_list[x] + '.fasta', 'w')
        temp.write('>' + virus_name_list[x] + '\n')
        temp.write(virus_dna_list[x])
        temp.close()
        strain = ''

        subprocess.call('blastn -query ' + virus_name_list[x] + '.fasta -db /Users/uwvirongs/Downloads/surpi-master/nt -num_threads 8 -num_descriptions 0 -num_alignments 1 -word_size 28 | tee ' + virus_name_list[x] + '/' + virus_name_list[x] + '.result', shell=True )
        for line in open(virus_name_list[x] + '/' + virus_name_list[x] + '.result'):
            if line[0] == '>':
                # TODO: this takes the name of the virus and removes everything to the right of strain or isolate, just add more words in the same format if you get some issues
                strain = line.split('|')[4].split('strain')[0].split('isolate')[0].split(',')[0]
                #strain = ' '.join(strain)
                break

        strain_list[virus_name_list[x]] = strain
    for x in range(0, len(virus_name_list)):
        #print(x)
        #print(virus_name_list[x])
        #print(virus_annot_map.keys())
        write_a_virus(pull_metadata(virus_name_list[x]), virus_name_list[x], virus_dna_list[x], virus_annot_map[virus_name_list[x]], strain_list[virus_name_list[x]])
        run_tbl(virus_name_list[x])
