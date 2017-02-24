
import subprocess
import argparse
import glob
import timeit
import sys
from tools import annotate_assembled

start = timeit.default_timer()

# IMPORTANT NOTE- THIS NEEDS TO BE RUN WITH PROKKA IN -c MODE
# NOTE- You also need to have a copy of tbl2asn in your home directory
# NOTE- Needs info.csv, template.sbt, annotation.cmt
# python virus_annotate filename basecutoff runname
# TODO: Instead of pulling from the info.csv file we should pull directly from the master spreadsheet

parser = argparse.ArgumentParser(description='Annotate a set of assemblies with prokka, HHpred, and blast')
parser.add_argument('filename', help='The input file, must be a .fasta can contain as many contigs as you please')
parser.add_argument('basecutoff', type=int, help='Minimum length of contig, all assemblies of length less than this '
                                                 'will not be processed Passing a zero at this argument signfies that '
                                                 'you would like to run this on a segmented genome that you have '
                                                 'assembled into one file')
parser.add_argument('runname', help='Give your run a unique name! Useful to track the results later or if you want to '
                                    'do more than one run back to back')
args = parser.parse_args()


min_base_cutoff = args.basecutoff
filename = args.filename
run_name = args.runname

fasta_file = open(filename)
info_list = []
dna_list = []
dna_string = ''

WORKING_DIRECTORY = '/home/ryan/Virus-Annotate/'
HHBLIT_DB = '/home/ryan/uniprot20_2015_06/uniprot20_2015_06'
BLAST_DB = '/home/ryan/db/nr'

# This hacks the preassembled fasta file
if min_base_cutoff == 0:
    annotate_assembled(filename, run_name)
    sys.exit()

# This reads the assembly file and saves each contig greater than min base cuttoff
for line in fasta_file:
    if line[0] == '>':
        info_list.append(line[:-1])
        if dna_string != '':
            dna_list.append(dna_string)
            dna_string = ''
    else:
        dna_string += line[:-1]
dna_list.append(dna_string)

pruned_info_list = []
pruned_dna_list = []
count = 0
for entry in info_list:
    if int(entry.split('_')[3]) >= min_base_cutoff:
        pruned_info_list.append(entry)
        pruned_dna_list.append(dna_list[count])
        count += 1

node_list = []
subprocess.call('mkdir ' + WORKING_DIRECTORY + run_name, shell=True)
run_name = run_name + '/'

# Create a new folder for each contig saved in the last step, and save the contig in it as a .fsa with
# Info taken from info.csv (a tab seperated list of information for genbank)
for x in range(0, len(pruned_dna_list)):
    # Write the file that prokka will use
    one = open(WORKING_DIRECTORY + 'temp_contig_dir/contig_' + str(x).zfill(2) + '.fasta', 'w')
    one.write(pruned_info_list[x] + '\n')
    one.write(pruned_dna_list[x])
    one.close()

    # Write the .fsa file for tbl2asn to use, complete with the comment information
    subprocess.call('mkdir ' + WORKING_DIRECTORY + run_name + 'contig' + str(x).zfill(2), shell=True)
    two = open(WORKING_DIRECTORY + run_name + 'contig' + str(x).zfill(2) + '/contig' + str(x).zfill(2) + '.fsa', 'w')
    data = pruned_info_list[x].split('_')
    newline = data[0] + '_' + data[1] + '_' + data[2] + '_' + data[3]
    node_list.append(newline)
    # Read all the comment information in and then create new header line.
    info = open(WORKING_DIRECTORY + 'info.csv')
    count = 1
    for line in info:
        if count == 1:
            title_list = line.split(',')
        if count == 2:
            data_list = line.split(',')
        count += 1
    for index in range(0, len(title_list)):
        newline = newline + ' [' + title_list[index].rstrip('\n') + '=' + data_list[index].rstrip('\n') + ']'
    newline = newline + ' [strain=' + run_name + str(x) + ']'
    two.write(newline + '\n')
    two.write(pruned_dna_list[x])
    two.close()
    info.close()
inc = 0
# Run Prokka in viral mode on each of the putative genomes and save em' all into the prokka_output folder
for contig in sorted(glob.glob(WORKING_DIRECTORY + 'temp_contig_dir/*.fasta')):
    subprocess.call('prokka ' + contig + ' --outdir ' + WORKING_DIRECTORY + 'prokka_output/ --prefix contig' +
                    str(inc).zfill(2) + ' --kingdom Virus --force --compliant', shell=True)
    inc += 1


# move and modify the .tbl files from prokka into our final directory
count = 0
for tbl in sorted(glob.glob(WORKING_DIRECTORY + 'prokka_output/*.tbl')):
    readnext=True
    original = open(tbl)
    newfile = open(WORKING_DIRECTORY + run_name + 'contig' + str(count).zfill(2) + '/contig' + str(count).zfill(2) + '.tbl', 'w')
    for line in original:
        if line[0] == '>' or line[0] == '<':
            newfile.write('>Feature ' + node_list[count][1:] + '\n')
        elif line[0] != '\t':
            if abs(int(line.split()[1]) - int(line.split()[0])) >= 240 and line.split()[2] != 'gene':
                newfile.write(line)
                readnext = True
            else:
                readnext = False
        elif line.split()[0] == 'product' and readnext:
            newfile.write(line)
            readnext = True
    newfile.close()
    original.close()
    count += 1


# Go through all the prokka outputs line by line and push the 'genes' that have 'hypothetical protein' through hhblits
# And blastp NOTE- Blastp is running in -remote mode so this needs to get changed during final deployment
count = 0
for prokka_output in sorted(glob.glob(WORKING_DIRECTORY + 'prokka_output/*.faa')):
    dna_string = ''
    header_list = []
    gene_list = []
    for line in open(prokka_output, 'r'):
        if line[0] == '>':
            header_list.append(line[:-1])
            if dna_string != '':
                gene_list.append(dna_string)
                dna_string = ''
        else:
            dna_string += line[:-1]
    gene_list.append(dna_string)
    gene_num = 0
    for x in range(0, len(header_list)):
        # Creates a file the contig folder named gene_num.hhr only hypothetical proteins are displayed
        if len(gene_list[x]) >= 80:
            gene_num += 1
            print(header_list[x])
            print(gene_list[x])
            d = open('/home/ryan/temp_gene.part', 'w')
            d.write(header_list[x] + '\n')
            d.write(gene_list[x])
            d.close()
            # Actually call hhblits and blastp, I'm not even sure if blastp is correct here, and it is remote, but the
            # heavy lifting for this step is done already
            subprocess.call('hhblits -i /home/ryan/temp_gene.part -o ' + WORKING_DIRECTORY + run_name + 'contig' +
                            str(count).zfill(2) + '/hh_result' + str(gene_num).zfill(2) +
                           '.hhr -d ' + HHBLIT_DB + ' -Z 10 -B 10', shell=True)
            #subprocess.call('blastp -query /home/ryan/temp_gene.dna -db nr -remote | tee ' + WORKING_DIRECTORY +
                            #run_name + 'contig' + str(count).zfill(2) + '/blast_result' + str(gene_num).zfill(2) + '.txt', shell=True)
            # TODO: make this not remote whenever possible
            # TODO: also doen't like proteins not starting with MET
            # TODO: comment this out for more accurate timing results
    count += 1

count = 0
for prokka_output in sorted(glob.glob(WORKING_DIRECTORY + 'prokka_output/*.ffn')):
    dna_string = ''
    header_list = []
    gene_list = []
    for line in open(prokka_output, 'r'):
        if line[0] == '>' or line[0] == '<':
            header_list.append(line[:-1])
            if dna_string != '':
                    gene_list.append(dna_string)
            dna_string = ''
        else:
            dna_string += line[:-1]
    gene_list.append(dna_string)
    gene_num = 0
    print()

    for x in range(0, len(header_list)):
        if len(gene_list[x]) >= 240:
            gene_num += 1
            print(header_list[x])
            print(gene_list[x])
            d = open('/home/ryan/temp_gene.dna', 'w')
            d.write(header_list[x] + '\n')
            d.write(gene_list[x])
            name = 'hypothetical protein'
            d.close()
            subprocess.call('blastn -query /home/ryan/temp_gene.dna -db nt -remote | tee ' + WORKING_DIRECTORY +
                            run_name + 'contig' + str(count).zfill(2) + '/blast_result' + str(gene_num).zfill(2) + '.txt', shell=True)


    count += 1

# TODO: Pause here and allow manual annotation? or at least provide a function to go through and re-tbl2asn everything
# Cleanup step
for x in range(0, count):
    subprocess.call('cp ' + WORKING_DIRECTORY + 'template.sbt ' + WORKING_DIRECTORY + run_name + 'contig' +
                    str(x).zfill(2), shell=True)
    subprocess.call('cp ' + WORKING_DIRECTORY + 'assembly.cmt ' + WORKING_DIRECTORY + run_name + 'contig' +
                    str(x).zfill(2), shell=True)
    # This means you need to have a copy of tbl2asn in your home directory
    subprocess.call('./tbl2asn -p ' + WORKING_DIRECTORY + run_name + 'contig' + str(x).zfill(2) +
                    ' -t ' + WORKING_DIRECTORY + run_name + 'contig' + str(x).zfill(2) +
                    '/template.sbt -Y ' + WORKING_DIRECTORY + run_name + 'contig' + str(x).zfill(2) + '/assembly.cmt -V vb',
                    shell=True)

subprocess.call('rm -rf ' + WORKING_DIRECTORY + 'prokka_output/* ' + WORKING_DIRECTORY + 'temp_contig_dir/*',
                shell=True)
stop = timeit.default_timer()
print('Annotations should be in folders named contig# all unknown genes will have blastp and hhpred outputs with lower '
      'numbers corresponding to earlier unknown genes in the prokka output')
print('total run time')
print(stop - start)
