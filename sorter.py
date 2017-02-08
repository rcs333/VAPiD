# We gon make a small version of our blaster and sorter here that just takes a fasta, cause this is as far back as
# I really want to go for right now

import subprocess
# Top line is ribosomal slippage
# bottom line is for rna editing, in case I want to handle these guys differently

hard_list = ['Human coronavirus 229E', 'Human coronavirus OC43', 'Human parainfluenza virus 3', 'astrovirus 1', 'west nile', 'HIV', 'HTLV', 'influenza A',
             'HPIV', 'measels', 'mumps', 'ebola', 'Human parainfluenza virus 4', 'Human parainfluenza virus 2']

fasta_file_loc = 'test.fasta'
BLAST_DB_LOCATION = '/Users/uwvirongs/Downloads/centrifuge/nt'

subprocess.call('blastn -query ' + fasta_file_loc + ' -db ' + BLAST_DB_LOCATION + ' -num_threads 8 -num_descriptions 0 '
                '-num_alignments 1 -word_size 28 | tee test.result', shell=True)

for line in open('test.result'):
    if line[0] == '>':
        strain = line.split('|')[4]

if strain in hard_list:
    print('hard')
else:
    print('easy')

