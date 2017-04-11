# This program will take a fasta file containing all the reads in a sample and map it to a reference genome
# Then hopefully I'll set this up to were I can just call all the functions and import it directly into the sorter code

import subprocess
import glob

input_dir = '/Users/uwvirongs/Documents/ClinSeq/ct-test/'

for fastq in glob.glob(input_dir + '*.fastq'):
    x = 0
    subprocess.call('centrifuge -x /Users/uwvirongs/Documents/ClinSeq/ct-test/nt/nt -q ' + fastq +
                    ' -S classification.txt --report-file report.txt -p 8', shell=True)

    # Go through the report and grab the most abundant hit for this sample (the whole line)
    # Need to rename this to the official report file that we get
    gi_to_freq = {}

    for line in open('report.txt'):
        if go:
            temp_gi = line.split()[1]
            if temp_gi in gi_to_freq.keys():
                gi_to_freq[temp_gi] += 1
            else:
                gi_to_freq[temp_gi] = 1
        else:
            go = True

    abund = 0
    for x in gi_to_freq.keys():
        if gi_to_freq[x] > abund:
            tru_gi = x
            abund = gi_to_freq[x]

    # Pull the GI from the map we loaded in memory
    # Download the reference .gbk and .fasta

    # maybe we should mkdir at this phase - but we don't even have a real FASTA yet which we NEED
    # save the results

