# This program will take a fasta file containing all the reads in a sample and map it to a reference genome
# Then hopefully I'll set this up to were I can just call all the functions and import it directly into the sorter code

import subprocess
import glob

input_dir = '/Users/uwvirongs/Documents/ClinSeq/ct-test/'

for fastq in glob.glob(input_dir + '*.fastq'):
    x = 0
    subprocess.call('centrifuge -x /Users/uwvirongs/Documents/ClinSeq/ct-test/nt/nt -q ' + fastq +
                    ' -S classification.txt --report-file report.txt -p 8', shell=True)
    abund = 0.0
    line_of_interest = ''
    # Go through the report and grab the most abundant hit for this sample (the whole line)
    for line in open('report.txt'):
        if float(line.split('\t')[6]) > abund:
            line_of_interest = line
            abund = float(line.split('\t')[6])
            ref_id = line.split('\t')[1]
            # Pull the GI from the map we loaded in memory

            # Download the reference .gbk and .fasta

    # maybe we should mkdir at this phase - but we don't even have a real FASTA yet which we NEED
    # save the results


def grab_ref(taxid, sample):

    tru_id = 'unknown'
    max = -999999999.0

    for line in open(input_dir + sample + '_report.txt'):
        id = line.split()[0]
        abundance = line.split()[6]
        if float(abundance) > max:
            tru_id = id
            max = float(abundance)
    return(tru_id)



# load the map of gi to taxid into memory -
# then when we pull the taxid we can immideatly download a .gbk file for annotations
