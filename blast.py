 #!/usr/bin python3

# --------------------------------------------------------------
# Sequence Alignment Module (2/5)
# --------------------------------------------------------------

# ------------------------- WARNINGS!---------------------------
# This program requires Python 3 to work——hence, the '.py' extension. 
# If you don't have Python 3, please crawl out of whichever rock you've 
#+been living under and install it.

# This program requires Biopython to work. If you don't have Biopython 
#+installed, what are you waiting for? A sign from God? This is it.

# This program requires BLAST to work. If you don't have BLAST installed 
#+locally, please reconsider your life choices and then download it from 
#+the NCBI's website. Since I'm nice, here's a link: 
#+https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).


# -------------------------- MODULES ---------------------------
from Bio import SeqIO
import os
import sys
import csv

# ------------------------- FUNCTIONS --------------------------

# This one does the BLASTin'! Then it filters the results according to 
#+Identity and Coverage values and generates files with the results 
#+in tsv and fasta format. It returns a dictionary with the fasta 
#+sequences and ID

def to_blast_or_not_to_blast(query_file, db_file, evalue_co, raw_results_file,
    filtered_results_file, fasta_results_file, iden_co, cov_co, query_dictionary):

    # Run BLAST 
    cmd = "blastp -query "+ query_file +" -subject "+ db_file + " -evalue "+ str(evalue_co) \
        +" -outfmt \"6 qseqid sseqid qcovs pident evalue sseq\" -out "+ raw_results_file
    os.system(cmd)

    # Open raw blast results and create a new file for the filtered results
    with open(raw_results_file,'r')as fin,\
         open (filtered_results_file,'w') as fout:
        # Print out a header
        print('Query_ID\tSubject_ID\tCov\tIden\teValue\tSequence', file=fout)

        # Set the counter
        hit_count=0

        # To make sure there are no repeated (spare) sequences, create a dictionary
        kill_the_spare={}

        # Open with csv module and set delimiter to tab
        # Parse though each row
        for row in csv.reader(fin, delimiter='\t'):

            # If the coverage and identity values are above the threshold
            if (float(row[2]) >= float(cov_co)) \
                and (float(row[3]) >= float(iden_co)) \
                and (row[1] not in kill_the_spare):

                with open(fasta_results_file, 'a') as fasta:

                    # Add sequence to blast outfile
                    print(*row, sep='\t', file=fout)

                    # Add sequence to fasta outfile
                    print(f'>{row[1]}\n{row[5]}', file=fasta)

                    # Increase counter
                    hit_count+=1

                    # Add sequence to dictionary (will be usefull later, promise)
                    kill_the_spare[row[1]+'('+row[0]+')']=row[5]

        # Add query sequence to fasta file
        with open(fasta_results_file, 'a') as fasta:
            for key, value in query_dictionary.items():
                print(f'>{key}\n{value}', file=fasta)

    # Add query sequences to all-sequences dictionary
    kill_the_spare.update(query_dictionary)

    # Print out some awesome summary of the blast
    print(f'''### RESULTS ###
    \nQuery file: {query_file} \
    \nSubject file: {db_file} \
    \nEvalue threshold: {evalue_co} \
    \nUnfiltered BLAST Results file: {raw_results_file} \
    
    \n### FILTERS ### 
    \nFiltered BLAST Results file: {filtered_results_file} \
    \nCoverage threshold: {iden_co} \
    \nIdentity threshold: {cov_co} \
    \nHits: {hit_count}''')

    if hit_count == 0:
        print('\n\033There are no hits.')
        sys.exit(ZERO)

    # Return dictionary to use later
    return kill_the_spare


