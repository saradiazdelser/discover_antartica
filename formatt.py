 #!/usr/bin python3

# --------------------------------------------------------------
# Format Verification Module (1/5)
# --------------------------------------------------------------

# ------------------------- WARNINGS!---------------------------
# This program requires Python 3 to work——hence, the '.py' extension. If you don't have Python 3, 
#+ please crawl out of whichever rock you've been living under and install it.

# This program requires Biopython to work. If you don't have Biopython installed, 
#+what are you waiting for? A sign from God? This is it.

# -------------------------- MODULES ---------------------------

from Bio import SeqIO
import re

# ------------------------- FUNCTIONS --------------------------

# This one checks the format of a given file
def the_format_is_strong_in_this_one(mystery_file):
    
    with open (mystery_file, 'r') as file:
    	# Save the first line as string
        line=file.readline()
        if line.startswith('>'):
            return 'FA'

        elif line.startswith('LOCUS'):
             # The first line contains 'bp'
            if re.search('bp', line):
                return 'GB'

            # The first line contains 'aa'
            elif re.search('aa', line):
                return 'GP'
        else:
            print(f'Good Sir, {mystery_file} is not a valid format!')


# This one converts the genbank database into a multifasta
def in_the_begining_there_was_genbank(gb_file, fasta_file):

	# Open the gb file and set the fasta file as appendable  
	with open(gb_file) as input_file, open(fasta_file, "a") as output_file:

		# Read the genbank file
		for records in SeqIO.parse(input_file, "genbank"):
			# in
		    for f in records.features :
		        if (f.type=="CDS"):
		        	try:
		        		output_file.write(f">{f.qualifiers['locus_tag'][0]}_{records.name}\
		        			\n{f.qualifiers['translation'][0]}\n")
		        	except KeyError as e:
		        		continue


	print('Success! %s added to the multifasta database' % gb_file)


# This one converts the genprot database into a multifasta
def in_the_begining_there_was_genprot(gp_file, fasta_file):
	
	# Open the gb file and set the fasta file as writeable  
	with open(gp_file) as input_file, open(fasta_file, "a") as output_file:
		# Read the genbank file
		sequences_gp = SeqIO.parse(input_file, "genbank")
		# Create a fasta file
		sequences_fasta = SeqIO.write(sequences_gp, output_file, "fasta")

	print('Success! %s added to the multifasta database' % gp_file)


# This one converts the genprot database into a multifasta
def in_the_begining_there_was_fasta(fa_file, fasta_file):
	
	# Open the gb file and set the fasta file as writeable  
	with open(fa_file) as input_file, open(fasta_file, "a") as output_file:
		# Read the fasta file
		sequences_fa = SeqIO.parse(input_file, "fasta")
		# Create a fasta file
		sequences_fasta = SeqIO.write(sequences_gp, output_file, "fasta")

	print('Success! %s added to the multifasta database' % fa_file)


# This one turns a given query fasta file into a dictionary. 
def the_perks_of_being_a_dictionary(query_file):
	# Open fasta file
	with open(query_file,'r') as file:
		# Create an empty dictionary
		seq_dictionary={}
		# Parse through the fasta file
		for record in SeqIO.parse(file, 'fasta'):
			# Add entry to dictionary
			seq_dictionary['query_'+str(record.id)]=str(record.seq)

	# Return the dictionary, so it can be used outside this function	
	return seq_dictionary
