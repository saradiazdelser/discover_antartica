
#!/usr/bin python3

# --------------------------------------------------------------
# Domain Finder Module (4/5)
# --------------------------------------------------------------

# ------------------------- WARNINGS!---------------------------
# This script requires Python 3 to work——hence, the '.py' extension. If you don't have Python 3, 
#+ please crawl out of whichever rock you've been living under and install it.

# This script also requires Biopython to work. If you don't have Biopython installed, 
#+what are you waiting for? A sign from God? This is it.

# -------------------------- MODULES ---------------------------
import re
from Bio import SeqIO
from Bio.ExPASy import Prosite,Prodoc


# ------------------------- FUNCTIONS --------------------------

# This one basically takes the fasta file given and runs it against the prosite
#+database to identify the protein's domains.
def its_a_wonderful_database(fasta_file, text_file, DOC, dictionary, prosite_dat_file, prosite_doc_file):
	
	# Open output file
	with open(text_file,'w') as output_file:

		# Import fasta dictionary
		for key,value in dictionary.items():
			print(f"### {key} ###", file=output_file)
				
			# (Re)set count
			hits=0

			# Open the prosite database 
			with open (prosite_dat_file, 'r') as prosite: 
				
				# Parse through the records
				for record in Prosite.parse(prosite):

					# Check if it exists
					if record.pattern:

						# Magic into valid regex
						pattern=one_regex_to_rule_them_all(record.pattern)

						# Search for pattern 
						if re.findall(pattern, value, re.IGNORECASE):
							
							# Count it
							hits +=1

							# Write summary of domains into text file
							print(f'\nName: {record.name}\
								\nAccession: {record.accession}\
								\nDescription: {record.description}\
								\nPattern: {record.pattern}\
								',file=output_file)

							if DOC == True:
								print(f'Documentation: \
									{whats_in_the_prodoc(id=record.accession,doc_file=prosite_doc_file)}\
									', file=output_file)

				# Let the user know how many domains were found
				print(f'{key} has {str(hits)} known domains.')


# This one turns a given prosite pattern into a python-aproved regex
# I doesn't necessarily have to be a separate function (let's face it, none of 
#+them do), but I think it's neater like this.
def one_regex_to_rule_them_all(pattern):

	pattern=pattern.replace('.', '')\
					.replace('x', '.')\
					.replace('-', '')\
					.replace('{', '[^')\
					.replace('}', ']')\
					.replace('(', '{')\
					.replace(')', '}')
	# Return the pattern, so it can be used outside this function
	return pattern


# This one generates the documentation from prosite.doc, according to the 
#+given record accession
def whats_in_the_prodoc(id, doc_file):
	# Open the prosite documentation
	# 'Encoding' is necessary because some characters are not in UTF-8
	with open(doc_file, 'r', encoding='cp1252') as documentation: 
		# Parse through it
		for record in Prodoc.parse(documentation):
			# Look for the domain accession we want
			if record.accession == id:
				# Return the documentation
				return record.text 

