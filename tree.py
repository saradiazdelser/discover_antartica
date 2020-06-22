 #!/usr/bin python3

# --------------------------------------------------------------
# Tree-Making Module (3/5)
# --------------------------------------------------------------

# ------------------------- WARNINGS!---------------------------
# This program requires Python 3 to work——hence, the '.py' extension. If you don't have Python 3, 
#+ please crawl out of whichever rock you've been living under and install it.

# This program also requires MUSCLE to work. If you don't have MUSCLE installed, 
# you should check out how to do it here: https://www.drive5.com/muscle/downloads.htm  

# -------------------------- MODULES ---------------------------
import os


# ------------------------- FUNCTIONS --------------------------

# This one runs a multiple sequence alignment with MUSCLE.
# Tecnically, you could run it with Biopython, but it won't let you create 
#+a tree, so you'd have to do that part with the os module anyway. 

# Then it uses MUSCLE maketree and produces a netwick tree. 

def run_muscle_run(fasta_file, aligned_file, tree_file):

	# Run multiple alignment command
	cmd1 = "muscle -in " + fasta_file + " -out "+ aligned_file
	os.system(cmd1)

	# Run the tree making command
	cmd2 = "muscle -maketree -in " + aligned_file + " -out " + tree_file + " -cluster neighborjoining"
	os.system(cmd2)
	


